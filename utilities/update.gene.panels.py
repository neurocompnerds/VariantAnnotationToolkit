#!/usr/bin/python3

# This script fetches gene panel information from several different sources and then updates the file in ANNOVAR.

import sys, getopt, csv, os, re, io, gzip, json, shutil, subprocess, tempfile, zipfile
from datetime import datetime
import requests

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
CONFIG_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, '..', 'config'))
PANEL_LIST = os.path.join(CONFIG_DIR, 'panel.list.json')
SUPPORTED_BUILDS = ('Hs38DH', 'hg38')
PANELAPP_BUILD_KEY = 'GRch38'  # Hs38DH and hg38 both correspond to "GRch38" in PanelApp JSON
PROMOTER_PAD = 1000  # bp added upstream of each gene (strand aware)
HTTP_TIMEOUT = 120
HTTP_HEADERS = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) gene-panel-updater'}
DOWNLOAD_DIR_NAME = 'panel_downloads'


def usage():
    print(
'''
# update.gene.panels.py a script to extract gene panels from specified sources for inclusion in our ANNOVAR database.
# This utility will fetch gene panel information from Genomics England or Australian PanelApp, SFARI Gene and the Genes4Epilepsy GitHub page
# information is converted to bed format or other formats used by ANNOVAR.
# Gene coordinates are taken from the GENCODE GTF set in config/environment.cfg (or environment.local.cfg) and extended 1kb upstream.
#
# Usage update.gene.panels.py -p "ID,ASD,Epilepsy,DDG2P,CP" | -p all [-b build] [-f] | [ -l | --list ] | [ -h | --help ]
#
# Options:
# -p | --panels  A list of gene panel names to fetch. Separate multiple names with a comma. If "all" is specified, then all panels will be fetched.
# -b | --build   OPTIONAL: Genome build, either Hs38DH or hg38. Default = Hs38DH
# -f | --force   OPTIONAL: Rebuild panels even if the source has not changed since the last update.
# -l | --list    OPTIONAL: Lists all available gene panel names.
# -h | --help    OPTIONAL: Displays usage information.
#
# Outputs (written to $AV_DB for the build):
#   ${BUILD}_${panel}_latest.bed          BED panels
#   ${BUILD}_${panel}.txt                 DDG2P format panels
#   ${BUILD}_panel_updates.log            One line appended per updated file
#   ${BUILD}_${panel}_unmatched.log       Genes not matched (or matched by symbol only), overwritten each run
#   panel_downloads/                      Latest raw download for each panel
#
# Script created by Mark Corbett with AI assistance on 13/11/2026
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
# Mark Corbett with AI assistance; 18/11/2026; Implement downloads, GENCODE coordinate lookup, BED / DDG2P outputs and update logging
#
'''
         )


def log(msg):
    sys.stderr.write(msg + '\n')
    sys.stderr.flush()


# ---------------------------------------------------------------- configuration

def load_panel_list():
    with open(PANEL_LIST) as f:
        return json.load(f)


def list_panels(panels):
    print('%-10s %-45s %-28s %s' % ('name', 'long_name', 'data_source_name', 'output_format'))
    for p in panels:
        print('%-10s %-45s %-28s %s' % (p['name'], p.get('long_name', ''), p.get('data_source_name', ''), p.get('output_format', '')))


def load_environment(build):
    """Source the bash config files with BUILD set and return the resolved variables."""
    script = (
        'BUILD="$1"; source "$2/environment.cfg" >/dev/null; '
        'if [ -f "$2/environment.local.cfg" ]; then source "$2/environment.local.cfg" >/dev/null; fi; '
        'printf "%s\\n%s\\n%s\\n%s\\n" "$BUILD" "$AV_DB" "$GencodeGTF" "$GencodeURL"'
    )
    out = subprocess.run(['bash', '-c', script, 'bash', build, CONFIG_DIR],
                         check=True, capture_output=True, text=True).stdout
    values = out.rstrip('\n').split('\n')[-4:]
    env = dict(zip(('BUILD', 'AV_DB', 'GencodeGTF', 'GencodeURL'), values))
    for key, val in env.items():
        if not val:
            raise RuntimeError('Config variable %s is empty after reading %s' % (key, CONFIG_DIR))
    return env


# ---------------------------------------------------------------- file helpers

def atomic_write(path, write_func, mode='w'):
    """Write to a temp file in the destination folder then rename into place."""
    folder = os.path.dirname(path) or '.'
    fd, tmp = tempfile.mkstemp(dir=folder, prefix='.' + os.path.basename(path) + '.', suffix='.tmp')
    try:
        with os.fdopen(fd, mode) as f:
            write_func(f)
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)  # mkstemp creates files readable by owner only
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def http_get(url, **kwargs):
    r = requests.get(url, headers=kwargs.pop('headers', HTTP_HEADERS), timeout=HTTP_TIMEOUT, **kwargs)
    r.raise_for_status()
    return r


def ordinal_date(dt):
    day = dt.day
    suffix = 'th' if 11 <= day % 100 <= 13 else {1: 'st', 2: 'nd', 3: 'rd'}.get(day % 10, 'th')
    return '%d%s %s at %s' % (day, suffix, dt.strftime('%B %Y'), dt.strftime('%H:%M'))


def clean_text(value):
    """Make a value safe for ANNOVAR outputs that are later combined as CSV."""
    if value is None:
        return '.'
    if isinstance(value, (list, tuple)):
        value = '; '.join(str(v) for v in value if v is not None and str(v).strip() != '')
    value = str(value).replace(',', ';').replace('"', '')
    value = re.sub(r'\s+', ' ', value).strip()
    value = re.sub(r'\s*;\s*', '; ', value).strip('; ').strip()
    return value if value else '.'


def strip_version(accession):
    """ENSG00000156531.17 -> ENSG00000156531; ENSG00000182378.15_PAR_Y -> ENSG00000182378"""
    m = re.match(r'^(ENS[A-Z]*[GT]\d+)', accession.strip())
    return m.group(1) if m else accession.strip().split('.')[0]


def split_accessions(value):
    if not value:
        return []
    return [strip_version(a) for a in re.split(r'[;,|\s]+', str(value)) if a.strip()]


# ---------------------------------------------------------------- GENCODE lookup

def ensure_gtf(env):
    path = env['GencodeGTF']
    if os.path.exists(path):
        return path
    log('GENCODE GTF not found at %s. Downloading from %s' % (path, env['GencodeURL']))
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with requests.get(env['GencodeURL'], headers=HTTP_HEADERS, timeout=HTTP_TIMEOUT, stream=True) as r:
        r.raise_for_status()
        atomic_write(path, lambda f: shutil.copyfileobj(r.raw, f, 1024 * 1024), mode='wb')
    log('GENCODE GTF saved to %s' % path)
    return path


class GeneLookup:
    """Gene coordinates from a GENCODE GTF keyed by versionless Ensembl gene ID."""

    ATTR_RE = re.compile(r'(\S+) "([^"]*)"')

    def __init__(self, gtf_path):
        self.genes = {}        # gene_id -> list of (chrom, start0, end, strand, symbol)
        self.transcripts = {}  # transcript_id -> gene_id
        self.symbols = {}      # symbol -> set(gene_id)
        opener = gzip.open if gtf_path.endswith('.gz') else open
        log('Reading gene coordinates from %s' % gtf_path)
        with opener(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t', 8)
                if len(cols) < 9 or cols[2] not in ('gene', 'transcript'):
                    continue
                attrs = dict(self.ATTR_RE.findall(cols[8]))
                gene_id = strip_version(attrs.get('gene_id', ''))
                if not gene_id:
                    continue
                if cols[2] == 'transcript':
                    if 'transcript_id' in attrs:
                        self.transcripts[strip_version(attrs['transcript_id'])] = gene_id
                    continue
                symbol = attrs.get('gene_name', '')
                # GTF is 1-based inclusive; store 0-based half-open
                self.genes.setdefault(gene_id, []).append((cols[0], int(cols[3]) - 1, int(cols[4]), cols[6], symbol))
                if symbol:
                    self.symbols.setdefault(symbol, set()).add(gene_id)
        log('Loaded %d genes and %d transcripts' % (len(self.genes), len(self.transcripts)))

    def find(self, accessions, symbol):
        """Return (gene_id, records, method) where method is 'ensembl', 'symbol' or None."""
        for acc in accessions:
            gene_id = acc if acc in self.genes else self.transcripts.get(acc)
            if gene_id in self.genes:
                return gene_id, self.genes[gene_id], 'ensembl'
        if symbol:
            ids = self.symbols.get(symbol.strip(), set())
            if len(ids) == 1:
                gene_id = next(iter(ids))
                return gene_id, self.genes[gene_id], 'symbol'
        return None, [], None


def padded_region(record):
    """Extend the gene 1kb upstream relative to strand. Coordinates are 0-based half-open."""
    chrom, start, end, strand, _ = record
    if strand == '-':
        end += PROMOTER_PAD
    else:
        start = max(0, start - PROMOTER_PAD)
    return chrom, start, end


def chrom_sort_key(chrom):
    c = chrom[3:] if chrom.lower().startswith('chr') else chrom
    if c.isdigit():
        return (0, int(c), '')
    order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
    if c.upper() in order:
        return (0, order[c.upper()], '')
    return (1, 0, chrom)


# ---------------------------------------------------------------- panel building

class PanelResult:
    def __init__(self):
        self.rows = set()       # tuples of output fields; first three after any prefix are coordinates
        self.genes_written = set()
        self.problems = []      # (status, symbol, accessions, note)


def add_gene(result, lookup, accessions, symbol, make_fields):
    gene_id, records, method = lookup.find(accessions, symbol)
    supplied = ','.join(accessions) if accessions else '.'
    if method is None:
        result.problems.append(('UNMATCHED', symbol or '.', supplied, 'No GENCODE gene found for accession or unique symbol'))
        return
    if method == 'symbol':
        result.problems.append(('SYMBOL_MATCH', symbol, supplied, 'Matched by gene symbol to %s' % gene_id))
    for rec in records:
        result.rows.add(make_fields(padded_region(rec)))
    result.genes_written.add(gene_id)


def panelapp_accessions(gene):
    builds = gene.get('gene_data', {}).get('ensembl_genes') or {}
    build_data = next((v for k, v in builds.items() if k.lower() == PANELAPP_BUILD_KEY.lower()), {}) or {}
    accs = []
    for release in sorted(build_data, key=lambda r: int(r) if str(r).isdigit() else -1, reverse=True):
        accs += split_accessions(build_data[release].get('ensembl_id', ''))
    return accs


def build_panelapp(data, panel, lookup):
    result = PanelResult()
    output_format = panel['output_format'].upper()
    for gene in data.get('genes', []):
        gd = gene.get('gene_data', {}) or {}
        symbol = gd.get('gene_symbol') or gene.get('entity_name') or ''
        level = clean_text(gene.get('confidence_level'))
        if output_format == 'DDG2P':
            omim = gd.get('omim_gene') or []
            extra = (clean_text(symbol), clean_text(gene.get('mode_of_inheritance')),
                     clean_text(gene.get('mode_of_pathogenicity')), level,
                     clean_text(gene.get('phenotypes')), clean_text(omim[0] if omim else None))
            make = lambda region, extra=extra: region + extra
        else:
            name = clean_text('%s_Level_%s' % (symbol, level))
            make = lambda region, name=name: region + (name,)
        add_gene(result, lookup, panelapp_accessions(gene), symbol, make)
    return result


def read_table(text, delimiter):
    reader = csv.reader(io.StringIO(text), delimiter=delimiter)
    rows = [r for r in reader if r and any(c.strip() for c in r)]
    return rows[0], rows[1:]


def column_index(header, name, fallback, panel_name):
    lowered = [h.strip().lower() for h in header]
    if name.lower() in lowered:
        return lowered.index(name.lower())
    log('WARNING: %s: column "%s" not found in header. Using field %d' % (panel_name, name, fallback + 1))
    return fallback


def build_table(text, delimiter, panel, lookup, columns, make_name):
    """Generic builder for Genes4Epilepsy and SFARI tables. columns maps role -> (header name, fallback index)."""
    result = PanelResult()
    header, rows = read_table(text, delimiter)
    idx = {role: column_index(header, name, fb, panel['name']) for role, (name, fb) in columns.items()}
    get = lambda row, role: row[idx[role]] if idx[role] < len(row) else ''
    for row in rows:
        name = make_name(lambda role: get(row, role))
        add_gene(result, lookup, split_accessions(get(row, 'accession')), get(row, 'symbol').strip(),
                 lambda region, name=name: region + (name,))
    return result


def build_epilepsy(text, panel, lookup):
    cols = {'symbol': ('Gene', 1), 'accession': ('Ensemble_ID', 2), 'inheritance': ('Inheritance', 5), 'phenotype': ('Phenotype(s)', 6)}
    return build_table(text, '\t', panel, lookup, cols,
                       lambda g: '%s_%s' % (clean_text(g('inheritance')), clean_text(g('phenotype'))))


def build_sfari(text, panel, lookup):
    cols = {'symbol': ('gene-symbol', 1), 'accession': ('ensembl-id', 3), 'category': ('genetic-category', 5), 'score': ('gene-score', 6)}
    return build_table(text, ',', panel, lookup, cols,
                       lambda g: '%s; gene-score:%s' % (clean_text(g('category')), clean_text(g('score'))))


# ---------------------------------------------------------------- source fetching

class Source:
    """What was fetched for a panel and how to store it."""
    def __init__(self, content, store_name, description, unchanged):
        self.content = content          # bytes
        self.store_name = store_name    # filename inside panel_downloads
        self.description = description  # text for the update log
        self.unchanged = unchanged      # True if identical to the last stored download


def stored_files(download_dir, panel_name):
    prefix = panel_name + '_'
    return sorted(f for f in os.listdir(download_dir) if f.startswith(prefix) and not f.endswith('.tmp'))


def fetch_panelapp(panel, download_dir):
    url = panel['URL'].rstrip('/') + '/'
    content = http_get(url, params={'format': 'json'}).content
    version = str(json.loads(content).get('version', 'unknown'))
    store = '%s_PanelApp_v%s.json' % (panel['name'], version)
    desc = '%s panel %s version %s' % (panel['data_source_name'], url.rstrip('/').split('/')[-1], version)
    return Source(content, store, desc, store in stored_files(download_dir, panel['name']))


def fetch_genes4epilepsy(panel, download_dir):
    release = http_get(panel['URL'], headers=dict(HTTP_HEADERS, Accept='application/vnd.github+json')).json()
    # The release asset may be the plain TSV or a zip containing it. Prefer the plain TSV when both exist.
    assets = [a for a in release.get('assets', []) if re.match(r'^EpilepsyGenes_v\d{4}-\d{2}\.tsv(\.zip)?$', a.get('name', ''))]
    if not assets:
        raise RuntimeError('No EpilepsyGenes_vYYYY-MM.tsv or .tsv.zip asset in the latest Genes4Epilepsy release')
    asset = max(assets, key=lambda a: (a['name'][:len('EpilepsyGenes_v0000-00')], a['name'].endswith('.tsv')))
    new_version = re.search(r'v(\d{4}-\d{2})', asset['name']).group(1)
    store = '%s_%s' % (panel['name'], asset['name'])
    old_versions = [m.group(1) for f in stored_files(download_dir, panel['name'])
                    for m in [re.search(r'_EpilepsyGenes_v(\d{4}-\d{2})\.tsv(\.zip)?$', f)] if m]
    desc = '%s %s' % (panel['data_source_name'], asset['name'])
    if old_versions and max(old_versions) >= new_version:
        return Source(None, store, desc, True)
    content = http_get(asset['browser_download_url']).content
    return Source(content, store, desc, False)


def epilepsy_tsv_text(src):
    """Return the TSV text from a Genes4Epilepsy download, unzipping it if needed."""
    content = src.content
    if src.store_name.endswith('.zip'):
        with zipfile.ZipFile(io.BytesIO(content)) as z:
            members = [n for n in z.namelist() if n.endswith('.tsv') and not os.path.basename(n).startswith('.')]
            if len(members) != 1:
                raise RuntimeError('Expected one .tsv file in %s but found %d' % (src.store_name, len(members)))
            content = z.read(members[0])
    return content.decode('utf-8-sig', errors='replace')


def fetch_sfari(panel, download_dir):
    content = http_get(panel['URL']).content
    if content.lstrip()[:1] in (b'<', b'{'):
        raise RuntimeError('SFARI download does not look like a CSV file (starts with %r)' % content.lstrip()[:20])
    unchanged = False
    for f in stored_files(download_dir, panel['name']):
        with open(os.path.join(download_dir, f), 'rb') as fh:
            if fh.read() == content:
                unchanged = True
    store = '%s_SFARI_%s.csv' % (panel['name'], datetime.now().strftime('%Y-%m-%d'))
    desc = '%s download of %s' % (panel['data_source_name'], datetime.now().strftime('%Y-%m-%d'))
    return Source(content, store, desc, unchanged)


FETCHERS = {
    'Genomics_England_PanelApp': fetch_panelapp,
    'PanelApp_Australia': fetch_panelapp,
    'Genes4Epilepsy': fetch_genes4epilepsy,
    'SFARI_GENE': fetch_sfari,
}


# ---------------------------------------------------------------- per panel processing

def output_path(env, panel):
    if panel['output_format'].upper() == 'DDG2P':
        name = '%s_%s.txt' % (env['BUILD'], panel['name'])
    else:
        name = '%s_%s_latest.bed' % (env['BUILD'], panel['name'])
    return os.path.join(env['AV_DB'], name)


def process_panel(panel, env, lookup_loader, force):
    source_name = panel.get('data_source_name')
    if source_name not in FETCHERS:
        raise RuntimeError('Unknown data_source_name "%s"' % source_name)
    if panel.get('output_format', '').upper() not in ('BED', 'DDG2P'):
        raise RuntimeError('Unknown output_format "%s"' % panel.get('output_format'))
    out_path = output_path(env, panel)
    download_dir = os.path.join(env['AV_DB'], DOWNLOAD_DIR_NAME)
    os.makedirs(download_dir, exist_ok=True)

    log('== %s: fetching from %s' % (panel['name'], source_name))
    src = FETCHERS[source_name](panel, download_dir)
    if src.unchanged and os.path.exists(out_path) and not force:
        log('%s: source unchanged (%s). %s not updated.' % (panel['name'], src.description, os.path.basename(out_path)))
        return
    if src.content is None:  # Genes4Epilepsy unchanged but rebuild needed: reuse stored copy
        with open(os.path.join(download_dir, src.store_name), 'rb') as fh:
            src.content = fh.read()

    lookup = lookup_loader()
    if source_name in ('Genomics_England_PanelApp', 'PanelApp_Australia'):
        result = build_panelapp(json.loads(src.content), panel, lookup)
    elif source_name == 'Genes4Epilepsy':
        result = build_epilepsy(epilepsy_tsv_text(src), panel, lookup)
    else:
        result = build_sfari(src.content.decode('utf-8-sig', errors='replace'), panel, lookup)
    if not result.rows:
        raise RuntimeError('No genes could be placed on the genome. Output not written.')

    is_ddg2p = panel['output_format'].upper() == 'DDG2P'
    rows = sorted(result.rows, key=lambda r: (chrom_sort_key(r[0]), r[1], r[2]) + tuple(r[3:]))
    def write_rows(f):
        for r in rows:
            fields = (('1',) if is_ddg2p else ()) + (r[0], str(r[1]), str(r[2])) + tuple(r[3:])
            f.write('\t'.join(fields) + '\n')
    atomic_write(out_path, write_rows)

    # Keep the raw download, then remove older downloads for this panel
    store_path = os.path.join(download_dir, src.store_name)
    if not os.path.exists(store_path) or not src.unchanged:
        atomic_write(store_path, lambda f: f.write(src.content), mode='wb')
    for old in stored_files(download_dir, panel['name']):
        if old != src.store_name:
            os.remove(os.path.join(download_dir, old))

    now = datetime.now()
    unmatched = sum(1 for p in result.problems if p[0] == 'UNMATCHED')
    def write_problems(f):
        f.write('# %s gene matching report generated on %s. Source: %s\n' % (panel['name'], ordinal_date(now), src.description))
        f.write('status\tgene_symbol\tsupplied_accession\tnote\n')
        for p in result.problems:
            f.write('\t'.join(p) + '\n')
    atomic_write(os.path.join(env['AV_DB'], '%s_%s_unmatched.log' % (env['BUILD'], panel['name'])), write_problems)
    with open(os.path.join(env['AV_DB'], '%s_panel_updates.log' % env['BUILD']), 'a') as f:
        f.write('%s updated on %s. Source: %s. Genes written: %d. Unmatched: %d.\n'
                % (os.path.basename(out_path), ordinal_date(now), src.description, len(result.genes_written), unmatched))
    log('%s: wrote %s (%d genes, %d rows, %d unmatched, %d symbol matches)'
        % (panel['name'], out_path, len(result.genes_written), len(rows), unmatched, len(result.problems) - unmatched))


# ---------------------------------------------------------------- main

def main(argv):
    try:
        opts, args = getopt.getopt(argv, 'hlp:b:f', ['help', 'list', 'panels=', 'build=', 'force'])
    except getopt.GetoptError as e:
        log('ERROR: %s' % e)
        usage()
        return 2
    if not opts:
        usage()
        return 0
    panels = load_panel_list()
    requested, build, force = None, 'Hs38DH', False
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            return 0
        if opt in ('-l', '--list'):
            list_panels(panels)
            return 0
        if opt in ('-p', '--panels'):
            requested = [p.strip() for p in arg.split(',') if p.strip()]
        elif opt in ('-b', '--build'):
            build = arg
        elif opt in ('-f', '--force'):
            force = True
    if not requested:
        log('ERROR: No panels requested. Use -p')
        usage()
        return 2
    if build not in SUPPORTED_BUILDS:
        log('ERROR: Build "%s" is not supported. Choose one of: %s' % (build, ', '.join(SUPPORTED_BUILDS)))
        return 2

    by_name = {p['name']: p for p in panels}
    if requested == ['all']:
        selected = panels
    else:
        unknown = [r for r in requested if r not in by_name]
        if unknown:
            log('ERROR: Unknown panel name(s): %s. Valid names: %s' % (', '.join(unknown), ', '.join(by_name)))
            return 2
        selected = [by_name[r] for r in requested]

    env = load_environment(build)
    log('Build %s. ANNOVAR database folder %s' % (env['BUILD'], env['AV_DB']))
    if not os.path.isdir(env['AV_DB']):
        log('ERROR: AV_DB folder %s does not exist' % env['AV_DB'])
        return 1

    cache = {}
    def lookup_loader():  # Only read the GTF if at least one panel needs rebuilding
        if 'lookup' not in cache:
            cache['lookup'] = GeneLookup(ensure_gtf(env))
        return cache['lookup']

    failed = []
    for panel in selected:
        try:
            process_panel(panel, env, lookup_loader, force)
        except Exception as e:
            log('ERROR: %s failed: %s' % (panel['name'], e))
            failed.append(panel['name'])
    if failed:
        log('Finished with failures: %s' % ', '.join(failed))
        return 1
    log('Finished')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
