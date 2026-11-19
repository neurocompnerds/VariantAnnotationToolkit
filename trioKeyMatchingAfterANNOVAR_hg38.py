#!/usr/bin/python3

# Script to filter trios for rare possibly disease causing alleles in the child, covers IBD, comp het, X-linked, autosomal dominant and clinVar flagged genotypes from a multisample ANNOVAR file
import csv
import getopt
import sys

import pandas as pd


def usage():
    print(
'''
# trioKeyMatchingAfterANNOVAR_hg38.py Script to filter trios for rare possibly disease causing alleles in the child, 
# covers IBD, comp het, X-linked, autosomal dominant and clinVar flagged genotypes from a multisample ANNOVAR file
# outputs various filtered tables for further analysis in excel.
#
# Usage: trioKeyMatchingAfterANNOVAR_hg38.py -i ANNOVAR.table.txt -c child_ID -m mother_ID -f father_ID | [ -h | --help ]
#
# Options:
# -i           /path/to/inputFile    REQUIRED: A multisample ANNOVAR table in tab delimited format
# -c           child_ID              REQUIRED: The ID of the affected child's sample as listed in the ANNOVAR table
# -m           mother_ID             REQUIRED: The ID of the mother's sample as listed in the ANNOVAR table
# -f           father_ID             REQUIRED: The ID of the father's sample as listed in the ANNOVAR table
# -h | --help  Displays help         OPTIONAL: Displays usage information.
#
# Script created by Mark Corbett on 20/12/2019
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Date; Name; Description)
# 08/12/2021; Mark; Add gnomADv3 geneotypes AF column to the 0.0001 filter list. Fix Gene.refGene. Change best gene candidate filter to whitelist.
# Mark Corbett; 06/12/2023; Add in phased genotypes
#
'''
         )

# Set initial values
inputFile = ''
geneTerms = ['exonic', 'splicing', 'UTR5', 'ncRNA_exonic', 'ncRNA_splicing']
notGeneTerms = ['downstream', 'intergenic', 'intronic', 'ncRNA_exonic', 'ncRNA_intronic', 'ncRNA_splicing', 'ncRNA_UTR3', 'ncRNA_UTR5', 'upstream', 'UTR3', 'UTR5']
filterTerms = ['.', 'PASS']
ncSpliceTerms = ['splicing', 'intronic']
filter005 = ['esp6500siv2_all', '1000g2015aug_all']
filter0001 = ['exac03', 'gnomad211_exome', 'gnomad312_genome', 'AF']
pathogenicFilter = ['Pathogenic', 'Likely_pathogenic']
nullAlelles = ['0/0', r'0\|0', r'\./\.']

# Read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],'hi:c:m:f:',['help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-i"):
        inputFile = arg
    elif opt in ("-c"):
        childID = arg
    elif opt in ("-m"):
        mumID = arg
    elif opt in ("-f"):
        dadID = arg

# Make sure you have what you need
if inputFile == '':
    usage()
    print('Hey, you forgot to tell me which ANNOVAR file to filter\n')
    sys.exit(2)

# Create the filter function
def bestGeneCandidatesFilter(df):
    df=df[df['FILTER'].isin(filterTerms)]
    df=df[(df[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
    df=df[(df[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]
    df=df[df[f"Func.{geneAnnotationSource}"].isin(geneTerms)]
    return df

# Count the number of columns in the ANNOVAR table
with open(inputFile) as f:
    reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
    first_row = next(reader)
    num_cols = len(first_row)-1 # use -1 to set the numbering to 0 based

# Open ANNOVAR table with pandas setting the chr-start-ref-obs column as the index
ANNOVARtable=pd.read_csv(inputFile, sep='\t', index_col = num_cols)

# Find out which annotation source was used
if 'Func.refGene' in ANNOVARtable.columns:
    geneAnnotationSource = 'refGene'
elif 'Func.ensGene' in ANNOVARtable.columns:
    geneAnnotationSource = 'ensGene'
elif 'Func.knownGene' in ANNOVARtable.columns:
    geneAnnotationSource = 'knownGene'
elif 'Func.gene' in ANNOVARtable.columns:
    geneAnnotationSource = 'gene'
else:
    print('ERROR: No valid gene annotation source found.')
    sys.exit(1)
samples = [mumID, dadID, childID]

# de novo
dnList=ANNOVARtable[ANNOVARtable[samples[0]].str.contains('|'.join(nullAlelles)) & ANNOVARtable[samples[1]].str.contains('|'.join(nullAlelles)) & ANNOVARtable[samples[2]].str.match(pat = r'(0/1)|(0\|1)|(1\|0)')]
dnList.to_csv(childID+".dn."+inputFile, sep='\t')
spliceCandidates=dnList[dnList[f"Func.{geneAnnotationSource}"].isin(ncSpliceTerms)]
spliceCandidates.to_csv(childID+".dn.SpliceCandidates."+inputFile, sep='\t')
dnList=bestGeneCandidatesFilter(df=dnList)
dnList.to_csv(childID+".dn.BestGeneCandidates."+inputFile, sep='\t')

# AR, identical by descent and X-linked 
homList=ANNOVARtable[~ANNOVARtable[samples[0]].str.match(pat = r'(1/1)|(1\|1)') & ~ANNOVARtable[samples[1]].str.match(pat = r'(1/1)|(1\|1)') & ANNOVARtable[samples[2]].str.match(pat = r'(1/1)|(1\|1)')]
homList.to_csv(childID+".ibdAndXl."+inputFile, sep='\t')
spliceCandidates=homList[homList[f"Func.{geneAnnotationSource}"].isin(ncSpliceTerms)]
spliceCandidates.to_csv(childID+".ibdAndXl.SpliceCandidates."+inputFile, sep='\t')
homList=bestGeneCandidatesFilter(df=homList)
homList.to_csv(childID+".ibdAndXl.BestGeneCandidates."+inputFile, sep='\t')

# Compound het calls
mNotfHets=ANNOVARtable[ANNOVARtable[samples[0]].str.match(pat = r'(0/1)|(0\|1)|(1\|0)') & ANNOVARtable[samples[1]].str.contains('|'.join(nullAlelles)) & ANNOVARtable[samples[2]].str.match(pat = r'(0/1)|(0\|1)|(1\|0)')]
mGenes=pd.unique(mNotfHets[f"Gene.{geneAnnotationSource}"])
fNotmHets=ANNOVARtable[ANNOVARtable[samples[0]].str.contains('|'.join(nullAlelles)) & ANNOVARtable[samples[1]].str.match(pat = r'(0/1)|(0\|1)|(1\|0)') & ANNOVARtable[samples[2]].str.match(pat = r'(0/1)|(0\|1)|(1\|0)')]
fGenes=pd.unique(fNotmHets[f"Gene.{geneAnnotationSource}"])
seriesCHgenes=pd.Series(mGenes.tolist() + fGenes.tolist())
chGenes=seriesCHgenes[seriesCHgenes.duplicated()]
compHets=pd.concat([mNotfHets, fNotmHets], axis=0, join='outer')
compHets=compHets[compHets[f"Gene.{geneAnnotationSource}"].isin(chGenes)] # All possible compHets
# Independently apply filters to mum and dad lists then filter the CH list
filtmNotfHets=bestGeneCandidatesFilter(df=mNotfHets)
filtfNotmHets=bestGeneCandidatesFilter(df=fNotmHets)
mGenes=pd.unique(filtmNotfHets[f"Gene.{geneAnnotationSource}"])
fGenes=pd.unique(filtfNotmHets[f"Gene.{geneAnnotationSource}"])
seriesCHgenes=pd.Series(mGenes.tolist() + fGenes.tolist())
chGenes=seriesCHgenes[seriesCHgenes.duplicated()]
compHets=compHets[compHets[f"Gene.{geneAnnotationSource}"].isin(chGenes)]
compHets.to_csv(childID+".ch."+inputFile, sep='\t')
spliceCandidates=compHets[compHets[f"Func.{geneAnnotationSource}"].isin(ncSpliceTerms)]
spliceCandidates.to_csv(childID+".ch.SpliceCandidates."+inputFile, sep='\t')
compHets=bestGeneCandidatesFilter(df=compHets)
compHets=compHets[compHets[f"Gene.{geneAnnotationSource}"].duplicated(keep=False)]  # Re-run the gene filter after the other filters
compHets.to_csv(childID+".ch.BestGeneCandidates."+inputFile, sep='\t')

# AD
hetList=ANNOVARtable[ANNOVARtable[samples[2]].str.match(pat = r'(0/1)|(0\|1)|(1\|0)')]
#hetList.to_csv("allHets."+inputFile, sep='\t') #Not likely to be worth writing out
spliceCandidates=hetList[hetList[f"Func.{geneAnnotationSource}"].isin(ncSpliceTerms)]
spliceCandidates.to_csv(childID+".allHets.SpliceCandidates."+inputFile, sep='\t')
hetList=bestGeneCandidatesFilter(df=hetList)
hetList.to_csv(childID+".allHets.BestGeneCandidates."+inputFile, sep='\t')

# ClinVar
cvList=ANNOVARtable[~ANNOVARtable[samples[2]].str.contains('|'.join(nullAlelles)) & ANNOVARtable['CLNSIG'].str.contains('|'.join(pathogenicFilter))]
cvList.to_csv(childID+".clinVar."+inputFile, sep='\t')
