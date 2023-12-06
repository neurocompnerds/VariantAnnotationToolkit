#!/usr/bin/python3

# Script to filter affected family members for matched genotypes in a multisample ANNOVAR file
import pandas as pd
import sys, getopt, csv

def usage():
    print(
'''
# familyKeyMatchingAfterANNOVAR.py a script to filter affected family members for matched genotypes in a multisample ANNOVAR file
# for rare possibly disease causing alleles. Covers IBD, comp het, X-linked, autosomal dominant and clinVar flagged genotypes
# outputs various filtered tables for further analysis in excel.
#
# Usage familyKeyMatchingAfterANNOVAR.py -i ANNOVAR.table.txt -s sampleList.txt | [ -h | --help ]
#
# Options:
# -i           /path/to/inputFile    REQUIRED: A multisample ANNOVAR table in tab delimited format
# -s           /path/to/sampleFile   OPTIONAL: A list of specific samples to extract. By default all samples are assumed affected and this might not be what you want
# -h | --help  Displays help                 OPTIONAL: Displays usage information.
#
# Script created by Mark Corbett on 15/08/2019
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
# Ali Gardner; 21/01/2021; Tweak to use with hg38 (no UK10, Wellderley, Exac.ro.1.filtered), change Func.gene to Func.refGene
# Thomas Litster; 27/04.2022; Added clinvar search (Will search for clinvar variants in all samples provided, even if variant is not shared)
# Mark Corbett; 06/12/2023; Add in phased genotypes
'''
         )

# Set initial values
inputFile = ''
sampleFile = ''
geneTerms = ['exonic', 'splicing', 'UTR5', 'ncRNA_exonic', 'ncRNA_splicing']
notGeneTerms = ['downstream', 'intergenic', 'intronic', 'ncRNA_exonic', 'ncRNA_intronic', 'ncRNA_splicing', 'ncRNA_UTR3', 'ncRNA_UTR5', 'upstream', 'UTR3', 'UTR5']
filterTerms = ['.', 'PASS']
ncSpliceTerms = ['splicing', 'intronic']
filter005 = ['esp6500siv2_all', '1000g2015aug_all']
filter0001 = ['exac03', 'gnomad211_exome', 'gnomad211_genome', 'AF']
pathogenicFilter = ['Pathogenic', 'Likely_pathogenic']
nullAlelles = ['0/0', '0|0', '\./\.']
# Read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],'hi:s:',['help'])
except getopt.GetoptError:
    usage
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-i"):
        inputFile = arg
    elif opt in ("-s"):
        sampleFile = arg

# Make sure you have what you need
if inputFile == '':
    usage()
    print('Hey, you forgot to tell me which ANNOVAR file to filter\n')
    sys.exit(2)    

# Count the number of columns in the ANNOVAR table
with open(inputFile) as f:
    reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
    first_row = next(reader)
    num_cols = len(first_row)-1 # use -1 to set the numbering to 0 based

# Open ANNOVAR table with pandas setting the chr-start-ref-obs column as the index
ANNOVARtable=pd.read_csv(inputFile, sep='\t', index_col = num_cols) 
coreTable=ANNOVARtable.loc[:,:'FORMAT']

if sampleFile =='':
    print('INFO: Using all samples available in the file\n')
    samples = ANNOVARtable.columns[ANNOVARtable.columns.get_loc("FORMAT")+1:]
elif sampleFile !='':
    samples = [line.rstrip() for line in open(sampleFile)]

dfCore=coreTable
for s in samples: # Maybe this loop could be an apply function?
    currentSampleList=ANNOVARtable[[s]]
    homList=currentSampleList[currentSampleList[s].str.match(pat = '(1/1)|(1\|1)')]
    dfCore = pd.concat([dfCore,homList], axis=1, join='inner') # Add , sort='False' once Ubuntu is upgraded

dfCore.to_csv("ibdAndXl."+inputFile, sep='\t')

#Generic filters for most likely pathogenic
dfCore=dfCore[dfCore['FILTER'].isin(filterTerms)]
dfCore=dfCore[(dfCore[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
dfCore=dfCore[(dfCore[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]

#BestGeneCandidates
bgc=dfCore[dfCore['Func.refGene'].isin(geneTerms)]
bgc.to_csv("ibdAndXl.BestGeneCandidates."+inputFile, sep='\t')

# Cadidates to test with spliceAI
spliceCandidates=dfCore[dfCore['Func.refGene'].isin(ncSpliceTerms)]
spliceCandidates.to_csv("ibdAndXl.SpliceCandidates."+inputFile, sep='\t')

# Reset and repeat for het calls
dfCore=coreTable 
for s in samples:  
    currentSampleList=ANNOVARtable[[s]]
    homList=currentSampleList[currentSampleList[s].str.match(pat = '(0/1)|(0\|1)|(1\|0)')]
    dfCore = pd.concat([dfCore,homList], axis=1, join='inner') # Add , sort='False' once Ubuntu is upgraded

dfCore.to_csv("het."+inputFile, sep='\t')
dfCore=dfCore[dfCore['FILTER'].isin(filterTerms)]
dfCore=dfCore[(dfCore[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
dfCore=dfCore[(dfCore[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]

#BestGeneCandidates
bgc=dfCore[dfCore['Func.refGene'].isin(geneTerms)]
bgc.to_csv("het.BestGeneCandidates."+inputFile, sep='\t')

# Find cadidates to test with spliceAI
spliceCandidates=dfCore[dfCore['Func.refGene'].isin(ncSpliceTerms)]
spliceCandidates.to_csv("het.SpliceCandidates."+inputFile, sep='\t')

# Find any ClinVar variants
SampleStr=''
CLNSIGStr=''
for s in samples:
    CLNSIGStr = CLNSIGStr + "~ANNOVARtable['" + s + "'].str.contains('|'.join(nullAlelles)) | "
    SampleStr = SampleStr + s + "_"

cvList=ANNOVARtable[eval(CLNSIGStr[:-3]) & ANNOVARtable['CLNSIG'].str.contains('|'.join(pathogenicFilter))]
cvList.to_csv(SampleStr+"clinVar."+inputFile, sep='\t')
