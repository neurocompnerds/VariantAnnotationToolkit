#!/usr/bin/python3

# Script to filter parents for recessive and X-linked genotypes from in a multisample ANNOVAR file
import pandas as pd
import sys, getopt, csv

def usage():
    print(
'''
# preConceptionTesting.py a script to filter affected family members for matched genotypes in a multisample ANNOVAR file
# also outputting a BestGeneCandidates file.  
#
# Usage preConceptionTesting.py -i ANNOVAR.table.txt -s sampleList.txt | [ -h | --help ]
#
# Options:
# -i           /path/to/inputFile    REQUIRED: A multisample ANNOVAR table in tab delimited format
# -m           mother_ID   REQUIRED: The ID of the mother's sample as listed in the ANNOVAR table
# -f           father_ID   REQUIRED: The ID of the father's sample as listed in the ANNOVAR table
# -h | --help  Displays help                 OPTIONAL: Displays usage information.
#
# Script created by Mark Corbett on 15/08/2019
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
#
'''
         )

# Set initial values
inputFile = ''
sampleFile = ''
notGeneTerms = ['downstream', 'intergenic', 'intronic', 'ncRNA_exonic', 'ncRNA_intronic', 'ncRNA_splicing', 'ncRNA_UTR3', 'ncRNA_UTR5', 'upstream', 'UTR3', 'UTR5']
filterTerms = ['.', 'PASS']
ncSpliceTerms = ['splicing', 'intronic']
filter005 = ['esp6500siv2_all', '1000g2015aug_all', 'UK10K-AF-all']
filter0001 = ['ExAC.r0.1.filtered', 'exac03', 'gnomad211_exome', 'gnomad211_genome']
nullAlelles = ['0/0', '\./\.']
# Read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],'hi:m:d:',['help'])
except getopt.GetoptError:
    usage
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-i"):
        inputFile = arg
    elif opt in ("-m"):
        mumID = arg
    elif opt in ("-f"):
        dadID = arg

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
samples = [mumID, dadID]

hetList=ANNOVARtable[ANNOVARtable[samples[0]].str.match('0/1') & ANNOVARtable[samples[1]].str.match('0/1')]
hetList.to_csv("allSharedHetCalls."+inputFile, sep='\t')

#Generic filters for most likely pathogenic
hetList=hetList[hetList['FILTER'].isin(filterTerms)]
hetList=hetList[(hetList[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
hetList=hetList[(hetList[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]
hetList=hetList[~hetList['Func.gene'].isin(notGeneTerms)]
hetList.to_csv("allSharedHetCalls.BestGeneCandidates."+inputFile, sep='\t')

# Compound het calls
mNotfHets=ANNOVARtable[ANNOVARtable[samples[0]].str.match('0/1') & ANNOVARtable[samples[1]].str.contains('|'.join(nullAlelles))]
fNotmHets=ANNOVARtable[ANNOVARtable[samples[0]].str.contains('|'.join(nullAlelles)) & ANNOVARtable[samples[1]].str.match('0/1')]
compHets=pd.concat([mNotfHets, fNotmHets], axis=0, join='outer')
compHets=compHets[compHets['gene.gene'].duplicated(keep=False)] # All possible compHets
compHets=compHets[compHets['FILTER'].isin(filterTerms)]
compHets=compHets[(compHets[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
compHets=compHets[(compHets[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]
compHets=compHets[~compHets['Func.gene'].isin(notGeneTerms)]
compHets=compHets[compHets['gene.gene'].duplicated(keep=False)]  # Re-run the gene filter after the other filters
compHets.to_csv("allcompHetCalls.BestGeneCandidates."+inputFile, sep='\t')

# X-linked
xList=mNotfHets[mNotfHets['chr'].str.match('|'.join(['chrX', 'X']))]
xList=xList[xList['FILTER'].isin(filterTerms)]
xList=xList[(xList[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
xList=xList[(xList[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]
xList=xList[~xList['Func.gene'].isin(notGeneTerms)]
xList.to_csv("allX-linked.BestGeneCandidates."+inputFile, sep='\t')