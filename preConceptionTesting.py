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
# Script created by Mark Corbett on 20/12/2019
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
# Mark Corbett; 06/12/2023; Add in pahsed genotypes and update ANNOVAR field names
#
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
nullAlelles = ['0/0', '0|0'. '\./\.']
# Read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],'hi:m:f:',['help'])
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

# Create the filter function
def bestGeneCandidatesFilter(df):
    df=df[df['FILTER'].isin(filterTerms)]
    df=df[(df[filter005].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.005)).all(axis=1)]
    df=df[(df[filter0001].apply(pd.to_numeric, errors='coerce').fillna(0).lt(0.0001)).all(axis=1)]
    df=df[~df['Func.gene'].isin(notGeneTerms)]
    return df

# Count the number of columns in the ANNOVAR table
with open(inputFile) as f:
    reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
    first_row = next(reader)
    num_cols = len(first_row)-1 # use -1 to set the numbering to 0 based

# Open ANNOVAR table with pandas setting the chr-start-ref-obs column as the index
ANNOVARtable=pd.read_csv(inputFile, sep='\t', index_col = num_cols)
samples = [mumID, dadID]

hetList=ANNOVARtable[ANNOVARtable[samples[0]].str.match(pat = '(0/1)|(0\|1)|(1\|0)') & ANNOVARtable[samples[1]].str.match(pat = '(0/1)|(0\|1)|(1\|0)')]
hetList.to_csv("allSharedHetCalls."+inputFile, sep='\t')

#Generic filters for most likely pathogenic
hetList=bestGeneCandidatesFilter(df=hetList)
hetList.to_csv("allSharedHetCalls.BestGeneCandidates."+inputFile, sep='\t')

# Compound het calls
mNotfHets=ANNOVARtable[ANNOVARtable[samples[0]].str.match(pat = '(0/1)|(0\|1)|(1\|0)') & ANNOVARtable[samples[1]].str.contains('|'.join(nullAlelles))]
mGenes=pd.unique(mNotfHets['Gene.refGene'])
fNotmHets=ANNOVARtable[ANNOVARtable[samples[0]].str.contains('|'.join(nullAlelles)) & ANNOVARtable[samples[1]].str.match(pat = '(0/1)|(0\|1)|(1\|0)')]
fGenes=pd.unique(fNotmHets['Gene.refGene'])
seriesCHgenes=pd.Series(mGenes.tolist() + fGenes.tolist())
chGenes=seriesCHgenes[seriesCHgenes.duplicated()]
compHets=pd.concat([mNotfHets, fNotmHets], axis=0, join='outer')
compHets=compHets[compHets['Gene.refGene'].isin(chGenes)] # All possible compHets
# Independently apply filters to mum and dad lists then filter the CH list
filtmNotfHets=bestGeneCandidatesFilter(df=mNotfHets)
filtfNotmHets=bestGeneCandidatesFilter(df=fNotmHets)
mGenes=pd.unique(filtmNotfHets['Gene.refGene'])
fGenes=pd.unique(filtfNotmHets['Gene.refGene'])
seriesCHgenes=pd.Series(mGenes.tolist() + fGenes.tolist())
chGenes=seriesCHgenes[seriesCHgenes.duplicated()]
compHets=compHets[compHets['Gene.refGene'].isin(chGenes)]
compHets=bestGeneCandidatesFilter(df=compHets)
compHets=compHets[compHets['Gene.refGene'].duplicated(keep=False)]  # Re-run the gene filter after the other filters
compHets.to_csv("allcompHetCalls.BestGeneCandidates."+inputFile, sep='\t')

# X-linked
xList=filtmNotfHets[filtmNotfHets['chr'].str.contains("X", na=False)]
xList.to_csv("allX-linked.BestGeneCandidates."+inputFile, sep='\t')

# ClinVar
cvList=ANNOVARtable[~ANNOVARtable[samples[0,1]].str.contains('|'.join(nullAlelles)) & ANNOVARtable['CLNSIG'].str.contains('|'.join(pathogenicFilter))]
cvList.to_csv("clinVar."+inputFile, sep='\t')
