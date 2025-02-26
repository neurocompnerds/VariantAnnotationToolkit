#!/usr/bin/python3

# Script to split multisample ANNOVAR file
import pandas as pd
import sys, getopt, csv

def usage():
    print(
'''
# splitMultiANNOVAR.py a script to take any ANNOVAR multisample output including our GenomeAnnotationsCombined 
# or BestGeneCandidate files and split them into individual samples.  
# The default behaviour is to remove reference "0/0" calls.
#
# Usage splitMultiANNOVAR.py -i ANNOVAR.table.txt -s sampleList.txt [-k] | [ -h | --help ]
#
# Options:
# -i           /path/to/inputFile    REQUIRED: A multisample ANNOVAR table in tab delimited format
# -s           /path/to/sampleFile   OPTIONAL: A list of specific samples to extract. By default all samples are split into new files
# -k           keep reference calls  OPTIONAL: Default is FALSE.  Add this key if you want to keep 0/0 genotypes
# -h | --help  Displays help                 OPTIONAL: Displays usage information.
#
# Script created by Mark Corbett on 14/03/2019
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
#
'''
         )

# Set initial values
inputFile = ''
sampleFile = ''
keepRefs = False

# Read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],'hi:s:k',['help'])
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
    elif opt in ("-k"):
        keepRefs = True

# Make sure you have what you need
if inputFile == '':
    usage()
    print('Hey, you forgot to tell me which ANNOVAR file to split\n')
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

for s in samples:
    currentSampleList=ANNOVARtable[[s]]
    if keepRefs == False :
        currentSampleList=currentSampleList[~currentSampleList[s].str.match(pat='(0/0) | (0\|0) | (./.)')]
        output = pd.concat([coreTable,currentSampleList], axis=1, join='inner') # Add , sort='False' once Ubuntu is upgraded
    elif keepRefs == True :
        output = pd.concat([coreTable,currentSampleList], axis=1, join='inner') # Add , sort='False' once Ubuntu is upgraded
    output.to_csv(s+".GenomeAnnotationsCombined.txt", sep='\t')

# Output some shell commands that can be cut and pasted into a script or directly to the shell to reorganise the files if needed
print('INFO: You may need to run the following shell commands to sort your tables:\n')
print('for id in ' + ' '.join(samples[0:]) + ' ; do')
print('(')
print('sed \'s,\\([^\\t]*\\) *\\(.*\\),\\2\\t\\1,\' $id.GenomeAnnotationsCombined.txt | sed \'s,^\\t,,g\' | head -n1 > $id.header.txt')
print('sed \'1d\' $id.GenomeAnnotationsCombined.txt | sed \'s,\\([^\\t]*\\) *\\(.*\\),\\2\\t\\1,\' | sed \'s,^\\t,,g\' | sort -k1,1 -k2,2n >> $id.header.txt')
print('mv $id.header.txt $id.GenomeAnnotationsCombined.txt')
print(') &')
print('done')
print('wait\n')

