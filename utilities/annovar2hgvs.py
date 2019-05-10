#!/usr/bin/python3

# Script to convert ANNOVAR cDNA c.N##N to HGVS c.##N>N format

import sys, getopt

def usage():
    print(
'''
# annovar2hgvs.py
# Converts ANNOVAR cDNA c.N##N to HGVS c.##N>N format.
# Needs a file with the cDNA coordinates from the AAchange column of ANNOVAR
#
# Usage annovar2hgvs.py -i variants.txt | [ -h | --help ]
#
# Options:
# -i           /path/to/inputFile    REQUIRED: A multisample ANNOVAR table in tab delimited format
# -h | --help  Displays help                 OPTIONAL: Displays usage information.
#
# Script created by Mark Corbett on 10/05/2019
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
#
'''
         )

# Set initial values
inputFile = ''

# Read command line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],'hi:',['help'])
except getopt.GetoptError:
    usage
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-i"):
        inputFile = arg

# Make sure you have what you need
if inputFile == '':
    usage()
    print('Hey, you forgot to tell me which variants to convert\n')
    sys.exit(2)

# Count the number of columns in the ANNOVAR table
with open(inputFile) as f:
    annovarList = f.readlines()

for annovar in annovarList:
    outFile = open("hgvs."+inputFile, 'a')
    hgvs = annovar[:2]+annovar[3:-2]+annovar[2]+">"+annovar[-2:]  # I think should include the \n character
    outFile.write(hgvs)
    
outFile.close()
