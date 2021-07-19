#!/bin/bash
usage()
{
echo "
# VCFSplitWithGATKv3.sh
# Script for extracting a single sample from a multi-sample VCF
# Requires: GATKv3.x
# Usage: $0 /path/to/VCF.file.vcf SAMPLE
#
# Variables:
# \$1	/path/to/VCF.file.vcf	REQUIRED. Your VCF file, it can be bgzipped with .gz extension
# \$2	SAMPLE	REQUIRED.Name of the sample you want to extract as it is in the VCF.  (See note below).
# \$3		Optional. change reference file from default hg19_1stM_unmask_ran_all.fa
#
# How can I know which samples are in my VCF?
# bcftools query -l /path/to/VCF.file.vcf (It can be bgzipped with .gz extension)
#
# Original: Mark Corbett, 12/02/2014
# Modified: (Date; Name; Description)
# 01/09/2014: Ali: So points to GATK ver 3
# 15/09/2015: Mark: Increase threads and memory
# 04/01/2016: Mark: Fix file naming so that \$VCF can contain the path. Add usage.
# 07/03/2016: Mark: Add provision for using a different GATK index.
# 16/06/2020: Ali: change for GATK4 & appropriate commands & ucsc ref."
}

if [ -z "$1" ]; then
	usage
	exit 1
fi
case $1 in
	-h | --help )	usage
			exit
			;;
esac

## Set Variables ##
VCF=$1 # list VCF file as first agument
OUTPREFIX=$2 #List sample name from VCF file to extract as the second argument
GATKINDEX=$3

# Variables that usually don't need changing once set for your system
GATKPATH=/home/neuro/Documents/Ali/gatk/build/libs # Where the GATK program.  Be mindful that GATK is under rapid development so things may change over time!
GATKREFPATH=~/Public/RefSeqIndexAllPrograms #Refseq index library locations
if [ -z "$GATKINDEX" ]; then
	GATKINDEX=ucsc.hg19.fasta # name of the genome reference
fi

## Start of the script ##
echo "GATK index is $GATKINDEX"
baseFile=$(basename $VCF)

# SelectVariants by sample name
java -Xmx16g -jar $GATKPATH/gatk.jar SelectVariants \
--variant $VCF \
-R $GATKREFPATH/$GATKINDEX \
-O $OUTPREFIX.$baseFile \
-sn $OUTPREFIX \
--exclude-non-variants > $OUTPREFIX.pipeline.log 2>&1

