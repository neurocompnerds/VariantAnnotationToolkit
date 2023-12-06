#!/bin/bash

usage()
{
echo "## AnnovarGenomeSummaryCombo.sh -c ComboFile.csv -g GenomeSummaryFile.csv -o OUTPREFIX [-h | --help]
# Script is specific to the annovar pipeline in neurogenetics this will eventually be merged with that script.
# For combining two output files:
# OUTPREFIX.annotated.genome_summmary.fixed.csv
# OUTPREFIX.snps.filtered.avinput.combo.txt
#
# Options:
# -c Combo file 
# -g ANNOVAR genome_summary or multianno file file
# -o OUTPREFIX		Usually a DNA number that forms the prefix of the file (required)
# -h | --help	Prints this message
#
# Example:
# AnnovarGenomeSummaryCombo.sh -c ComboFile.csv -g multianno.csv
#
# History:
# Created by: Mark Corbett, 16/07/2013, mark.corbett@adelaide.edu.au
# Modified: (Date; Name; Modification)
# 02/12/2013; Mark Corbett; Added generic file specification.  Now calculates number of cols in genome summary file
# 31/12/2013; Mark Corbett; Changed evaluation of total line numbers was working on previous version of Ubuntu
# 23/04/2014; Mark Corbett; Minor shuffle of column numbers to account for addition of new databases in ANNOVAR
# 14/05/2014; Mark Corbett; Added generation of a variant key
# 28/08/2014; Mark Corbett; Update to suit July 2014 ANNOVAR. Move script to version 2
# 08/09/2014; Thuong Ha; Minor shuffle of column numbers to account for addition of new databases in ANNOVAR
# 31/10/2014; Mark Corbett; Minor shuffle of column numbers to account for addition of new databases in ANNOVAR
# 31/10/2014; Mark Corbett; Fix missing VCF fields
# 27/05/2015; Mark Corbett; Update to suit March 2015 ANNOVAR
# 04/01/2016; Mark Corbett; Check if input VCF was bgzipped before extracting missing VCF fields
# 11/04/2016; Ali Gardner; Update for use with new databases in Annovar
# 11/07/2018; Mark Corbett; Change ordering and headers of new ANNOVAR output
# 15/03/2019; Mark corbett; Updated script to handle extra Otherinfo columns that get added using table_annovar.pl -allsamples flag
# 11/06/2019; Ali Gardner; Updated some databases in Annovar
# 27/03/2020; Mark Corbett; Fix the annoying otherInfo labels
"
}

while [ "$1" != "" ]; do
	case $1 in
		-c )	shift
								ComboFile=$1
								;;
		-g )	shift
								GenomeFile=$1
								;;
		-o )	shift						
								OUTPREFIX=$1
								;;
		-h | --help )			usage
								exit 1
								;;
		* )						usage
								exit 1
	esac
	shift
done

if [ -z "$ComboFile" ]; then #If combo file is empty then do not proceed
	echo "#ERROR: You need to specify the name of the ANNOVAR combo file"
	usage
	exit 1
fi

if [ -z "$GenomeFile" ]; then #If genome summary file is empty then do not proceed
	echo "#ERROR: You need to specify the name of the ANNOVAR genome summary file"
	usage
	exit 1
fi

if [ -z "$OUTPREFIX" ]; then #If no outprefix specified then do not proceed
	echo "#ERROR: You need to specify a prefix for the file output"
	usage
	exit 1
fi

# Sanity Check
LinesA=$(wc -l < $ComboFile | sed 's/[^0-9]*//g')
LinesB=$(wc -l < $GenomeFile | sed 's/[^0-9]*//g')

if [ "$LinesA" -eq "$LinesB" ]; then
	echo "# OK We're good to go. Now combining the files"
	else
	echo "#ERROR: The number of lines in both files is not the same.
If the VCF was from complete genomics try running the output of awk -F "," '\$5 !~/[\<\[\]]/' $ComboFile to a new combo file and use that instead"
	wc -l $ComboFile
	wc -l $GenomeFile
	exit 1
fi

# Fix up field delimiters 
awk -F\" '{for(i=1;i<=NF;i+=2) {gsub(",", "\t", $i)}}1' OFS= $ComboFile > $ComboFile.txt
awk -F\" '{for(i=1;i<=NF;i+=2) {gsub(",", "\t", $i)}}1' OFS= $GenomeFile > $GenomeFile.txt
NColsGenomeFile=$(awk -F "\t" 'FNR == 2 {print NF}' $GenomeFile.txt)

# Combine columns (specific to how ANNOVAR is run so may or may not work perfectly)
cut -f 1-5 $ComboFile.txt > $OUTPREFIX.tmp.1.txt # chr,start,end,ref,obs
cut -f 6-19 $GenomeFile.txt > $OUTPREFIX.tmp.2a.txt # Func.gene,Gene,GeneDetail,ExonicFunc,AAChange,Conserved,SegDup,ESP6500siv2_ALL,ExAC.r0.1.filtered,1000g2014oct_ALL,UK10K,cg69,Wellderly,PopFreqMax
cut -f 14-16 $ComboFile.txt > $OUTPREFIX.tmp.2b.txt # exac03,gnomad_exome,gnomad_genome
cut -f 20-26 $GenomeFile.txt > $OUTPREFIX.tmp.3a.txt # snp150,snp138NonFlagged, clinvar (5cols)
cut -f 6,7,8,10,12 $ComboFile.txt > $OUTPREFIX.tmp.3b.txt # DDG2P,EpilepsyGene,CPGene,IDGene,MCDGene
cut -f 9,11,13 $ComboFile.txt > $OUTPREFIX.tmp.3c.txt # GDIScores,LoFToolScores,RVISExACscores
cut -f 27-119 $GenomeFile.txt > $OUTPREFIX.tmp.4.txt # All the other annotations
cut -f 120-$NColsGenomeFile $GenomeFile.txt > $OUTPREFIX.tmp.5.txt # VCF info

# Replace Otherinfo labels with VCF fields"
if [ "${OUTPREFIX:(-2)}" = "gz" ]; then
	VCFFIELDS=$(zcat $OUTPREFIX | grep -m1 \#CHROM) #Assumes $OUTPREFIX = a bgzipped vcf file which might be the case if fed from ANNOVARv3.sh
        NcolsVCFfields=$(zcat $OUTPREFIX | grep -m1 \#CHROM | awk -F "\t" 'FNR == 1 {print NF}')
else
	VCFFIELDS=$(grep -m1 \#CHROM $OUTPREFIX) #Assumes $OUTPREFIX = a vcf file which might be the case if fed from ANNOVARv3.sh
        NcolsVCFfields=$(grep -m1 \#CHROM $OUTPREFIX | awk -F "\t" 'FNR == 1 {print NF}')
fi

NColsTmp5=$(awk -F "\t" 'FNR == 2 {print NF}' $OUTPREFIX.tmp.5.txt)

if [ "$NcolsVCFfields" -eq "$NColsTmp5" ]; then
    printf "$VCFFIELDS\n" > $OUTPREFIX.tmp.6.txt
    sed '1d' $OUTPREFIX.tmp.5.txt >> $OUTPREFIX.tmp.6.txt
    mv $OUTPREFIX.tmp.6.txt $OUTPREFIX.tmp.5.txt
else # A bit risky but should work unless something weird happened
    printf "Allele_Fraction\t$VCFFIELDS\n" > $OUTPREFIX.tmp.6.txt
    cut -f 1,4-$NColsTmp5 $OUTPREFIX.tmp.5.txt | sed '1d'>> $OUTPREFIX.tmp.6.txt # Assume -withzyg info has been added. Keep the Allele fraction but chuck out the Depth and QUAL data
    mv $OUTPREFIX.tmp.6.txt $OUTPREFIX.tmp.5.txt
fi

# Create a key for variant lookup chr-start-ref-obs
cut -f 1,2,4,5 $OUTPREFIX.tmp.1.txt | sed 's,\t,-,g' > $OUTPREFIX.tmp.key.txt  
paste $OUTPREFIX.tmp.1.txt $OUTPREFIX.tmp.2a.txt $OUTPREFIX.tmp.2b.txt $OUTPREFIX.tmp.3a.txt $OUTPREFIX.tmp.3b.txt $OUTPREFIX.tmp.3c.txt $OUTPREFIX.tmp.4.txt $OUTPREFIX.tmp.5.txt $OUTPREFIX.tmp.key.txt > $OUTPREFIX.GenomeAnnotationsCombined.txt

# Clean up temp files
rm $OUTPREFIX.tmp.1.txt $OUTPREFIX.tmp.2a.txt $OUTPREFIX.tmp.2b.txt $OUTPREFIX.tmp.3a.txt $OUTPREFIX.tmp.3b.txt $OUTPREFIX.tmp.3c.txt $OUTPREFIX.tmp.4.txt $OUTPREFIX.tmp.5.txt $ComboFile.txt $GenomeFile.txt
