#!/bin/bash
# ANNOVARv3.sh my.vcf

# Variables that usually don't need changing once set for your system
AnnovarPATH=/opt/annovar # Where the Annovar program is
SCRIPTPATH=/home/neuro/scripts/VariantAnnotationToolkit # Where the python & perl scripts for Annovar or other general scripts are
BUILD=Hs38DH # Genome build used by ANNOVAR either hg18 or $BUILD. This will also be incorporated into file names
AV_INPUT=$1.avinput
AV_DB=/opt/annovar/humandb/hg38/

usage()
{
echo "# $0 A script to annotate variants using ANNOVAR
# Usage $0 VCF.filename.vcf | -h
#
#	-h --help	Prints this message
#
# Example:
# $0 GATK.filtered.snps.vcf
#
# Paths currently set
# AnnovarPATH=$AnnovarPATH
# SCRIPTPATH=$SCRIPTPATH
# BUILD=$BUILD
# AV_DB=$AV_DB
#
# History:
# Original: Mark Corbett
# mark.corbett@adelaide.edu.au
# Modified (date; name; description)
# 25/11/2013; Mark Corbett; Extract Annovar script from BWA-GATK-Annovar pipeline and update to latest version
# 29/11/2013; Mark Corbett; Standardise variable assignments. Fix clean up script creation.  Add usage.
# 17/03/2014; Mark Corbett; Add call to AnnovarGenomeSummaryCombo.sh
# 22/04/2014; Mark Corbett; Update to include new databases
# 04/08/2014; Mark Corbett; Change splicing threshold from default 2 to 12
# 04/08/2014; Mark Corbett; Change to using table_annovar.pl as summarise_annovar.pl is deprecated
# 08/09/2014; Thuong Ha; Update to include new databases
# 31/10/2014; Mark Corbett; Add scary number of variants in ExAC database v0.1
# 27/05/2015; Mark Corbett; Bring up to date with the March 2015 version of ANNOVAR
# 04/01/2016; Mark corbett; Replace vcf4old with vcf4 as new GATK outputs this format
# 11/04/2016; Ali Gardner; Update new databases used with Annovar
# 11/07/2018; Mark Corbett; Organise table order see also AnnovarGenomeSummaryCombo.v2.sh
# 03/09/2018; Ali Gardner; Update most recent ID gene list
# 01/03/2019; Mark Corbett; Update to match Ali's December 2018 version.  Now annotates multisample VCF.
# 15/03/2019; Mark Corbett; Remove Eye gene list and replace with CP gene bed file.
# 11/06/2019; Ali Gardner; Updated some databases in Annovar
# 02/06/2020; Ali Gardner; Updated for use with hg38 vcfs
" 
}

if [ -z "$1" ] ; then
	usage
	exit
fi

case $1 in
	-h | --help )	usage
			exit
			;;
esac

# Create clean up script
cp $SCRIPTPATH/BWA-Picard-GATK-CleanUp.sh $1.CleanUp.sh # Items are added to this script that can then be run later to delete all redundant files

#Annovar - Conversion to Annovar file format. Will skip if file exists so if you think the file is dodgy need to delete it
if [ ! -f $AV_INPUT ]; then
	perl $AnnovarPATH/convert2annovar.pl --format vcf4 --allsample --withfreq --includeinfo $1 > $AV_INPUT
fi

# Run Annovar
# This section is starting to slow down might need to consider running some filters seperately and in parallel in future
perl $AnnovarPATH/table_annovar.pl -thread 12 $AV_INPUT $AV_DB/ \
--buildver $BUILD \
--remove \
--protocol gene,phastConsElements100way,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,gnomad30_genome,avsnp151,clinvar_20240611,dbnsfp41a,\
dbnsfp47a_interpro,dbscsnv11,regsnpintron,spliceai_filtered,gene4denovo201907,\
dgvMerged,evoCpg,cpgIslandExt,evofold,gwasCatalog,switchDbTss,targetScanS,vistaEnhancers,wgRna \
--operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,r \
--otherinfo \
--nastring . \
--csvout \
-out $1.snps_annotated > $1.pipeline.log 2>&1

# ExAC and gnomAD
perl $AnnovarPATH/annotate_variation.pl --filter --buildver $BUILD --thread 8 --dbtype exac03 $AV_INPUT $AV_DB/ >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl --filter --buildver $BUILD --thread 8 --dbtype gnomad211_exome $AV_INPUT $AV_DB/ >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl --filter --buildver $BUILD --thread 8 --dbtype gnomad312_genome $AV_INPUT $AV_DB/ >> $1.pipeline.log 2>&1

# Generate special regionanno files
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile Epilepsy_hg38_genes_Mar22.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","EpilepsyGene",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with Epilepsy_gene
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_EpilepsyGene  #Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile ID_hg38_genes_sept24.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","IDGene",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with ID_gene
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_IDGene  # Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile hg38_CPgenesMay2020.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","CPgene",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with CPgene
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_CPgene  #Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile CMHg38_04062020AG_50bp.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","MCDGene",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with MCDgene
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_MCDGene  # Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile Kat_hg38_steroids_bed.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","steroidGenes",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with steroidGenes
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_steroidGenes  # Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile burden_panel_strict_kj.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","burdenAsdpanel",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with burdenAsdpanel
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_burdenAsdpanel  # Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed -bedfile mega_panel_asd_id_panelapp_kj.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
awk '{gsub("bed","megaAsdpanel",$1)}1' $AV_INPUT.$BUILD\_bed > $AV_INPUT.$BUILD.bed_new   # replace bed with megaAsdpanel
awk -v OFS="\t" '$1=$1' $AV_INPUT.$BUILD.bed_new > $AV_INPUT.$BUILD\_megaAsdpanel  # Replace whitespaces with tabs in linux
rm $AV_INPUT.$BUILD.bed_new $AV_INPUT.$BUILD\_bed

perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype DDG2P --colsWanted 5,6,7,8,9 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype LoFToolScores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype RVISExACscores --colsWanted 4,5 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype GDIScores --colsWanted 4,5,6,7 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
#perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype oe_upper_mis_scores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
#perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype oe_upper_LoF_scores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype pLIscores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype Zscore --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1

# Combine all variant files
python $SCRIPTPATH/annovar_combine_csv.py $AV_INPUT $AV_INPUT.$BUILD* > $AV_INPUT.combo.csv
#perl $SCRIPTPATH/vcfCSVFix.pl $AV_INPUT.combo.csv > $AV_INPUT.combo.txt

# Create genome summary combo file
$SCRIPTPATH/AnnovarGenomeSummaryCombo.v2_for_hg38.sh -c $AV_INPUT.combo.csv -g $1.snps_annotated.$BUILD\_multianno.csv -o $1

# Add files to clean up script
echo "rm $AV_INPUT" >> $1.CleanUp.sh
echo "rm $AV_INPUT.$BUILD\_*" >> $1.CleanUp.sh
	
