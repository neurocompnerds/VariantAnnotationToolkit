#!/bin/bash
# ANNOVARv3_for_hg38.sh my.vcf(.gz)

# Variables that usually don't need changing once set for your system
SCRIPTPATH="$(dirname "$(readlink -f "$0")")" # Where the python & perl scripts for Annovar or other general scripts are
BUILD=Hs38DH # Genome build used by ANNOVAR either hg18 or $BUILD. This will also be incorporated into file names
source $SCRIPTPATH/config/environment.cfg
source $SCRIPTPATH/config/environment.local.cfg 2>/dev/null # Load local settings if they exist
AV_INPUT=$1.avinput
AV_DB=/opt/annovar/humandb/hg38/

usage()
{
echo "# $0 A script to annotate variants using ANNOVAR
# Usage $0 VCF.filename.vcf.gz [-b build -g| -h
#
#	-h --help	Prints this message
#
# Example:
# $0 GATK.filtered.snps.vcf.gz
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
# 26/02/2025; Mark Corbett; Bring changes made to other test / draft scripts to this script
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

#Annovar - Conversion to Annovar file format. Will skip if file exists so if you think the file is dodgy you need to delete it
if [ ! -f $AV_INPUT ]; then
	perl $AnnovarPATH/convert2annovar.pl --format vcf4 --allsample --withfreq --includeinfo $1 > $AV_INPUT
fi

# Run Annovar
# This section is starting to slow down might need to consider running some filters seperately and in parallel in future
perl $AnnovarPATH/table_annovar.pl -thread 8 $AV_INPUT $AV_DB/ \
--buildver $BUILD \
--remove \
--protocol ensGene,\
phastConsElements100way,genomicSuperDups,\
esp6500siv2_all,1000g2015aug_all,gnomad30_genome,avsnp151,clinvar_latest,dbnsfp54a,dbscsnv11,regsnpintron,spliceai_filtered,gene4denovo201907,\
dgvMerged,evoCpg,cpgIslandExt,evofold,gwasCatalog,switchDbTss,targetScanS,vistaEnhancers,wgRna \
--operation g,r,r,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,r \
--otherinfo \
--nastring . \
--csvout \
-out $1.snps_annotated > $1.pipeline.log 2>&1

# ExAC and gnomAD
perl $AnnovarPATH/annotate_variation.pl --filter --buildver $BUILD --thread 8 --dbtype exac03 $AV_INPUT $AV_DB/ >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl --filter --buildver $BUILD --thread 8 --dbtype gnomad211_exome $AV_INPUT $AV_DB/ >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl --filter --buildver $BUILD --thread 8 --dbtype gnomad312_genome $AV_INPUT $AV_DB/ >> $1.pipeline.log 2>&1

# Generate special regionanno files
for Panel in Epilepsy ID CP ASD; do
  perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype bed  --colsWanted 4 -bedfile ${BUILD}_${Panel}_latest.bed $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
  sed 's,^bed,'${Panel}',g' "${AV_INPUT}.${BUILD}_bed" > "${AV_INPUT}.${BUILD}_${Panel}"
  rm "${AV_INPUT}.${BUILD}_bed" 
done

perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype DDG2P --colsWanted 5,6,7,8,9 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype LoFToolScores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype RVISExACscores --colsWanted 4,5 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype GDIScores --colsWanted 4,5,6,7 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype oeUpperMisScores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
perl $AnnovarPATH/annotate_variation.pl -regionanno -buildver $BUILD -dbtype oeUpperLoFScores --colsWanted 4 $AV_INPUT $AV_DB >> $1.pipeline.log 2>&1
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
