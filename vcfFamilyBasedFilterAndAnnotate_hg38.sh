#!/bin/bash
usage ()
{
echo "# For processing VCFs and doing some intial filtering of variants
# Requires:
# gnu parallel, vcftools
# Example:
# $0 -f SMITH -i ~/MyVCF/Folder -v MyVcf -a DNA1,DNA2 -c ~/MyFolder/ControlList.csv [-e DNA3] | [-h | --help]
# 
# Options:
# -f FAMILY		Text of name of family to sort (a subdirectory of this name will be created in your input directory)
# -i /path/to/input_dir		Where your main VCF is stored
# -v VCF		The name of the VCF file to query
# -a list,of,affected		Samples that you would list under + in vcf contrast
# -c /path/to/list/of/control.txt		Samples you would list under - in vcf contrast listed in a file separated by commas
# -e extra control samples to annotate	List of samples from the controls to include annotations for
# -h | --help		Displays this message
#
# Mark Corbett; 16/05/2014; mark.corbett at adelaide.edu.au
# Modified (Date; Name; Description)
# 26/05/2014; Mark Corbett; Fixed up to workable reasonably generic script
# 15/09/2014; Mark Corbett; Update to suit ANNOVARv3.sh script
# 03/02/2015; Mark Corbett; Update to put in strict variant filter
# 15/09/2015; Mark Corbett; Parallelize variant annotation for loop
# 21/01/2021; Ali Gardner; change columns & values for filtering line 126 awk as no UK10K or Wellderly, now gnomad_exome & genome
#
"
}

# Set variables
while [ "$1" != "" ]; do
	case $1 in
		-f )			shift
						Family=$1
						;;
		-i )			shift
						inputDir=$1
						;;
		-v )			shift
						VCF=$1
						;;
		-a )			shift
						Affected=$1
						arrAffected=$(echo $Affected | tr "," "\n") # Translate comma to newlines to create array
						;;
		-c )			shift
						ControlFile=$1
						Controls=$(cat $ControlFile)
						;;
		-e )			shift
						Extras=$1
						arrExtras=$(echo $Extras | tr "," "\n")
						arrAnnotate=$(echo $arrAffected $arrExtras)
						;;

		-h | --help )	usage
						exit 1
						;;
		* )				usage
						exit 1
	esac
	shift
done
if [ -z "$Family" ]; then #If no family specified then do not proceed
	usage
	echo "#ERROR: You need to specify a family name to test."
	exit 1
fi
if [ -z "$inputDir" ]; then #If no input directory specified then do not proceed
	usage
	echo "#ERROR: You need to tell me where to find your VCF file."
	exit 1
fi
if [ -z "$VCF" ]; then #If VCF file not specified then do not proceed
	usage
	echo "#ERROR: You need to tell me the name of your VCF file."
	exit 1
fi
if [ -z "$Affected" ]; then #If list of affected samples is not specified then do not proceed
	usage
	echo "#ERROR: You need to tell me which DNA id to test as affected in a comma separated list."
	exit 1
fi
if [ -z "$ControlFile" ]; then #If no list of control file then do not proceed
	usage
	echo "#ERROR: You need to tell me which DNA id to test as controls in a comma separated file."
	exit 1
fi
if [ -z "$arrAnnotate" ]; then #If extras not specified then only annotate Affecteds
	arrAnnotate=$arrAffected
fi

famDir=$inputDir/$Family
if [ ! -d $famDir ]; then # If directory doesn't exist then create it.
	mkdir -p $famDir # -p forces the entire path to be created otherwise epic failure ensues
fi

## Start script ##

cd $famDir
vcf-contrast -n +$Affected -$Controls $inputDir/$VCF > $Family.common.vcf
for Sample in $arrAnnotate; do
	( 
	mkdir $Sample
	~/Documents/Scripts/gitHub/VariantAnnotationToolkit/VCFSplitWithGATKv4_hg38.sh $Family.common.vcf $Sample
	mv $Sample.$Family.common.* $famDir/$Sample/
	cd $Sample
	~/Documents/Scripts/gitHub/VariantAnnotationToolkit/ANNOVARv3_for_hg38.sh $Sample.$Family.common.vcf
	./$Sample.$Family.common.vcf.CleanUp.sh
	# Need this next line as normally key is removed in the ANNOVAR script
	cut -f 1,2,4,5  $Sample.$Family.common.vcf.GenomeAnnotationsCombined.txt | sed 's,\t,-,g' > $Sample.key.txt
	cd $famDir
	) &
done
wait

# Implementing nested for loop to reduce samples to common lists in all affecteds against all other affecteds
# This next bit rocks for creating the intersection of both samples - way faster than doing it with excel / libreoffice
cd $famDir
for Sample in $arrAffected; do # Test each sample
cp $famDir/$Sample/$Sample.$Family.common.vcf.GenomeAnnotationsCombined.txt $Family.GenomeAnnotationsCombined.$Sample.txt
	for Aff in $arrAffected; do # Against all other samples (and self but this doesn't matter)
		parallel --pipe --block 10M grep -F -w -f $famDir/$Aff/$Aff.key.txt < $Family.GenomeAnnotationsCombined.$Sample.txt > tmp.$Sample.txt
		mv tmp.$Sample.txt $Family.GenomeAnnotationsCombined.$Sample.txt
		head -n1 $Family.GenomeAnnotationsCombined.$Sample.txt > $Family.BestGeneCandidates.$Sample.txt # Save the header row
		awk -F"\t" '$13 < 0.005 && $14 < 0.005 && $15 < 0.001 && $16 < 0.001 && $17 < 0.001 {print}' $Family.GenomeAnnotationsCombined.$Sample.txt |\
		grep -Fwvf ~/Documents/Scripts/gitHub/VariantAnnotationToolkit/Func.Gene.Remove.txt - >> $Family.BestGeneCandidates.$Sample.txt # Variant frequency filters, exonic, splicing and PASS
	done
done

