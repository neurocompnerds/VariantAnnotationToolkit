#!/bin/bash

usage()
{
echo "# $0 A script to append keys to best candidate files for manual filtering
# This script is specifically for CP twin analysis
# Requirements: Run ANNOVARv3.sh on a multisample r single sample vcf for each family member.
# use the splitMultiANNOVAR.py script to break GenomeAnnotationsCombined.txt files out to individuals
# You may need to generate the key files.
#
# Usage $0 -f familyID | [-h | --help]
#
#	-f		    REQUIRED: ID of the family (included in the file name eg. V2038
#	-h | --help	Prints this message
#
# Example:
# $0 -f V2038
#
# History:
# Original: forked from trioKeyMatchingAfterANNOVAR.sh by Mark Corbett and Nandini 05/03/2019
# mark.corbett at adelaide.edu.au
# Modified (date; name; description)
# 16/03/2019; Mark Corbett; Minor bug fixes to ake teh script more tolerant of file names and missing key files
#
"
}

# Read flags from the command line
while [ "$1" != "" ]; do
    case $1 in
        -f )             shift
                         familyID=$1
                         ;;
        -h | --help )    usage
                         exit 0
                         ;;
        * )              usage
                         exit 1
        esac
        shift
done

# Check the script has everything it needs
if [ -z "$familyID" ] ; then
    usage
    echo "You need to specify the ID of the family you want to filter."
    exit 1
fi

workDir=$(pwd)

# Assign family ID
probandID=$familyID\-1
dadID=$familyID\-2
mumID=$familyID\-3
sibID=$familyID\-4

## Start the script
if [ -f $probandID.*GenomeAnnotationsCombined.txt ] ; then
    awk '{OFS="\t"; print $NF}' $probandID.*GenomeAnnotationsCombined.txt > $probandID.keys.txt &
    awk '{OFS="\t"; print $NF}' $sibID.*GenomeAnnotationsCombined.txt > $sibID.keys.txt &
    awk '{OFS="\t"; print $(NF-1), $NF}' $mumID.*GenomeAnnotationsCombined.txt > $mumID.keys.txt &
    awk '{OFS="\t"; print $(NF-1), $NF}' $dadID.*GenomeAnnotationsCombined.txt > $dadID.keys.txt
    wait
else
    usage
    echo "I couldn't find the file $probandID.*GenomeAnnotationsCombined.txt.
    Sorry I'm a bit particular about file names and I need to be executed from the folder that contains the proband's GenomeAnnotationsCombined file.
    Your current directory has these files:"
    ls
    exit 1
fi

## Assign proband and sib to names - to run them one after another

for name in $probandID $sibID ; do
    (
    echo "Finding matching keys for $name"
    grep -Fwf $name.keys.txt $mumID.keys.txt | sort -k2 > $mumID.$name.keys.txt &
    grep -Fwf $name.keys.txt $dadID.keys.txt | sort -k2 > $dadID.$name.keys.txt

    wait

    awk '{OFS="\t"; print $0, NR}' $name.keys.txt | sort -k1,1 > $name.sort.keys.txt
    join -a 1 -e "." -o '0,2.1,1.2' -1 1 -2 2 $name.sort.keys.txt $mumID.$name.keys.txt | sort -k3 -n | cut -f2 -d " " > $mumID.$name.column.txt &
    join -a 1 -e "." -o '0,2.1,1.2' -1 1 -2 2 $name.sort.keys.txt $dadID.$name.keys.txt | sort -k3 -n | cut -f2 -d " " > $dadID.$name.column.txt
    wait

    paste $name.*GenomeAnnotationsCombined.txt $mumID.$name.column.txt $dadID.$name.column.txt > $name.parentKeys.GenomeAnnotationsCombined.txt
    head -n1 $name.parentKeys.GenomeAnnotationsCombined.txt > $name.parentKeys.GenomeAnnotationsCombined.header.txt

    # Filter based on inheritance models
    cp $name.parentKeys.GenomeAnnotationsCombined.header.txt $name.parentKeys.GenomeAnnotationsCombined.dn.txt
    cp $name.parentKeys.GenomeAnnotationsCombined.header.txt $name.parentKeys.GenomeAnnotationsCombined.ch.txt
    cp $name.parentKeys.GenomeAnnotationsCombined.header.txt $name.parentKeys.GenomeAnnotationsCombined.ibd.txt
    cp $name.parentKeys.GenomeAnnotationsCombined.header.txt $name.parentKeys.GenomeAnnotationsCombined.Xl.txt

    # de novo variants
    awk '$(NF-3) ~ /^0\/1/ && $(NF-1) ~ /^\./ && $NF ~ /^\./ {OFS="\t"; print}' $name.parentKeys.GenomeAnnotationsCombined.txt >> $name.parentKeys.GenomeAnnotationsCombined.dn.txt
    cut -f7 $name.parentKeys.GenomeAnnotationsCombined.dn.txt | sort | uniq > $name.CH.gene.list.1.txt  # These variants also need to be saved for possible CH model

    #X-linked
    awk '$1 == "X" && $(NF-3) ~ /^1\/1/ && $(NF-1) !~ /^1/ && $NF ~ /^\./ {OFS="\t"; print}' $name.parentKeys.GenomeAnnotationsCombined.txt >> $name.parentKeys.GenomeAnnotationsCombined.Xl.txt

    # IBD
    awk '$(NF-3) ~ /^1\/1/ && $(NF-1) !~ /^1/ && $NF !~ /^1/ {OFS="\t"; print}' $name.parentKeys.GenomeAnnotationsCombined.txt >> $name.parentKeys.GenomeAnnotationsCombined.ibd.txt

    # Potential compound heterozygous
    # Het variants from mum
    awk '$(NF-3) ~ /^0\/1/ && $(NF-1) ~ /^0\/1/ && $NF ~ /^\./ {OFS="\t"; print}' $name.parentKeys.GenomeAnnotationsCombined.txt > $name.parentKeys.GenomeAnnotationsCombined.chM.txt
    cut -f7 $name.parentKeys.GenomeAnnotationsCombined.chM.txt | sort | uniq >> $name.CH.gene.list.1.txt
    # Het variants from dad
    awk '$(NF-3) ~ /^0\/1/ && $(NF-1) ~ /^\./ && $NF ~ /^0\/1/ {OFS="\t"; print}' $name.parentKeys.GenomeAnnotationsCombined.txt > $name.parentKeys.GenomeAnnotationsCombined.chF.txt
    cut -f7 $name.parentKeys.GenomeAnnotationsCombined.chF.txt | sort | uniq >> $name.CH.gene.list.1.txt
    grep -v Gene.gene $name.CH.gene.list.1.txt | sort | uniq -d > $name.CH.gene.list.2.txt
    cat $name.parentKeys.GenomeAnnotationsCombined.dn.txt $name.parentKeys.GenomeAnnotationsCombined.chM.txt $name.parentKeys.GenomeAnnotationsCombined.chF.txt | sed '1d' | parallel --will-cite --pipe --block 10M grep -wf $name.CH.gene.list.2.txt >> $name.parentKeys.GenomeAnnotationsCombined.ch.txt

    # Clean up
    rm $mumID.$name.column.txt $dadID.$name.column.txt $mumID.$name.keys.txt $dadID.$name.keys.txt  $name.parentKeys.GenomeAnnotationsCombined.header.txt $name.parentKeys.GenomeAnnotationsCombined.chM.txt $name.parentKeys.GenomeAnnotationsCombined.chF.txt $name.CH.gene.list.1.txt $name.CH.gene.list.2.txt

    # Pull-out column 3 (from right) from each file created based on the inheritance models
    awk '{print $(NF-2)}' $name.parentKeys.GenomeAnnotationsCombined.dn.txt > $name.AllKeys.dn.txt
    awk '{print $(NF-2)}' $name.parentKeys.GenomeAnnotationsCombined.ch.txt > $name.AllKeys.ch.txt
    awk '{print $(NF-2)}' $name.parentKeys.GenomeAnnotationsCombined.ibd.txt > $name.AllKeys.ibd.txt
    awk '{print $(NF-2)}' $name.parentKeys.GenomeAnnotationsCombined.Xl.txt > $name.AllKeys.Xl.txt

    echo $name Key Matching done
    ) &
done
wait

rm $mumID.keys.txt $dadID.keys.txt

## Sort the keys into unique and shared

for model in dn ch ibd Xl ; do
    echo sorting $familyID of $model to unique and shared keys
    cat $probandID.AllKeys.$model.txt $sibID.AllKeys.$model.txt | sort | uniq -u > $familyID.uniqueKeys.$model.txt
    cat $probandID.AllKeys.$model.txt $sibID.AllKeys.$model.txt | sort | uniq -d > $familyID.sharedKeys.$model.txt
    echo sorting $familyID of $model
done

## Uniquekeys vs Annotated Genome followed by variant filtering

for name in $probandID $sibID ; do
    for model in dn ch ibd Xl  ; do
        (
        head -n1 $name.parentKeys.GenomeAnnotationsCombined.$model.txt > $name.uniqueKeys.GenomeAnnotationsCombined.$model.txt #header
        grep -Fwf $familyID.uniqueKeys.$model.txt $name.parentKeys.GenomeAnnotationsCombined.$model.txt >> $name.uniqueKeys.GenomeAnnotationsCombined.$model.txt
        ##variant filtering
        head -n1 $name.uniqueKeys.GenomeAnnotationsCombined.$model.txt > $name.uniqueKeys.BestGeneCandidates.$model.txt #header
        awk -F"\t" '$13 < 0.005 && $14 < 0.001 && $15 < 0.005 && $16 < 0.005 && $18 < 0.01 && $20 < 0.001 && $21 < 0.001 && $22 < 0.001 && $(NF-6)=="PASS" {print}' $name.uniqueKeys.GenomeAnnotationsCombined.$model.txt |\
        grep -Fwvf ~/Documents/Scripts/local/Func.Gene.Remove.txt - >> $name.uniqueKeys.BestGeneCandidates.$model.txt # Variant frequency filters, exonic, splicing and PASS
		) &
	done
done
wait

## SharedKeys vs Annotated Genome followed by variant filtering

for model in dn ch ibd Xl  ; do
    (
    grep -Fwf $familyID.sharedKeys.$model.txt $probandID.parentKeys.GenomeAnnotationsCombined.$model.txt > $probandID.sharedKeys.GenomeAnnotationsCombined.$model.txt # Extracting data/keys that has same pattern in both sharedKey (proband and sib) and original parent genome
    awk '{OFS="\t"; print $(NF-3), $(NF-2)}' $sibID.parentKeys.GenomeAnnotationsCombined.$model.txt | grep -Fwf $familyID.sharedKeys.$model.txt - | sort -k2 > $sibID.sharedKeys.$model.txt  #generating 4th and 5th column from right of sharedKeys into newly generated data
    awk '{OFS="\t"; print $(NF-2), NR}' $probandID.sharedKeys.GenomeAnnotationsCombined.$model.txt | sort -k1,1 > $familyID.sharedKeys.numbered.$model.txt
    join -a 1 -e "." -o '0,2.1,1.2' -1 1 -2 2 $familyID.sharedKeys.numbered.$model.txt $sibID.sharedKeys.$model.txt | sort -k3 -n | cut -f2 -d " " > $sibID.sharedKeys.column.$model.txt #combine both new and original shared data to generate a third new file
    paste $probandID.sharedKeys.GenomeAnnotationsCombined.$model.txt $sibID.sharedKeys.column.$model.txt > $familyID.sharedKeys.GenomeAnnotationsCombined.$model.txt

    ##variant filtering
    head -n1 $familyID.sharedKeys.GenomeAnnotationsCombined.$model.txt > $familyID.BestGeneCandidates.$model.txt #header
    awk -F"\t" '$13 < 0.005 && $14 < 0.001 && $15 < 0.005 && $16 < 0.005 && $18 < 0.01 && $20 < 0.001 && $21 < 0.001 && $22 < 0.001 && $(NF-7)=="PASS" {print}' $familyID.sharedKeys.GenomeAnnotationsCombined.$model.txt |\
    grep -Fwvf ~/Documents/Scripts/local/Func.Gene.Remove.txt - >> $familyID.BestGeneCandidates.$model.txt # Variant frequency filters, exonic, splicing and PASS
    ) &
done
wait

# Clean up
rm *.numbered.* *.column.*
