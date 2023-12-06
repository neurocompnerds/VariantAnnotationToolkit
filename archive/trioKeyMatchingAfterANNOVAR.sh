#!/bin/bash

usage()
{
echo "# $0 A script to append keys to best candidate files for manual filtering
# Requirements: The script needs to be executed from the folder that contains the BestGeneCandidates file of the proband
# The script expects to find parent files in folders with the same name as their ID side by side with the folder that contains the proband's files
# Most of our ANNOVAR scripts will set this up so you don't have to worry about it.
#
# Usage $0 -p probandID -m motherID -f fatherID | [-h | --help]
#
#	-p 			ID of proband
#	-m			ID of mother
#	-f			ID of father
#	-h --help	Prints this message
#
# Example:
# $0 -p 001P -m 001M -f 001F
#
# History:
# Original: Mark Corbett 08/02/2016
# mark.corbett@adelaide.edu.au
# Modified (date; name; description)
# 30/03/2016; Mark Corbett; Add in inheritance model filters
# 02/08/2016; Mark Corbett; Make file recognition patterns a little more permissive
# 06/03/2019; Mark Corbett; Found error with adding de novo variants to the CH gene list was only capturing duplicates not all unique gene hits
#
"
}

while [ "$1" != "" ]; do
	case $1 in
		-p )		shift
					probandID=$1
					;;
		-m )		shift
					mumID=$1
					;;
		-f )		shift
					dadID=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done

if [ -z "$probandID" ] ; then
	usage
	echo "You need to specify a proband ID with -p \"probandID\""
	exit 1
fi
if [ -z "$mumID" ] ; then
	usage
	echo "You need to specify mum's ID with -m \"motherID\""
	exit 1
fi
if [ -z "$dadID" ] ; then
	usage
	echo "You need to specify dad's ID with -f \"fatherID\""
	exit
fi

workDir=$(pwd)

## Start the script
if [ -f $probandID.BestGeneCandidates.txt ] ; then
	awk '{print $NF}' $probandID.BestGeneCandidates.txt > meh.txt
	awk '{OFS="\t"; print $(NF-1), $NF}' ../$mumID/*.GenomeAnnotationsCombined.txt > $mumID.keys.txt &
	awk '{OFS="\t"; print $(NF-1), $NF}' ../$dadID/*.GenomeAnnotationsCombined.txt > $dadID.keys.txt
	wait
else
	usage
	echo "I couldn't find the file $probandID.BestGeneCandidates.txt.
	Sorry I'm a bit particular about file names and I need to be executed from the folder that contains the proband's BestGeneCandidates file"
	exit 1
fi

grep -Fwf meh.txt $mumID.keys.txt | sort -k2 > $mumID.BC.keys.txt &
grep -Fwf meh.txt $dadID.keys.txt | sort -k2 > $dadID.BC.keys.txt
wait

awk '{OFS="\t"; print $0, NR}' meh.txt | sort -k1,1 > $probandID.BC.keys.txt

join -a 1 -e "." -o '0,2.1,1.2' -1 1 -2 2 $probandID.BC.keys.txt $mumID.BC.keys.txt | sort -k3 -n | cut -f2 -d " " > $mumID.column.txt &
join -a 1 -e "." -o '0,2.1,1.2' -1 1 -2 2 $probandID.BC.keys.txt $dadID.BC.keys.txt | sort -k3 -n | cut -f2 -d " " > $dadID.column.txt
wait

paste $probandID.BestGeneCandidates.txt $mumID.column.txt $dadID.column.txt > $probandID.parentKeys.BestGeneCandidates.txt
head -n1 $probandID.parentKeys.BestGeneCandidates.txt > $probandID.parentKeys.BestGeneCandidates.header.txt

# Filter based on inheritance models
cp $probandID.parentKeys.BestGeneCandidates.header.txt $probandID.parentKeys.BestGeneCandidates.dn.txt
cp $probandID.parentKeys.BestGeneCandidates.header.txt $probandID.parentKeys.BestGeneCandidates.ch.txt
cp $probandID.parentKeys.BestGeneCandidates.header.txt $probandID.parentKeys.BestGeneCandidates.ibd.txt
cp $probandID.parentKeys.BestGeneCandidates.header.txt $probandID.parentKeys.BestGeneCandidates.Xl.txt

awk '$(NF-3) ~ /^0\/1/ && $(NF-1) == "." && $NF == "." {OFS="\t"; print}' $probandID.parentKeys.BestGeneCandidates.txt >> $probandID.parentKeys.BestGeneCandidates.dn.txt
cut -f7 $probandID.parentKeys.BestGeneCandidates.dn.txt | sort | uniq > $probandID.CH.gene.list.1.txt
awk '$1 == "chrX" && $(NF-3) ~ /^1\/1/ && $(NF-1) !~ /^1/ && $NF == "." {OFS="\t"; print}' $probandID.parentKeys.BestGeneCandidates.txt >> $probandID.parentKeys.BestGeneCandidates.Xl.txt
awk '$(NF-3) ~ /^1\/1/ && $(NF-1) !~ /^1/ && $NF !~ /^1/ {OFS="\t"; print}' $probandID.parentKeys.BestGeneCandidates.txt >> $probandID.parentKeys.BestGeneCandidates.ibd.txt

awk '$(NF-3) ~ /^0\/1/ && $(NF-1) ~ /^0\/1/ && $NF == "." {OFS="\t"; print}' $probandID.parentKeys.BestGeneCandidates.txt > $probandID.parentKeys.BestGeneCandidates.chM.txt
cut -f7 $probandID.parentKeys.BestGeneCandidates.chM.txt | sort | uniq >> $probandID.CH.gene.list.1.txt
awk '$(NF-3) ~ /^0\/1/ && $(NF-1) == "." && $NF ~ /^0\/1/ {OFS="\t"; print}' $probandID.parentKeys.BestGeneCandidates.txt > $probandID.parentKeys.BestGeneCandidates.chF.txt
cut -f7 $probandID.parentKeys.BestGeneCandidates.chF.txt | sort | uniq >> $probandID.CH.gene.list.1.txt
grep -v Gene.gene $probandID.CH.gene.list.1.txt | sort | uniq -d > $probandID.CH.gene.list.2.txt
cat $probandID.parentKeys.BestGeneCandidates.dn.txt $probandID.parentKeys.BestGeneCandidates.chM.txt $probandID.parentKeys.BestGeneCandidates.chF.txt | sed '1d' | grep -wf $probandID.CH.gene.list.2.txt >> $probandID.parentKeys.BestGeneCandidates.ch.txt

rm $mumID.column.txt $dadID.column.txt meh.txt $mumID.BC.keys.txt $dadID.BC.keys.txt $mumID.keys.txt $dadID.keys.txt $probandID.parentKeys.BestGeneCandidates.header.txt $probandID.parentKeys.BestGeneCandidates.chM.txt $probandID.parentKeys.BestGeneCandidates.chF.txt $probandID.CH.gene.list.1.txt $probandID.CH.gene.list.2.txt
