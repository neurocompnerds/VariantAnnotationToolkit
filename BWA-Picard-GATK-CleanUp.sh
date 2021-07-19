#!/bin/bash
# This is an automatically generated script that can be used to remove redundant files from your current directory and link your reduced read bam to the bam file library
# To use this you can place the following lines in your script:
# cp /path/to/this/header/BWA-Picard-GATK-CleanUp.txt $OUTPREFIX.CleanUp.sh # on Neurogenetics:Brain the default location is /home/neuro/Documents/Pipelines/BWA-Picard-GATK-CleanUp.txt
# echo "rm /path/to/<file.to.delete>" >> $OUTPREFIX.CleanUp.sh # Use this pattern for as many files as you like
# For any file you might want to keep and delete later such as the SAM file you can use the following syntax
# echo "deleteLevel=1; if [ \$Keep < \$deleteLevel ] ; then ; rm <my.precious.file.name> ; fi"
# Header and template created by Mark Corbett on 08/05/2013
# Modified:
usage()
{
echo "# Usage $0 [-k --keep value]
#	-h --help	Prints this message
#	-k --keep [integer value]	Default=0 which means all files will be removed unless you set the value of deleteLevel greater than 0 in the following line:
#			echo \"deleteLevel=1 ; if [ \\\$Keep < \$deleteLevel ] ; then rm <my.precious.file.name> ; fi\"
# Example:
# $0 -k 1 
# For any file where the line echo \"deleteLevel=1 ; if [ \\\$Keep < \$deleteLevel ] ; then rm <my.precious.file.name> ; fi\" has been added will not be deleted"
}

Keep=0 # Default is true which will DELETE the sam file.  Yes I know this is counter intuitive but it is easier to script it this way

while [ "$1" != "" ]; do
	case $1 in
		-k | --keep )	shift
								Keep=$1
								;;
		-h | --help )			usage
								exit
								;;
		* )						usage
								exit 1
	esac
	shift
done

# If there is nothing below this line then this is just the header file and it doesn't really do much
