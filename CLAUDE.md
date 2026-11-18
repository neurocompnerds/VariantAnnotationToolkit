# The Variant Annotation Toolkit Repository

## Background:
This repository was initially just a collection of scripts in a folder not built within GitHub from the start.  It is not perfect by any means but it "works" so care needs to be taken before making edits.  The repository contains scripts that coordinate annotation of variant call format (VCF) files using the ANNOVAR program (https://annovar.openbioinformatics.org/en/latest/).  The subject matter is human medical genetics and genomics.

Typical usage of the repository is to annotate a VCF file using the `ANNOVARv3_for_hg38.sh` script and then downstream filtering of the output using either the `trioKeyMatchingAfterANNOVAR_hg38.py` or `familyKeyMatchingAfterANNOVAR_hg38.py` scripts.
Note that the `ANNOVARv3_for_hg38.sh` script calls `annovar_combine_csv.py` and the `AnnovarGenomeSummaryCombo.v2_for_hg38.sh` scripts.
Read these scripts and familiarise yourself with the key activities of the repository.

## Languages:
There are a mix of bash, python and perl scripts here.  The prefered language is python for most operations that require manipulation of data tables and bash for coordinating the annotation activity.

## Rules:
Show your plans before making any edits.
You should expect and request a written work plan from the user for new work.
Do not push changes to github.
Ignore the `archive` folder and never make any edits there unless specifically told to.
