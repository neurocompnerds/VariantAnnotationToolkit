#!/usr/bin/python3

# This script fetches gene panel information from several different sources and then updates the file in ANNOVAR.

import sys, getopt, csv, os
import requests
import json

def usage():
    print(
'''
# update.gene.panels.py a script to extract gene panels from specified sources for inclusion in our ANNOVAR database.
# This utility will fetch gene panel information from Genomics England or Australian PanelApp and the Genes4Epilepsy GitHub page
# information is converted to bed format or other formats used by ANNOVAR.
#
# Usage update.gene.panels.py -p "ID | ASD | Epilepsy | DDG2P | CP" | "all" | [ -l | --list ] | [ -h | --help ]
#
# Options:
# -p             OPTIONAL: A list of gene panel IDs to fetch. Separate multiple IDs with a comma. If "all" is specified, then all panels will be fetched.
# -l | --list    OPTIONAL: Lists all available gene panel IDs.
# -h | --help    OPTIONAL: Displays usage information.
#
# Script created by Mark Corbett with AI assistance on 13/11/2026 
# Contact: mark.corbett at adelaide.edu dot au
# Edit History (Name; Date; Description)
'''
         )


