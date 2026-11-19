#!/usr/bin/python

'''

SA Cancer Genomics Facility
Created on Oct 18, 2012

@author: Dave Lawrence

Edits:
11/07/2018: Mark Corbett: Fix for library names that have underscores
'''

import csv
import os
import re
import sys

REGEX_VARIANT = r"\S+\s+\S+\s+\S+\s+\S+\s+\S+"

def readVariants(fileName, extension):
    pattern = f"^{extension}\t([^\t]+)\t({REGEX_VARIANT})"

    variants = {}
    with open(fileName) as f:
        for line in f:
            match = re.search(pattern, line)
            if match is None:
                raise ValueError(f"'{line}' has no match for '{pattern}'")

            variants[match.group(2)] = match.group(1)

    return variants


def loadColumnsForVariants(annovarOutputFiles):
    data = {}
    for fileName in annovarOutputFiles:
        extension = os.path.splitext(fileName)[1][1:]
        if extension.endswith(("log", "_filtered")):
            continue

        try:
            components = extension.split("_")
            if len(components) < 2:
                raise ValueError("Expected extension ('%s') to have 2 underscore separated parts")
            annovarDB = components[1]
            if len(components) > 2:
                annovarDB = '_'.join(components[1:-1])
            data[annovarDB] = readVariants(fileName, annovarDB)
        except Exception as e:
            sys.stderr.write(f"Problem parsing file: '{fileName}'\n")
            raise e
    return data

def parseAnnovarInputLine(line):
    pattern = rf"^({REGEX_VARIANT})\s+(.*)"
    match = re.search(pattern, line)
    if not match:
        raise ValueError(f"line '{line}' did not match '{pattern}'")

    variant = match.group(1)
    other = match.group(2)
    startColumns = variant.split('\t')
    otherColumns = other.split('\t')

    return (variant, startColumns, otherColumns)


def combineAnnovarOutputs(annovarInput, annovarOutputFiles):
    data = loadColumnsForVariants(annovarOutputFiles)
    dbColumnNames = sorted(data.keys())
    header = ["chr", "start", "end", "ref", "obs"] + dbColumnNames + ["other..."]

    writer = csv.writer(sys.stdout)
    writer.writerow(header)

    with open(annovarInput) as f:
        for line in f:
            (variant, startColumns, otherColumns) = parseAnnovarInputLine(line)

            row = startColumns
            for c in dbColumnNames:
                row.append(data[c].get(variant, ""))

            row += otherColumns
            writer.writerow(row)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.stderr.write(f"Usage {os.path.basename(sys.argv[0])} avinput <list of files to combine>\n")
        sys.exit(1)
    
    combineAnnovarOutputs(sys.argv[1], sys.argv[2:])

    
