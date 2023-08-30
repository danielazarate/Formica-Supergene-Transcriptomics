#!/usr/bin/python3

# Python script for extracting specific chromosomes (or scaffolds or contigs...) \
# and their alignment information from a SAM file. 

# Takes as input a SAM file (generated from BWA MEM2 and an argument (list of \
# names of regions one wishes to extract. We also have a neat argument (--invert) \
# which allows you to extract everything BUT what is specified. 
# Outputs: a file with the rows extracted from the SAM ".scaffExtract"

# Script used in publication - \
# Zarate et al. In Prep. Differential Gene Expression of a supergene regulating \
# social behavior in an ant (Formica neoclara). 

# Authors: Daniela Zarate, PHD UC Riverside & Stephen Brennan, UC San Diego 

import os
import argparse

# Set up a very helpful argument parser ecosystem to read in arguments passed through on the command line
parser = argparse.ArgumentParser(
    description='script for filtering out individuals with high proportions of missing data (VCF files)',
    epilog="Good Luck with Your Bioinformatics! This script is written and directed by Daniela Zarate PHD")
parser.add_argument('--debug', help='print debug output', action='store_true')
parser.add_argument('--input', help='file to read from')
parser.add_argument('--scaffs', help='list of scaffolds to extract/exclude')
parser.add_argument('--invert', help='extract all but named scaffold (y/n)', default='n')


## Global Variables

args = parser.parse_args()
output = args.input + ".scaffExtract"
scaffList = args.scaffs.split()

## Main

# Read input file, write header and all desired lines to file
with open(args.input, 'r') as file, open(output, 'w+') as outFile:
    for line in file:
        if args.invert == 'y':
            if line.split()[0] == "@SQ" or line.split()[0] == "@PG" or line.split()[2] not in scaffList:
                outFile.write(line)
        else:    
            if line.split()[0] == "@SQ" or line.split()[0] == "@PG" or line.split()[2] in scaffList:
                outFile.write(line)
    
