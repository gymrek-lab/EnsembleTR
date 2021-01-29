#!/usr/bin/env python3

"""
Usage:

statSTR --vcf VCF --afreq --out stdout | ./get_numalleles.py > output

Parses statSTR output to add a column with number of common alleles
This is a fix until I add this stat to statSTR
"""

import sys

MINFREQ = 0.01

def GetAlleleCount(afreqs):
    allele_list = afreqs.split(",")
    num_alleles = 0
    for item in allele_list:
        try:
            allele, freq = item.split(":")
        except: return 0
        if float(freq) >= MINFREQ: num_alleles += 1
    return num_alleles

line = sys.stdin.readline()
header_items = line.strip().split()
afreq_col = header_items.index("afreq")
sys.stdout.write("\t".join(header_items+["numalleles"])+"\n")

line = sys.stdin.readline()
while line.strip() != "":
    items = line.strip().split()
    count = GetAlleleCount(items[afreq_col])
    sys.stdout.write("\t".join(items+[str(count)])+"\n")
    line = sys.stdin.readline()



