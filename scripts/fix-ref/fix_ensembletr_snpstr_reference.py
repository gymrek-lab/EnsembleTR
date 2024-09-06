#!/usr/bin/env python

"""
Script to clean up SNP/TR reference haplotype panels

TODO: remove alleles with AC=0
"""

import argparse
import cyvcf2
import pyfaidx
import sys
import subprocess
import os


MAX_ALLOWED_DUPS = 1000

def IsTRRecord(record_id):
    return record_id is None or record_id.strip() == "."

def GetTRRecordID(record, allids):
    locid = "EnsTR:{chrom}:{pos}"%.format(chrom=record.CHROM, pos=record.POS)
    if locid not in allids:
        return locid
    for i in range(1, MAX_ALLOWED_DUPS):
        newlocid = "{locid}:{i}"
        if newlocid not in allids:
            sys.stderr.write("Adding duplicate locus {newlocid}\n")
            return newlocid
    raise ValueError("Error: too many duplicates of {locid}")

def CheckReference(record, refgenome):
    # REF in VCF should match what is in the reference genome
    refseq = refgenome[record.CHROM][record.POS-1:record.POS-1+len(record.REF)]
    return (refseq == record.REF)

def run():
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Input SNP/TR VCF file", type=str, required=True)
    parser.add_argument("--ref", help="Reference fasta file", type=str, required=True)
    parser.add_argument("--out", help="Prefix for output VCF file", type=str, required=True)
    parser.add_argument("--max-alleles", help="Ignore loci with more than this many alleles", type=int, default=-1)
    parser.add_argument("--min-alleles", help="Ignore loci with fewer than this many alleles", type=int, default=-1)
    parser.add_argument("--max-records", help="Quit after processing this many records (for debug)", \
        default=-1, type=int)
    args = parser.parse_args()
    return args

def main(args):
    # Set up VCF reader
    reader = cyvcf2.Reader(args.vcf)

    # Set up the reference
    refgenome = pyfaidx.Fasta(args.ref)

    # Set up writer, adding the missing header
    reader.add_to_header('##command=hipstr;note this is a dummy header line')
    reader.add_to_header('##INFO=<ID=VT,Number=1,Type=String,Description="Type of variant">')
    writer = cyvcf2.Writer(args.out + ".vcf", reader)

    # Go through each record
    # If SNP: just print it
    # if STR: filter records with incorrect reference,
    #   and modify ID
    allids = set()
    num_records_processed = 0
    for record in reader:
        num_records_processed += 1
        if args.max_records > 0 and num_records_processed > args.max_records:
            break # for debug
        if not IsTRRecord(record.ID):
            record.INFO["VT"] = "OTHER"
            writer.write_record(record)
        else:
    		# Check reference
            if not CheckReference(record, refgenome):
                continue # skip this record
            # TODO - remove alleles with AF=0
            # Check if too many alleles
            if args.max_alleles != -1 and (1+len(record.ALT)) > args.max_alleles:
                sys.stderr.write("Skipping {chrom}:{pos} with {numalt} ALT alleles\n".format(
                    chrom=record.CHROM, pos=record.POS, numalt=len(record.ALT))
                )
                continue 
            # Check if too few alleles
            # Include cases where there is an ALT listed but the AF is 0
            if (args.min_alleles != -1 and (1+len(record.ALT)) < args.min_alleles) or \
                (record.INFO["AF"]==0):
                sys.stderr.write("Skipping {chrom}:{pos} with {numalt} ALT alleles, AF={alleles}\n".format(
                    chrom=record.CHROM, pos=record.POS, numalt=len(record.ALT), alleles=str(record.INFO["AF"]))
                )
                continue
    		# Modify ID
            record.ID = GetTRRecordID(record, allids)
            allids.add(record.ID)
            record.INFO["VT"] =  "TR"
    		# Write to file
            writer.write_record(record)

    reader.close()
    writer.close()
    sys.exit(0)

if __name__ == "__main__":
    run()