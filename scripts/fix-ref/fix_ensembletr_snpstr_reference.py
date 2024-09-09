#!/usr/bin/env python

"""
Script to clean up SNP/TR reference haplotype panels
"""

import argparse
import cyvcf2
import pyfaidx
import sys
import subprocess
import os
import tempfile

MAX_ALLOWED_DUPS = 1000

def GetWriter(fname, reader):
    """
    Get regular file writer, not cyvcf2.Writer
    Since we are going to update INFO/FORMAT downstream
    """
    if fname == "stdout.vcf":
        outf = sys.stdout
    else: outf = open(fname, "w")
    # Get header from reader, but remove DS/GP
    header = reader.raw_header
    for line in header.split("\n"):
        if line.startswith("##FORMAT=<ID=DS,") or \
            line.startswith("##FORMAT=<ID=GP,") or \
            line.strip() == "":
            continue
        outf.write(line.strip() + "\n")
    return outf

def GetAlleleCounts(record):
    """
    Returns:
    - allele_counts: dictionary of allele_string:allele_count
    - allele_order: list of alleles, in same order as original 
      but with those with AC=0 removed
    """
    alleles = [record.REF] + record.ALT
    all_allele_calls = []
    for genotype in record.genotypes:
        all_allele_calls.append(alleles[genotype[0]])
        all_allele_calls.append(alleles[genotype[1]])
    allele_counts = {}
    for a in alleles:
        acount = all_allele_calls.count(a)
        if acount > 0: allele_counts[a] = acount
    # Keep order same as original
    allele_order = [item for item in alleles if item in allele_counts.keys()]
    return allele_counts, allele_order

def WriteRecord(writer, record):
    """simple function to write the record as is """
    writer.write(str(record).strip()+"\n")

def UpdateINFO(INFO, allele_counts, allele_order):
    """
    Keep INFO as is, except update AC/AF
    """
    info_items = []
    for infokey in dict(INFO).keys():
        if infokey in ["AC", "AF"]:
            continue
        else:
            # Note all other INFO fields are single values
            infoval = INFO.get(infokey)
            info_items.append(f"{infokey}={infoval}")
    # Add AC/AF
    num_chroms = sum(allele_counts.values())
    ac_vals = [allele_counts[a] for a in allele_order[1:]]
    af_vals = [allele_counts[a]/num_chroms for a in allele_order[1:]]
    info_items.append("AC=" + ",".join([str(item) for item in ac_vals]))
    info_items.append("AF=" + ",".join(["%.4f"%item for item in af_vals]))
    return ";".join(info_items)

def GetGT(sample, orig_alleles, allele_order):
    """Update the GT based on new allele order"""
    if sample[2] == True:
        sep = "|"
    else: sep = "/"
    a1 = allele_order.index(orig_alleles[sample[0]])
    a2 = allele_order.index(orig_alleles[sample[1]])
    return f"{a1}{sep}{a2}"

def IsTRRecord(record_id):
    """Check if the record is a TR"""
    return record_id is None or record_id.strip() == "."

def GetTRRecordID(record, allids):
    """Get new uniquified record ID"""
    locid = "EnsTR:{chrom}:{pos}".format(chrom=record.CHROM, pos=record.POS)
    if locid not in allids:
        return locid
    for i in range(1, MAX_ALLOWED_DUPS):
        newlocid = f"{locid}:{i}"
        if newlocid not in allids:
            sys.stderr.write(f"Adding duplicate locus {newlocid}\n")
            return newlocid
    raise ValueError(f"Error: too many duplicates of {locid}")

def CheckReference(record, refgenome):
    """ REF in VCF should match what is in the reference genome"""
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
    parser.add_argument("--out", help="Prefix for output VCF file. or 'stdout' for standard output", type=str, default="stdout")
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
    writer = GetWriter(args.out + ".vcf", reader)

    # Go through each record
    # If SNP: just print it
    # if STR: (1) Filter records with incorrect reference
    #         (2) Filter too many/too few alleles
    #         (3) Modify record ID
    #         (4) Write record, without DP/GS, with updated AF/AC and without AC=0
    allids = set()
    num_records_processed = 0
    num_snps_keeped = 0
    num_strs_keeped = 0
    num_str_failed_ref = 0
    num_str_failed_max_allele = 0
    num_str_failed_min_allele = 0

    for record in reader:
        num_records_processed += 1
        if args.max_records > 0 and num_records_processed > args.max_records:
            break # for debug
        if not IsTRRecord(record.ID):
            record.INFO["VT"] = "OTHER"
            WriteRecord(writer, record)
            num_snps_keeped += 1
        else:
    		# (1) Filter records with incorrect reference
            if not CheckReference(record, refgenome):
                num_str_failed_ref += 1
                continue # skip this record

            # (2) Filter too many/too few alleles
            # Note: GetAlleleCounts() doesn't include things with AC=0
            allele_counts, allele_order = GetAlleleCounts(record)
            num_alleles = len(allele_counts.keys())
            if args.max_alleles != -1 and num_alleles > args.max_alleles:
                num_str_failed_max_allele += 1
                sys.stderr.write(f"Skipping {record.CHROM}:{record.POS} with {num_alleles} alleles\n")
                continue # skip this record
            if args.min_alleles != -1 and num_alleles < args.min_alleles:
                num_str_failed_min_allele += 1
                sys.stderr.write(f"Skipping {record.CHROM}:{record.POS} with {num_alleles} alleles\n")
                continue # skip this record

            # (3) Modify record ID
            record.ID = GetTRRecordID(record, allids)
            allids.add(record.ID)
            record.INFO["VT"] = "TR"
            num_strs_keeped += 1
            # (4) Write record to file, update AC/AF, only include GT, exclude AC=0
            orig_alleles = [record.REF] + record.ALT
            updated_info = UpdateINFO(record.INFO, allele_counts, allele_order)
            out_items = [record.CHROM, record.POS, record.ID, \
                allele_order[0], ",".join(allele_order[1:]), \
                ".", "PASS", updated_info, "GT"]
            for sample in record.genotypes:
                out_items.append(GetGT(sample, orig_alleles, allele_order))
            writer.write("\t".join([str(item) for item in out_items])+"\n")
    sys.stdout.write(f"All {num_records_processed:,} records processed: keep {num_snps_keeped:,} snps, {num_strs_keeped:,} STRs; remove {num_str_failed_ref:,} STRs mismatch with reference, {num_str_failed_max_allele:,} STRs with more than {args.max_alleles:,} alleles, {num_str_failed_min_allele:,} STRs less than {args.min_alleles} alleles.\n")
    reader.close()
    writer.close()
    sys.exit(0)

if __name__ == "__main__":
    run()
