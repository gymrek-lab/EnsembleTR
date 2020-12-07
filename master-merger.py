#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Test
./master-merger.py --vcfs hipstr.chr21.sorted.vcf.gz,advntr.chr21.sorted.vcf.gz,gangstr.chr21.sorted.vcf.gz
"""

import argparse
import vcf
import trtools.utils.tr_harmonizer as trh
import trtools.utils.mergeutils as mergeutils
import trtools.utils.common as common
import networkx

def MergeCalls(current_records, vcftypes, is_overlap_min, chrom, start_pos, end_pos):
    if sum(is_overlap_min) <= 1:
        return None # TODO. boring for debugging. just a single genotyper    

    print("############")
    print("%s:%s-%s"%(chrom, start_pos, end_pos))
    allele_list = []
    sample_calls = {}
    # Print out info
    for i in range(len(current_records)):        
        if is_overlap_min[i] and current_records[i] is not None:
            hm = trh.HarmonizeRecord(vcftypes[i], current_records[i])
            ref = hm.ref_allele
            repunit = hm.motif
            allele_list.append((ref, "*", vcftypes[i], "GT:0", repunit))
            altnum = 1
            for alt in hm.alt_alleles:
                allele_list.append((alt, len(alt)-len(ref), vcftypes[i], "GT:%s"%altnum, repunit))
                altnum += 1
            for sample in current_records[i]:
                if sample.sample not in sample_calls: sample_calls[sample.sample] = []
                if sample is None or not sample.called: 
                    sample_calls[sample.sample].append(([None, None], vcftypes[i]))
                    continue
                sample_calls[sample.sample].append((sample["GT"], vcftypes[i]))
    if len(allele_list) == sum(is_overlap_min): return None # TODO boring for debugging. all ref

    for a in allele_list:
        print("%s:%s\t%s (%s) repunit=%s"%(str(a[2]).split(".")[1],a[3], a[0], a[1], a[4]))
    for s in sample_calls:
        print(s)
        for call in sample_calls[s]:
            print("%s: %s"%(call[1], call[0]))

def GetOverlapMinRecords(current_records, vcftypes, is_min):
    # Get range of records
    start_pos = None
    end_pos = -1
    chrom = None
    for i in range(len(is_min)):
        if is_min[i]:
            chrom = current_records[i].CHROM
            start_pos = current_records[i].POS
            end = start_pos + len(trh.HarmonizeRecord(vcftypes[i], current_records[i]).ref_allele)
            if end > end_pos: end_pos = end

    # Set is_overlap_min for anything overlapping
    is_overlap_min = []
    for i in range(len(current_records)):
        if current_records[i] is None:
            is_overlap_min.append(False)
            continue
        if current_records[i].CHROM != chrom:
            is_overlap_min.append(False)
            continue
        start_pos = current_records[i].POS
        end = start_pos + len(trh.HarmonizeRecord(vcftypes[i], current_records[i]).ref_allele)        
        if end <= end_pos:
            is_overlap_min.append(True)
        else: is_overlap_min.append(False)

    return is_overlap_min, chrom, start_pos, end_pos

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcfs", help="Comma-separated list of VCFs to merge. Must be sorted/indexed", type=str, required=True)
    args = parser.parse_args()

    vcfreaders = []
    vcftypes = []

    # Load all VCFs, make sure we can infer type
    for invcf in args.vcfs.split(","):
        vcffile = vcf.Reader(open(invcf, "rb"))
        hm = trh.TRRecordHarmonizer(vcffile)
        vcfreaders.append(vcffile)
        vcftypes.append(hm.vcftype)

    # Get chroms
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)

    # Check samples same in each VCF
    # TODO

    # Walk through sorted readers
    current_records = [next(reader) for reader in vcfreaders]
    done = mergeutils.DoneReading(current_records)
    while not done:
        # Check chrom in contigs
        for r, reader in zip(current_records, vcfreaders):
            if r is None: continue
            if not r.CHROM in chroms:
                common.WARNING((
                    "Error: found a record in file {} with "
                    "chromosome '{}' which was not found in the contig list "
                    "({})").format(reader.filename, r.CHROM, ", ".join(chroms)))
                return 1
        # Get "is min"
        is_min = mergeutils.GetMinRecords(current_records, chroms)

        # Determine overlap with min record, candidates for merging
        is_overlap_min, chrom, start_pos, end_pos = GetOverlapMinRecords(current_records, vcftypes, is_min)

        # TODO merge
        was_merged = MergeCalls(current_records, vcftypes, is_overlap_min, chrom, start_pos, end_pos)

        # TODO will need to update is_min if we merged something

        # Move on
        current_records = mergeutils.GetNextRecords(vcfreaders, current_records, is_min)
        done = mergeutils.DoneReading(current_records)

    
if __name__ == "__main__":
    main()
