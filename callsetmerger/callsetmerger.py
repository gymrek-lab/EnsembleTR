#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
python callsetmerger.py --vcfs hipstr.chr21.sorted.vcf.gz,advntr.chr21.sorted.vcf.gz,gangstr.chr21.sorted.vcf.gz
"""
import argparse
import vcf
import trtools.utils.tr_harmonizer as trh
import trtools.utils.mergeutils as mergeutils
import trtools.utils.common as common

class VCFWrapper:
    def __init__(self, reader, vcftype):
        self.vcftype = vcftype
        self.vcfreader = reader

class Readers:
    def __init__(self, vcfpaths):
        self.vcfwrappers = []
        # Load all VCFs, make sure we can infer type
        for invcf in vcfpaths:
            vcffile = vcf.Reader(open(invcf, "rb"))
            hm = trh.TRRecordHarmonizer(vcffile)
            self.vcfwrappers.append(VCFWrapper(vcffile, hm.vcftype))
        # Get chroms
        contigs = self.vcfwrappers[0].vcfreader.contigs
        self.chroms = list(contigs)
        self.current_records = [next(wrapper.vcfreader) for wrapper in self.vcfwrappers]
        self.done = mergeutils.DoneReading(self.current_records)
        if not self.AreChromsValid():
            return 1
        self.is_min_pos_list = mergeutils.GetMinRecords(self.current_records, self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos =\
            self.GetCurrentRange()
        

    def AreChromsValid(self):
        for r, wrapper in zip(self.current_records, self.vcfwrappers):
            if r is None: continue
            if not r.CHROM in self.chroms:
                common.WARNING((
                "Error: found a record in file {} with "
                "chromosome '{}' which was not found in the contig list "
                "({})").format(wrapper.vcfreader.filename, r.CHROM, ", ".join(self.chroms)))
                return False
        return True

    def GetIsMinPosList(self):
        return mergeutils.GetMinRecords(self.current_records, self.chroms)

    def GetCurrentRange(self):
        start_pos = None
        end_pos = -1
        chrom = None
        for i in range(len(self.is_min_pos_list)):
            if self.is_min_pos_list[i]:
                chrom = self.current_records[i].CHROM
                start_pos = self.current_records[i].POS
                end = start_pos + len(trh.HarmonizeRecord(self.vcfwrappers[i].vcftype, \
                    self.current_records[i]).ref_allele)
                if end > end_pos: end_pos = end
        return chrom, start_pos, end_pos


    def GetOverlapMinRecords(self):
        # TODO move to increment?
        # Set is_overlap_min for anything overlapping
        is_overlap_min = []
        for i in range(len(self.current_records)):
            if self.current_records[i] is None:
                is_overlap_min.append(False)
                continue
            if self.current_records[i].CHROM != self.cur_range_chrom:
                is_overlap_min.append(False)
                continue
            start_pos = self.current_records[i].POS
            end = start_pos + len(trh.HarmonizeRecord(self.vcfwrappers[i].vcftype, \
                self.current_records[i]).ref_allele)        
            if end <= self.cur_range_end_pos:
                is_overlap_min.append(True)
            else: is_overlap_min.append(False)

        return is_overlap_min

    def MergeCalls(self, is_overlap_min):
        if sum(is_overlap_min) <= 1:
            return None # TODO. boring for debugging. just a single genotyper

        print("############")
        print("%s:%s-%s"%(self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos))
        allele_list = []
        sample_calls = {}
        # Print out info
        for i in range(len(self.current_records)):        
            if is_overlap_min[i] and self.current_records[i] is not None:
                hm = trh.HarmonizeRecord(self.vcfwrappers[i].vcftype, self.current_records[i])
                ref = hm.ref_allele
                repunit = hm.motif
                allele_list.append((ref, "*", self.vcfwrappers[i].vcftype, "GT:0", repunit))
                altnum = 1
                for alt in hm.alt_alleles:
                    allele_list.append((alt, len(alt)-len(ref), self.vcfwrappers[i].vcftype, "GT:%s"%altnum, repunit))
                    altnum += 1
                for sample in self.current_records[i]:
                    if sample.sample not in sample_calls: sample_calls[sample.sample] = []
                    if sample is None or not sample.called: 
                        sample_calls[sample.sample].append(([None, None], self.vcfwrappers[i].vcftype))
                        continue
                    sample_calls[sample.sample].append((sample["GT"], self.vcfwrappers[i].vcftype))
        if len(allele_list) == sum(is_overlap_min): return None # TODO boring for debugging. all ref

        for a in allele_list:
            print(a)
        for s in sample_calls:
            print(s)
            for call in sample_calls[s]:
                print("%s: %s"%(call[1], call[0]))

    def GoToNext(self):
        self.current_records = mergeutils.GetNextRecords([w.vcfreader for w in self.vcfwrappers], self.current_records, self.GetIsMinPosList())
        # after increment, check if chroms valid
        if not self.AreChromsValid():
            return 1
        self.done = mergeutils.DoneReading(self.current_records)
        self.is_min_pos_list = mergeutils.GetMinRecords(self.current_records, self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos =\
            self.GetCurrentRange()

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcfs", help="Comma-separated list of VCFs to merge. Must be sorted/indexed", type=str, required=True)
    args = parser.parse_args()

    readers = Readers(args.vcfs.split(","))


    # Check samples same in each VCF
    # TODO

    # Walk through sorted readers
    while not readers.done:

        # Determine overlap with min record, candidates for merging
        is_overlap_min = readers.GetOverlapMinRecords()

        # TODO merge
        was_merged = readers.MergeCalls(is_overlap_min)

        # TODO will need to update is_min_pos_list if we merged something

        # Move on
        readers.GoToNext()

if __name__ == "__main__":
    main()