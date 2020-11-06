#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
"""
import argparse

import trtools.utils.common as common
import trtools.utils.mergeutils as mergeutils
import trtools.utils.tr_harmonizer as trh
import vcf



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
        self.current_tr_records = [trh.HarmonizeRecord(wrapper.vcftype, next(wrapper.vcfreader))
                                   for wrapper in self.vcfwrappers]
        self.done = all([item.vcfrecord is None for item in self.current_tr_records])
        if not self.areChromsValid():
            raise ValueError('Invalid CHROM detected in record.')
        self.is_min_pos_list = mergeutils.GetMinRecords([item.vcfrecord for item in self.current_tr_records], self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos = \
            self.getCurrentRange()
        self.is_overlap_min = self.getOverlapMinRecords()

    def areChromsValid(self):
        for r, wrapper in zip(self.current_tr_records, self.vcfwrappers):
            if r.vcfrecord is None:
                continue
            if r.vcfrecord.CHROM not in self.chroms:
                common.WARNING((
                                   "Error: found a record in file {} with "
                                   "chromosome '{}' which was not found in the contig list "
                                   "({})").format(wrapper.vcfreader.filename, r.vcfrecord.CHROM,
                                                  ", ".join(self.chroms)))
                return False
        return True

    def getIsMinPosList(self):
        return mergeutils.GetMinRecords([item.vcfrecord for item in self.current_tr_records], self.chroms)

    def getCurrentRange(self):
        start_pos = None
        end_pos = -1
        chrom = None
        for i in range(len(self.is_min_pos_list)):
            if self.is_min_pos_list[i]:
                chrom = self.current_tr_records[i].vcfrecord.CHROM
                start_pos = self.current_tr_records[i].vcfrecord.POS
                end = start_pos + len(trh.HarmonizeRecord(self.vcfwrappers[i].vcftype,
                                                          self.current_tr_records[i].vcfrecord).ref_allele)
                if end > end_pos:
                    end_pos = end
        return chrom, start_pos, end_pos

    def getOverlapMinRecords(self):
        # Set is_overlap_min for anything overlapping
        is_overlap_min = []
        for i in range(len(self.current_tr_records)):
            if self.current_tr_records[i] is None:
                is_overlap_min.append(False)
                continue
            if self.current_tr_records[i].vcfrecord.CHROM != self.cur_range_chrom:
                is_overlap_min.append(False)
                continue
            start_pos = self.current_tr_records[i].vcfrecord.POS
            end = start_pos + len(trh.HarmonizeRecord(self.vcfwrappers[i].vcftype,
                                                      self.current_tr_records[i].vcfrecord).ref_allele)
            if end <= self.cur_range_end_pos:
                is_overlap_min.append(True)
            else:
                is_overlap_min.append(False)

        return is_overlap_min

    def getMergableCalls(self):
        if sum(self.is_overlap_min) <= 1:
            return None  # TODO. boring for debugging. just a single genotyper

        print("############")
        print("%s:%s-%s" % (self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos))
        allele_list = []
        sample_calls = {}
        # Print out info
        for i in range(len(self.current_tr_records)):
            if self.is_overlap_min[i] and self.current_tr_records[i] is not None:
                hm = trh.HarmonizeRecord(self.vcfwrappers[i].vcftype, self.current_tr_records[i].vcfrecord)
                ref = hm.ref_allele
                repunit = hm.motif
                allele_list.append((ref, "*", self.vcfwrappers[i].vcftype, "GT:0", repunit))
                altnum = 1
                for alt in hm.alt_alleles:
                    allele_list.append(
                        (alt, len(alt) - len(ref), self.vcfwrappers[i].vcftype, "GT:%s" % altnum, repunit))
                    altnum += 1
                for sample in self.current_tr_records[i].vcfrecord:
                    if sample.sample not in sample_calls:
                        sample_calls[sample.sample] = []
                    if sample is None or not sample.called:
                        sample_calls[sample.sample].append(([None, None], self.vcfwrappers[i].vcftype))
                        continue
                    sample_calls[sample.sample].append((sample["GT"], self.vcfwrappers[i].vcftype))
        if len(allele_list) == sum(self.is_overlap_min):
            return None  # TODO boring for debugging. all ref

        for a in allele_list:
            print(a)
        for s in sample_calls:
            print(s)
            for call in sample_calls[s]:
                print("%s: %s" % (call[1], call[0]))
        return allele_list, sample_calls

    def getCurrentRecords(self):
        ret = []
        for rec in self.current_tr_records:
            ret.append(rec.vcfrecord)
        return ret

    def goToNext(self):
        prev_records = self.current_tr_records
        new_records = []
        for idx, rec in enumerate(prev_records):
            if self.getIsMinPosList()[idx]:
                try:
                    new_records.append(
                        trh.HarmonizeRecord(self.vcfwrappers[idx].vcftype,
                                            next(self.vcfwrappers[idx].vcfreader)))
                except StopIteration:
                    new_records.append(None)
            else:
                new_records.append(
                    trh.HarmonizeRecord(self.vcfwrappers[idx].vcftype,
                                        prev_records[idx].vcfrecord)
                )
        self.current_tr_records = new_records
        self.updateObject()

    def updateObject(self):
        if not self.areChromsValid():
            raise ValueError('Invalid CHROM detected in record.')
        self.done = all([item.vcfrecord is None for item in self.current_tr_records])
        self.is_min_pos_list = mergeutils.GetMinRecords(self.getCurrentRecords(), self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos = \
            self.getCurrentRange()
        self.is_overlap_min = self.getOverlapMinRecords()


class MergeObject:
    def __init__(self, allele_list, sample_calls):
        self.allele_list = allele_list
        self.sample_calls = sample_calls

    def merge(self):
        pass
