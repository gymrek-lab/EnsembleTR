#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
"""
import argparse

import trtools.utils.common as common
import trtools.utils.utils as utils
import trtools.utils.mergeutils as mergeutils
import trtools.utils.tr_harmonizer as trh
import cyvcf2
from callsetmerger.recordcluster import RecordObj, RecordCluster, OverlappingRegion


class VCFWrapper:
    def __init__(self, reader, vcftype):
        self.vcfreader = reader
        self.vcftype = vcftype


class Readers:
    def __init__(self, vcfpaths):
        self.vcfwrappers = []
        # Load all VCFs, make sure we can infer type
        for invcf in vcfpaths:
            vcffile = cyvcf2.VCF(invcf)
            hm = trh.TRRecordHarmonizer(vcffile)
            self.vcfwrappers.append(VCFWrapper(vcffile, hm.vcftype))
        # Get chroms
        self.chroms = utils.GetContigs(self.vcfwrappers[0].vcfreader)
        self.current_tr_records = [trh.HarmonizeRecord(wrapper.vcftype, next(wrapper.vcfreader))
                                   for wrapper in self.vcfwrappers]
        self.done = all([item.vcfrecord is None for item in self.current_tr_records])
        if not self.areChromsValid():
            raise ValueError('Invalid CHROM detected in record.')
        self.is_min_pos_list = mergeutils.GetMinRecords(self.getCurrentRecordVCFRecs(),
                                                        self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos = \
            self.getCurrentRange()
        self.is_overlap_min = self.getOverlapMinRecords()

    def areChromsValid(self):
        for r, wrapper in zip(self.current_tr_records, self.vcfwrappers):
            if r is None or r.vcfrecord is None:
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
        return mergeutils.GetMinRecords(self.getCurrentRecordVCFRecs(), self.chroms)

    def getCurrentRecordVCFRecs(self):
        ret = []
        for item in self.current_tr_records:
            if item is None or item.vcfrecord is None:
                ret.append(None)
            else:
                ret.append(item.vcfrecord)
        return ret

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

            # a = self.current_tr_records[i].vcfrecord
            if self.current_tr_records[i] is None:
                is_overlap_min.append(False)
                continue
            if self.current_tr_records[i].vcfrecord.CHROM != self.cur_range_chrom:
                is_overlap_min.append(False)
                continue
            # TODO Discuss
            start_pos = self.current_tr_records[i].vcfrecord.POS
            # end = start_pos + len(trh.HarmonizeRecord(self.vcfwrappers[i].vcftype,
            #                                           self.current_tr_records[i].vcfrecord).ref_allele)
            # if end <= self.cur_range_end_pos:
            #     is_overlap_min.append(True)
            # else:
            #     is_overlap_min.append(False)

            end = start_pos + len(self.current_tr_records[i].vcfrecord.REF)
            if (self.cur_range_start_pos <= start_pos <= self.cur_range_end_pos) or \
                    (self.cur_range_start_pos <= end <= self.cur_range_end_pos):
                is_overlap_min.append(True)
            else:
                is_overlap_min.append(False)
        return is_overlap_min

    def getMergableCalls(self):
        # if sum(self.is_overlap_min) <= 1:
        #     return None  # TODO. boring for debugging. just a single genotyper

        # print("############")
        # print("%s:%s-%s" % (self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos))
        record_cluster_list = []
        # TODO remove (moving to another class)
        # allele_list = []
        # sample_calls = {}

        # Print out info
        for i in range(len(self.current_tr_records)):
            if self.is_overlap_min[i] and self.current_tr_records[i] is not None:
                curr_ro = RecordObj(self.vcfwrappers[i].vcftype, self.current_tr_records[i].vcfrecord)
                added = False
                for rc in record_cluster_list:
                    if rc.canonical_motif == curr_ro.canonical_motif:
                        rc.AppendRecordObject(curr_ro)
                        added = True
                if not added:
                    record_cluster_list.append(RecordCluster([curr_ro]))
        ov_region = OverlappingRegion(record_cluster_list)
        return ov_region

    def getCurrentRecords(self):
        ret = []
        for rec in self.current_tr_records:
            if rec is not None:
                ret.append(rec.vcfrecord)
            else:
                ret.append(None)
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
                if prev_records[idx] is None:
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
        self.done = all([item is None or item.vcfrecord is None for item in self.current_tr_records])
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
