#!/usr/bin/env python3

"""
Merge Module includes code and classes for merging record clusters.
Work in progress

"""
from callsetmerger.recordcluster import ClusterGraph, RecordResolver
from vcf.model import _Record,_Substitution,_Call, make_calldata_tuple
import traceback


# class MergeOverlappingRegion:
#     def __init__(self, overlapping_region):
#         self.ov_region = overlapping_region
#         self.rc_merge_objects = []
#         for rc in overlapping_region.RecordClusters:
#             self.rc_merge_objects.append(RecordClusterMerger(rc))

#     def GetVCFLines(self):
#         ret = []
#         for rc_merge_obj in self.rc_merge_objects:
#             line = rc_merge_obj.GetVCFLine()
#             if line is not None:
#                 ret.append(line)
#         return ret



def GetGTForAlleles(resolved_call):
    # resolved_call is a list of two alleles
    if resolved_call[0] == None and resolved_call[1] == None:
        return '.'
    if resolved_call[0] == None or resolved_call[1] == None:
        raise ValueError('Only one genotype set to None!')
    GTs = []
    for al in resolved_call:
        GTs.append(str(al.genotype_idx))
    return '/'.join(GTs)

class RecordClusterMerger:
    def __init__(self, rc, samples):
        self.record_cluster = rc
        self.record_resolver = RecordResolver(rc)
        self.samples = samples
        self.record_template = rc.record_objs[0].cyvcf2_record
        
    def ResolveAllSampleCalls(self):
        res_calls = {}
        for sample in self.samples:
            samp_call = self.record_cluster.GetSampleCall(sample)   # Get calls for this sample
            resolved_connected_comps = self.record_resolver.GetConnectedCompForSingleCall(samp_call)      # Resolve call for this sample (pick which caller(s) call we use)
            res_calls[sample] = self.record_resolver.ResolveSequenceForSingleCall(resolved_connected_comps, samp_call)
        return res_calls

    

    def GetCyVCFRecord(self):
        # TODO update template with info for resolved genotypes
        return self.record_template

    def GetPyVCFRecord(self):
        ALTS = []
        for alt_allele in self.record_template.ALT:
            ALTS.append(_Substitution(alt_allele))
        INFO = {}
        for info in self.record_template.INFO:
            INFO[info[0]] = self.record_template.INFO.get(info[0])
        # FORMAT = self.record_template.FORMAT
        FORMAT = ['GT','SRC']
        samp_fmt = make_calldata_tuple(FORMAT)

        res_calls = self.ResolveAllSampleCalls()
        print(res_calls)
        SAMPLES=[]
        for sample in self.record_cluster.samples:
            SAMPLES.append(_Call(self, sample, samp_fmt(GT=GetGTForAlleles(res_calls[sample]),SRC='hipstr')))
        return _Record(self.record_template.CHROM, 
            self.record_template.POS, 
            self.record_template.ID,
            self.record_template.REF,
            ALTS,
            self.record_template.QUAL,
            self.record_template.FILTER,
            INFO,
            ':'.join(FORMAT),
            dict([(x,i) for (i,x) in enumerate(self.record_cluster.samples)]),
            SAMPLES)


# For each sample, we have 1 object that we pass the graph to -> call for that sample
