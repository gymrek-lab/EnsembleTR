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


class OutVCFRecord:
    def __init__(self, res_pas, rc):
        self.pre_alleles = res_pas
        self.record_cluster = rc
        self.ref = ""
        self.alts = []
        self.sample_to_GT = {}
        self.sample_to_NCOPY = {}
        self.sample_to_SRC = {}
        self.sample_to_INPUTS = rc.GetAllInputs()
        # Add INFO fields as well

        for sample in self.pre_alleles:
            GT_list = []
            NCOPY_list  = []
            SRC_list = []
            for pa in self.pre_alleles[sample]:
                if self.ref == "":
                    self.ref = pa.reference_sequence
                # TODO Fix assertion
                #assert self.ref == pa.ref_seq

                if pa.allele_sequence != self.ref:
                    # we have alt allele
                    if pa.allele_sequence not in self.alts:
                        self.alts.append(pa.allele_sequence)
                    GT_list.append(str(self.alts.index(pa.allele_sequence) + 1))
                    NCOPY_list.append(str(pa.allele_ncopy))
                else:
                    # ref allele
                    GT_list.append('0')
                    NCOPY_list.append(str(pa.reference_ncopy))
                for caller in pa.support:
                    if caller.name not in SRC_list:
                        SRC_list.append(caller.name)
            if len(GT_list) == 0:
                GT_list = ['.']
            if len(SRC_list) == 0:
                SRC_list = ['.']
            if len(NCOPY_list) == 0:
                NCOPY_list = ['.']
            self.sample_to_GT[sample] = '/'.join(GT_list)
            self.sample_to_NCOPY[sample] = ','.join(NCOPY_list)
            self.sample_to_SRC[sample] = ','.join(SRC_list)


def get_info_string(data):
    out_recs = []
    for key in data:
        out_recs.append(str(key) + '=' + str(data[key]))
    return ';'.join(out_recs)

class RecordClusterOutput:
    def __init__(self, rc, samples):
        self.record_cluster = rc
        self.record_resolver = RecordResolver(rc)
        self.samples = samples
        self.record_template = rc.record_objs[0].cyvcf2_record
        
    def ResolveAllSampleCalls(self):
        res_pas = {}
        res_cer = {}
        for sample in self.samples:
            samp_call = self.record_cluster.GetSampleCall(sample)   # Get calls for this sample
            resolved_connected_comp_ids, certain_cc = self.record_resolver.GetConnectedCompForSingleCall(samp_call)      # Resolve call for this sample (pick which caller(s) call we use)
            res_cer[sample] = certain_cc
            res_pas[sample] = self.record_resolver.ResolveSequenceForSingleCall(resolved_connected_comp_ids, samp_call)
        return res_pas, res_cer

    

    def GetCyVCFRecord(self):
        # TODO update template with info for resolved genotypes
        return self.record_template

    # def GetPyVCFRecord(self):
    #     INFO = {}
    #     for info in self.record_template.INFO:
    #         INFO[info[0]] = self.record_template.INFO.get(info[0])
    #     # FORMAT = self.record_template.FORMAT
    #     FORMAT = ['GT','SRC']
    #     samp_fmt = make_calldata_tuple(FORMAT)

    #     # Create list of pre-alleles
    #     res_pas = self.ResolveAllSampleCalls()
    #     print(res_pas)
    #     out_rec = OutVCFRecord(res_pas)
    #     SAMPLES=[]
    #     for sample in self.record_cluster.samples:
    #         SAMPLES.append(_Call(self, sample, samp_fmt(GT=out_rec.sample_to_GT[sample],SRC=out_rec.sample_to_SRC[sample])))
        
    #     ALTS = []
    #     for alt_allele in out_rec.alts:
    #         ALTS.append(_Substitution(alt_allele))
    #     return _Record(self.record_template.CHROM, 
    #         self.record_template.POS, 
    #         self.record_template.ID,
    #         out_rec.ref,
    #         ALTS,
    #         self.record_template.QUAL,
    #         self.record_template.FILTER,
    #         INFO,
    #         ':'.join(FORMAT),
    #         dict([(x,i) for (i,x) in enumerate(self.record_cluster.samples)]),
    #         SAMPLES)

    def GetRawVCFRecord(self):
        FORMAT = ['GT', 'NCOPY', 'SRC','CERT','INPUTS']
        # Create list of pre-alleles
        res_pre_alleles, res_cert = self.ResolveAllSampleCalls()

        INFO_DICT = {'START': self.record_cluster.first_pos,
            'END': self.record_cluster.last_end,
            'PERIOD': len(self.record_cluster.motif),
            'RU': self.record_cluster.motif}
        # print(res_pas)
        out_rec = OutVCFRecord(res_pre_alleles, self.record_cluster)
        SAMPLE_DATA=[]
        for sample in self.record_cluster.samples:
            SAMPLE_DATA.append(':'.join(
                [out_rec.sample_to_GT[sample],
                out_rec.sample_to_NCOPY[sample],
                out_rec.sample_to_SRC[sample],
                str(res_cert[sample]),
                out_rec.sample_to_INPUTS[sample]
                ]
                ))
        
        ALTS = []
        INFO = get_info_string(INFO_DICT)
        for alt_allele in out_rec.alts:
            ALTS.append(_Substitution(alt_allele))
        # TODO remove record template and get information from pre allele
        return '\t'.join([str(self.record_template.CHROM), 
            str(self.record_template.POS), 
            '.',
            out_rec.ref,
            ','.join(out_rec.alts),
            '.',
            '.',
            INFO,
            ':'.join(FORMAT),
            '\t'.join(SAMPLE_DATA)]) + '\n'

# For each sample, we have 1 object that we pass the graph to -> call for that sample
