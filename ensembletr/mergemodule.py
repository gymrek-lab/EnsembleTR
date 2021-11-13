#!/usr/bin/env python3

"""
Merge Module includes code and classes for merging record clusters.
Work in progress

"""
from . import recordcluster as recordcluster

class RecordClusterOutput:
    def __init__(self, rc, samples):
        self.record_cluster = rc
        self.record_resolver = recordcluster.RecordResolver(rc)
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

    def GetRawVCFRecord(self):
        FORMAT = ['GT', 'NCOPY', 'SRC','CERT','INPUTS']
        # Create list of pre-alleles
        res_pre_alleles, res_cert = self.ResolveAllSampleCalls()

        INFO_DICT = {'START': self.record_cluster.first_pos,
                     'END': self.record_cluster.last_end,
                     'PERIOD': len(self.record_cluster.canonical_motif),
                     'RU': self.record_cluster.canonical_motif,
                     'METHODS': "|".join([str(int(item)) for item in self.record_cluster.vcf_types])}

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
        
        ALTS = out_rec.alts
        if len(ALTS) == 0: ALTS.append(".")
        INFO = get_info_string(INFO_DICT)
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

