#!/usr/bin/env python3

"""
RecordCluster object includes mergeable records. Mergeable records are defined as records that overlap and share a motif
Work in progress

"""
from typing import List

import trtools.utils.tr_harmonizer as trh
from trtools.utils.utils import GetCanonicalMotif

from enum import Enum
import networkx as nx


class AlleleType(Enum):
    Reference = 0
    Alternate = 1


ColorDict = {trh.VcfTypes.advntr: "tab:blue",
             trh.VcfTypes.eh: "tab:orange",
             trh.VcfTypes.hipstr: "tab:red",
             trh.VcfTypes.gangstr: "tab:green",
             # trh.VcfTypes.popstr: "tab:yellow",
             }
convert_type_to_idx = {trh.VcfTypes.advntr: 0,
                       trh.VcfTypes.eh: 1,
                       trh.VcfTypes.hipstr: 2,
                       trh.VcfTypes.gangstr: 3,
                       # trh.VcfTypes.popstr: 4,
                       }


def GetVcfTypesKey():
    return convert_type_to_idx.keys()


class Allele:
    def __init__(self, vcf_type, atype, allele_sequence, diff_from_ref, genotype_idx):
        if atype not in [AlleleType.Reference, AlleleType.Alternate]:
            raise ValueError('Unknown allele type: ' + atype)
        self.allele_type = atype
        self.allele_sequence = allele_sequence
        self.allele_size = diff_from_ref
        self.genotype_idx = genotype_idx
        self.vcf_type = vcf_type

    def GetLabel(self):
        return self.vcf_type.name + str(self.allele_size)


class RecordObj:
    def __init__(self, vcf_type, rec):
        self.record = rec
        self.vcf_type = vcf_type
        self.hm_record = trh.HarmonizeRecord(vcf_type, rec)
        self.ref = self.hm_record.ref_allele
        self.motif = self.hm_record.motif
        self.canonical_motif = GetCanonicalMotif(self.motif)


class RecordCluster:
    def __init__(self, recobjs):
        self.motif = recobjs[0].motif
        self.canonical_motif = GetCanonicalMotif(self.motif)
        self.vcf_types = [False] * len(convert_type_to_idx.keys())
        for rec in recobjs:
            self.vcf_types[convert_type_to_idx[rec.vcf_type]] = True
            if rec.canonical_motif != self.canonical_motif:
                raise ValueError('RecordObjs in a cluster must share canonical motif.\n' +
                                 'Cannot add RO with ' + rec.canonical_motif + ' canonical motif to cluster' +
                                 'with ' + self.canonical_motif + 'canonical motif.')
        self.record_objs = recobjs

    def AppendRecordObject(self, ro):
        if self.canonical_motif != ro.canonical_motif:
            raise ValueError('Canonical motif of appended record object mush match record cluster.')
        self.record_objs.append(ro)
        self.vcf_types[convert_type_to_idx[ro.vcf_type]] = True

    def GetVcfTypesTuple(self):
        return tuple(self.vcf_types)

    def GetAlleleList(self):
        alist = []
        for ro in self.record_objs:
            ref = ro.hm_record.ref_allele
            alist.append(Allele(ro.vcf_type, AlleleType.Reference, ref, 0, 0))
            altnum = 1
            for alt in ro.hm_record.alt_alleles:
                alist.append(Allele(ro.vcf_type, AlleleType.Alternate, alt, len(alt) - len(ref), altnum))
                altnum += 1
        return alist


class ClusterGraph:
    def __init__(self, allele_list):
        self.allele_list = allele_list
        self.graph = nx.Graph()
        self.labels = {}
        for al in allele_list:
            self.graph.add_node(al)
            self.labels[al] = al.GetLabel()
        self.colors = []
        for node in self.graph.nodes():
            self.colors.append(ColorDict[node.vcf_type])
        for nd1 in self.graph.nodes():
            for nd2 in self.graph.nodes():
                if nd1 == nd2:
                    continue
                if nd1.allele_size == nd2.allele_size \
                        and not self.graph.has_edge(nd1, nd2) \
                        and nd1.vcf_type != nd2.vcf_type:
                    self.graph.add_edge(nd1, nd2)
