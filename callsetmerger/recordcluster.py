#!/usr/bin/env python3

"""
RecordCluster object includes mergeable records. Mergeable records are defined as records that overlap and share a motif
Work in progress

"""
from typing import List

import trtools.utils.tr_harmonizer as trh
from trtools.utils.utils import GetCanonicalMotif
from collections import Counter
from enum import Enum
import networkx as nx
import numpy as np


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
        if self.allele_type == AlleleType.Reference:
            return self.vcf_type.name + '_*'
        else:
            return self.vcf_type.name + '_' + str(self.allele_size)


class RecordObj:
    def __init__(self, vcf_type, VCF, rec):
        self.cyvcf2_record = rec
        self.VCF = VCF
        self.samples = VCF.samples
        self.vcf_type = vcf_type
        self.hm_record = trh.HarmonizeRecord(vcf_type, rec)
        self.ref = self.hm_record.ref_allele
        self.motif = self.hm_record.motif
        self.canonical_motif = GetCanonicalMotif(self.motif)

    def GetSamples(self):
        return self.cyvcf2_record.samples

    def GetVcfRegion(self):
        return str(self.cyvcf2_record.CHROM) + ':' + str(self.cyvcf2_record.POS)
    
    def GetROSampleCall(self, sample):
        samp_idx = self.samples.index(sample)
        return self.cyvcf2_record.genotypes[samp_idx]


class RecordCluster:
    def __init__(self, recobjs):
        self.motif = recobjs[0].motif
        self.canonical_motif = GetCanonicalMotif(self.motif)
        self.vcf_types = [False] * len(convert_type_to_idx.keys())
        self.samples = recobjs[0].samples
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

    def GetVcfRegions(self):
        ret_list = []
        for ro in self.record_objs:
            ret_list.append(ro.GetVcfRegion())

        return ret_list

    def GetSampleCall(self, sample):
        ret_dict = {}
        for rec in self.record_objs:
            if rec.vcf_type in ret_dict:
                raise ValueError("Multiple records with same VCF type: " + str(rec.vcf_type))
            ret_dict[rec.vcf_type] = rec.GetROSampleCall(sample)
        return ret_dict

# OverlappingRegion includes 1 or more record clusters.
# These RCs can have different motifs.
class OverlappingRegion:
    def __init__(self, rcs):
        self.RecordClusters = rcs

    def GetCanonicalMotifs(self):
        ret = []
        for rc in self.RecordClusters:
            ret.append(rc.canonical_motif)
        return ret


class ClusterGraph:
    def __init__(self, record_cluster):
        allele_list = record_cluster.GetAlleleList()
        self.allele_list = allele_list
        self.graph = nx.Graph()
        self.labels = {}
        self.vcf_types = [False] * len(convert_type_to_idx.keys())
        for al in allele_list:
            self.graph.add_node(al)
            self.labels[al] = al.GetLabel()
            self.vcf_types[convert_type_to_idx[al.vcf_type]] = True
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
        self.sorted_connected_comps = sorted(nx.algorithms.components.connected_components(self.graph), key=len, reverse=True)
        # Identify representative nodes:
        # These are the nodes that are used as representatives for their connected comp
        # for component in self.GetSortedConnectedComponents():
        #     num_unique_callers_in_component = 0
        #     # If we see
        #     callers_seen = []
        #     for allele in component:
        #         if allele.vcf_type not in callers_seen:
        #             callers_seen.append(allele.vcf_type)
        #             num_unique_callers_in_component += 1
        #     list_unique_caller_nodes_in_conn_comp.append(num_unique_callers_in_component)

    def GetNodeObject(self, vcf_type, genotype_idx):
        for allele in self.graph.nodes:
            if allele.vcf_type == vcf_type and allele.genotype_idx == genotype_idx:
                return allele
        return None

    def GetConnectedCompForNode(self, node):
        if node is None:
            return None
        for component in self.GetSortedConnectedComponents():
            if node in component:
                return component
        return None

    def GetSortedConnectedComponents(self):
        return self.sorted_connected_comps

    def GetSingularityScore(self):
        # 1 means all components at least 1-to-1 (could be 2-to-1)
        # Lower than 1 means there are singular nodes
        num_callers_in_graph = sum(self.vcf_types)
        list_unique_caller_nodes_in_conn_comp = []
        for component in self.GetSortedConnectedComponents():
            num_unique_callers_in_component = 0
            callers_seen = []
            for allele in component:
                if allele.vcf_type not in callers_seen:
                    callers_seen.append(allele.vcf_type)
                    num_unique_callers_in_component += 1
            list_unique_caller_nodes_in_conn_comp.append(num_unique_callers_in_component)
        return np.mean(list_unique_caller_nodes_in_conn_comp) / float(num_callers_in_graph)

    def GetConfusionScore(self):
        # Measure of 1-to-1 ness
        # average over all connected components:
        # For each connected component: number of nodes / number of unique caller nodes
        list_comp_confusion_score = []
        for component in self.GetSortedConnectedComponents():
            num_nodes = 0
            num_unique_callers_in_component = 0
            callers_seen = []
            for allele in component:
                num_nodes += 1
                if allele.vcf_type not in callers_seen:
                    callers_seen.append(allele.vcf_type)
                    num_unique_callers_in_component += 1
            list_comp_confusion_score.append(float(num_nodes) / float(num_unique_callers_in_component))
        return np.mean(list_comp_confusion_score)

    # Could be independent (take graph as input)
    # TODO a function to assign ref and alt alleles for a given sample
    # Boundaries? if there is discrepancy
    # -> Pos filed to be consistent


def add_cc_support(cc_support, cc):
    if cc is not None:
        if cc not in cc_support.keys():
            cc_support[cc] = 1
        else:
            cc_support[cc] += 1
    return cc_support


class RecordResolver:
    def __init__(self, rc):
        self.record_cluster = rc
        self.rc_graph = ClusterGraph(rc)

    def GetConnectedCompForSingleCall(self, samp_call):
        r"""

        Parameters
        ----------
        samp_call: dict(vcftype:allele)
                allele: [al1, al2, BOOL] or [-1, BOOL] for no calls

        Returns
        -------
        list of connected components corresponding to resolved call
        """
        # Assign connected components to alleles
        conn_comp_dict = {}
        cc_support = {}
        num_valid_methods = 0
        for method in samp_call:
            # check for no calls
            if samp_call[method][0] == -1:
                cc0 = None
                cc1 = None
            else:
                cc0 = self.rc_graph.GetConnectedCompForNode(self.rc_graph.GetNodeObject(method, samp_call[method][0]))
                cc_support = add_cc_support(cc_support, cc0)
                cc1 = self.rc_graph.GetConnectedCompForNode(self.rc_graph.GetNodeObject(method, samp_call[method][1]))
                cc_support = add_cc_support(cc_support, cc1)
                num_valid_methods += 1
            conn_comp_dict[method] = [cc0, cc1]

        sorted_cc_support = dict(sorted(cc_support.items(), key=lambda item: item[1]))
        ret_ccs = []
        for cc in sorted_cc_support:
            # If all alleles in all methods point to one cc
            if sorted_cc_support[cc] <= num_valid_methods * 2 and sorted_cc_support[cc] > (num_valid_methods - 1) * 2:
                return [cc, cc]
            if len(ret_ccs) < 2:
                ret_ccs.append(cc)
            else:
                # We have already added 2 top supported ccs, but we still have another cc
                print("WARNING: extra cc", cc)

        return ret_ccs


    def ResolveSequenceForSingleCall(self, conn_comp_list, samp_call):
        node_dict = {}
        for method in samp_call:
            if samp_call[method][0] == -1:
                nd0 = None
                nd1 = None
            else:
                nd0 = self.rc_graph.GetNodeObject(method, samp_call[method][0])
                nd1 = self.rc_graph.GetNodeObject(method, samp_call[method][1])
            node_dict[method] = [nd0, nd1]
        
        # First iteration, just return the first allele in cc that has samp_call

        # for cc in self.rc_graph.
        # Check if 1-to-1-to-1
        # TODO
        #

