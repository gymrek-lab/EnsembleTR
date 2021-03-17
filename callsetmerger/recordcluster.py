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

class PreAllele:
    def __init__(self, ref_seq, seq, callers):
        self.ref_seq = ref_seq
        self.seq = seq
        self.support = callers
    
    def add_support(self, callers):
        for caller in callers:
            if caller not in self.support:
                self.support.append(caller)

    def __str__(self):
        return "Ref: " + self.ref_seq + "\nSeq: " + self.seq + "\nSupp: " + str(self.support)
    def __repr__(self):
        return "Ref: " + self.ref_seq + "\nSeq: " + self.seq + "\nSupp: " + str(self.support)
    
class Allele:
    def __init__(self, vcf_type, atype, allele_sequence, reference_sequence, diff_from_ref, genotype_idx):
        if atype not in [AlleleType.Reference, AlleleType.Alternate]:
            raise ValueError('Unknown allele type: ' + atype)
        self.allele_type = atype
        self.allele_sequence = allele_sequence
        self.reference_sequence = reference_sequence
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
            alist.append(Allele(ro.vcf_type, AlleleType.Reference, ref, ref, 0, 0))
            altnum = 1
            for alt in ro.hm_record.alt_alleles:
                alist.append(Allele(ro.vcf_type, AlleleType.Alternate, alt, ref, len(alt) - len(ref), altnum))
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
        self.subgraphs = [self.graph.subgraph(c).copy() for c in self.sorted_connected_comps]
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

    def GetSubgraphForNode(self, node):
        if node is None:
            return None
        for subg in self.subgraphs:
            if subg in self.GetConnectedComponentSubgraphs():
                if node in subg.nodes():
                    return subg
        return None

    def GetSortedConnectedComponents(self):
        return self.sorted_connected_comps

    def GetConnectedComponentSubgraphs(self):
        return self.subgraphs

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


def add_ccsg_support(sg_support, sg):
    if sg is not None:
        if sg not in sg_support.keys():
            sg_support[sg] = 1
        else:
            sg_support[sg] += 1
    return sg_support

def get_seq_to_pa_idx_dict(pa_list):
    seq_to_pa = {}
    for idx, pa in enumerate(pa_list):
        seq_to_pa[pa.seq] = idx
    return seq_to_pa
        

def add_preallele_support(pa_list, pa):
    if pa is not None:
        seq_to_pa_idx = get_seq_to_pa_idx_dict(pa_list)
        if pa.seq not in seq_to_pa_idx.keys():
            pa_list.append(pa)
        else:
            pa_list[seq_to_pa_idx[pa.seq]].add_support(pa.support)
    return pa_list

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
        list of connected component subgraphs corresponding to resolved call
        """
        # Assign connected component subgraphs to alleles
        conn_comp_sg_dict = {}
        ccsg_support = {}
        num_valid_methods = 0
        for method in samp_call:
            # check for no calls
            if samp_call[method][0] == -1:
                ccsg0 = None
                ccsg1 = None
            else:
                ccsg0 = self.rc_graph.GetSubgraphForNode(self.rc_graph.GetNodeObject(method, samp_call[method][0]))
                ccsg_support = add_ccsg_support(ccsg_support, ccsg0)
                ccsg1 = self.rc_graph.GetSubgraphForNode(self.rc_graph.GetNodeObject(method, samp_call[method][1]))
                ccsg_support = add_ccsg_support(ccsg_support, ccsg1)
                num_valid_methods += 1
            conn_comp_sg_dict[method] = [ccsg0, ccsg1]

        sorted_ccsg_support = dict(sorted(ccsg_support.items(), key=lambda item: item[1], reverse=True))
        ret_sgs = []
        certain = True
        for sg in sorted_ccsg_support:
            # If all alleles in all methods point to one cc
            # for Hom require support to be perfect (2x num_valide) if not, report not certain
            # for Het do the same for support two ccs each with num_valid
            if sorted_ccsg_support[sg] == num_valid_methods * 2: 
                return [sg, sg], certain
            elif sorted_ccsg_support[sg] == num_valid_methods:
                if len(ret_sgs) < 2:
                    ret_sgs.append(sg)
                else:
                    certain = False
                    # We have already added 2 top supported ccs, but we still have another cc
                    print("WARNING: extra cc", sg)
                    
            else:
                certain = False
                print("WARNING: at least one CC is not fully supported by either one or both alleles", sg)
        return ret_sgs, certain
        


    def ResolveSequenceForSingleCall(self, ccsg_list, samp_call):
        # Next update: Melissa's idea
        # Pre resolve possible sequences for each CC. Only check if samp_call contains the caller and genot_idx of the resolved seq
        # Add support based on overlap between samp_call and methods in cc
        pre_allele_list = []
        for ccsg in ccsg_list:
            # If number of nodes == number of callers: 1-1-1
            uniq_callers = []
            caller_to_nodes = {}
            for node in ccsg.nodes():
                if node.vcf_type not in caller_to_nodes:
                    caller_to_nodes[node.vcf_type] = [node]
                else:
                    caller_to_nodes[node.vcf_type].append(node)
                
                if node.vcf_type not in uniq_callers:
                    uniq_callers.append(node.vcf_type)
            if len(uniq_callers) == len(ccsg.nodes):
                tmp_node = list(ccsg.nodes())[0]
                pre_allele_list.append(PreAllele(tmp_node.reference_sequence, tmp_node.allele_sequence, uniq_callers))
            else:
                # if hipstr exists:
                if trh.VcfTypes.hipstr in uniq_callers:
                    # only one hipstr node
                    if len(caller_to_nodes[trh.VcfTypes.hipstr]) == 1:
                        tmp_node = caller_to_nodes[trh.VcfTypes.hipstr][0]
                        pa = PreAllele(tmp_node.reference_sequence, tmp_node.allele_sequence, [trh.VcfTypes.hipstr])
                        for node in ccsg:
                            if node != tmp_node and node.allele_sequence == tmp_node.allele_sequence:
                                pa.add_support([node.vcf_type])
                        
                    else:
                        # More than one hipstr node: decide which one matches the samp_call
                        tmp_node = None
                        for allele in ccsg:
                            if allele.vcf_type == trh.VcfTypes.hipstr and allele.genotype_idx in samp_call[trh.VcfTypes.hipstr]:
                                tmp_node = allele
                        if tmp_node is None:
                            raise ValueError("Could not find HipSTR allele matching samp_call in this ccsg!")
                        pa = PreAllele(tmp_node.reference_sequence, tmp_node.allele_sequence, [trh.VcfTypes.hipstr])
                        for node in ccsg:
                            if node != tmp_node and node.allele_sequence == tmp_node.allele_sequence:
                                pa.add_support([node.vcf_type])
                else:
                    raise ValueError("HipSTR doesn't exist and we have a discrepancy!")

        if len(pre_allele_list) >= 2:
            return pre_allele_list[0:2]
        elif len(pre_allele_list) == 1:
            return [pre_allele_list[0], pre_allele_list[0]]
        else:
            # raise ValueError("Too few pre_alleles")
            print('Warning: too few pre alleles: ', pre_allele_list)
            return pre_allele_list
             
        # node_dict = {}
        # for method in samp_call:
        #     if samp_call[method][0] == -1:
        #         nd0 = None
        #         nd1 = None
        #     else:
        #         nd0 = self.rc_graph.GetNodeObject(method, samp_call[method][0])
        #         nd1 = self.rc_graph.GetNodeObject(method, samp_call[method][1])
        #         if nd0 in ccsg_list[0].nodes() or nd0 in ccsg_list[1].nodes:
        #             add_preallele_support(pre_allele_list, PreAllele(nd0.reference_sequence, nd0.allele_sequence, [nd0.vcf_type]))
        #         if nd1 in ccsg_list[0].nodes() or nd1 in ccsg_list[1].nodes:
        #             add_preallele_support(pre_allele_list, PreAllele(nd1.reference_sequence, nd1.allele_sequence, [nd1.vcf_type]))
        #     node_dict[method] = [nd0, nd1]
        
        # # First iteration, just return the first two pre alleles
        # if len(pre_allele_list) >= 2:
        #     return pre_allele_list[0:2]
        # elif len(pre_allele_list) == 1:
        #     return [pre_allele_list[0], pre_allele_list[0]]
        # else:
        #     # raise ValueError("Too few pre_alleles")
        #     print('Warning: too few pre alleles: ', pre_allele_list)
        #     return pre_allele_list

        # for cc in self.rc_graph.
        # Check if 1-to-1-to-1
        # TODO
        #

