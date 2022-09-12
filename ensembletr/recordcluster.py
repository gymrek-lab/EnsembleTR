"""
Classes to keep track of mergeable records
"""

import trtools.utils.tr_harmonizer as trh
from trtools.utils.utils import GetCanonicalMotif
from collections import defaultdict
from enum import Enum
import networkx as nx
import numpy as np
import math

from . import utils as utils

CC_PREFIX = 'cc'

convert_type_to_idx = {trh.VcfTypes.advntr: 0,
                       trh.VcfTypes.eh: 1,
                       trh.VcfTypes.hipstr: 2,
                       trh.VcfTypes.gangstr: 3,
                       }

class RecordObj:
    r"""
    Main object to store a VCF record and associated metadata

    Parameters
    ----------
    rec : cyvcf2.VCF.vcfrecord
       VCF record for the object
    vcf_type: trh.TRRecordHarmonizer.vcftype
       Type of the VCF file
    """
    def __init__(self, rec, vcf_type, vcf_samples):
        self.cyvcf2_record = rec
        self.vcf_type = vcf_type
        self.hm_record = trh.HarmonizeRecord(vcf_type, rec)
        self.canonical_motif = GetCanonicalMotif(self.hm_record.motif)
        self.prepend_seq = ''
        self.append_seq = ''
        self.vcf_samples = vcf_samples

    def GetCalledAlleles(self):
        r"""
        Get the set of all alternate alleles
        with at least one calls

        Returns
        -------
        al_idx : set of int
            Set of called alleles (based on REF/ALT fields)
        """
        al_idx = set()
        for call in self.hm_record.vcfrecord.genotypes:
            if call[0] != -1 and len(call) == 3:
                al_idx.add(call[0])
                al_idx.add(call[1])
        return al_idx

    def GetROSampleCall(self, sample):
        r"""
        Get the genotype of a single sample

        Parameters
        ----------
        samp_idx : int
           Index of the sample

        Returns
        -------
        sample_gt : [int, int, bool]
           Corresponds to the genotype alleles and phasing
           (based on cyvcf2 representation)
        """
        samp_idx = self.vcf_samples.index(sample)
        return self.cyvcf2_record.genotypes[samp_idx]

    def GetSampleString(self, sample):
        r"""
        Get a user-readable string of a sample's genotype

        Parameters
        ----------
        samp_idx : int
           Index of the sample

        Returns
        -------
        callstr : str
            Format is "caller=allele1,allele2"
            Caller is one of gangstr/hipstr/eh/advntr
            Alleles are given in copy number        
        """
        samp_call = self.GetROSampleCall(sample)
        if samp_call is None or samp_call[0] == -1:
            sampdata = "."
        else:
            call_ncopy_list = []
            for idx in samp_call[0:2]:
                if idx == 0:
                    call_ncopy_list.append(str(round(self.hm_record.ref_allele_length,2)))
                else:
                    call_ncopy_list.append(str(round(self.hm_record.alt_allele_lengths[idx - 1],2)))
            sampdata = ",".join(call_ncopy_list)
        callstr = "%s=%s"%(self.vcf_type.name, sampdata)
        return callstr

    def GetScore(self, sample):
        r"""
        Get the score of a sample's genotype
        For HipSTR/GangSTR, use "FORMAT/Q"
        For adVNTR: use "FORMAT/ML"
        For ExpansionHunter, use a custom score based on
            FORMAT/REPCN and FORMAT/REPCI

        Parameters
        ----------
        samp_idx : int
           Index of the sample

        Returns
        -------
        score : float
           Indicates confidence in the call (0=low, 1=high)      
        """
        samp_idx = self.vcf_samples.index(sample)
        if self.vcf_type == trh.VcfTypes.advntr:
            return self.cyvcf2_record.format('ML')[samp_idx][0]
        elif self.vcf_type in [trh.VcfTypes.hipstr, trh.VcfTypes.gangstr]:
            return self.cyvcf2_record.format('Q')[samp_idx][0]
        elif self.vcf_type == trh.VcfTypes.eh:
            REPCI = self.cyvcf2_record.format('REPCI')[samp_idx]
            REPCN = self.cyvcf2_record.format('REPCN')[samp_idx]
            if REPCI == "." or REPCN == ".":
                return 0
            length = len(self.cyvcf2_record.INFO['RU'])
            return utils.GetEHScore(REPCI, REPCN, length)
        else:
            return 0 # shouldn't happen

class RecordCluster:
    r"""
    Class to keep track of a list of mergeable records

    Parameters
    ----------
    recobjs : list of RecordObj
       list of record objects to be merged
    ref_genome : pyfaidx.Fasta
       reference genome
    canon_motif : str
       canonical repeat motif
    samples : list of str
       List of samples to analyze

    """
    def __init__(self, recobjs, ref_genome, canon_motif, samples):
        self.canonical_motif = canon_motif
        self.vcf_types = [False] * len(convert_type_to_idx.keys())
        self.samples = samples
        self.fasta = ref_genome        
        self.record_objs = recobjs
        self.first_pos = -1
        self.last_pos = -1
        self.chrom = recobjs[0].cyvcf2_record.CHROM
        self.update()

    def AppendRecordObject(self, ro):
        r"""
        Add a record object to the RecordCluster
        and update associated metadata

        Parameters
        ----------
        ro : RecordObj
           the Record Object to be added
        """
        self.record_objs.append(ro)
        self.update()

    def update(self):
        r"""
        Update prepend/append sequences
        of individual record objects so they
        start and end at the same location.
        Extends all alleles to the maximum region
        spanned by all records in the cluster.
        """
        self.first_pos = min([rec.cyvcf2_record.POS for rec in self.record_objs])
        self.last_end = max([rec.cyvcf2_record.end for rec in self.record_objs])

        for rec in self.record_objs:
            self.vcf_types[convert_type_to_idx[rec.vcf_type]] = True
            chrom = rec.cyvcf2_record.CHROM
            if rec.cyvcf2_record.POS > self.first_pos:
                # Found a record that starts after
                # Should prepend the record
                rec.prepend_seq = self.fasta[chrom][self.first_pos : rec.cyvcf2_record.POS].seq.upper()
            
            if rec.cyvcf2_record.end < self.last_end:
                # Found a record that ends before last end
                # Should append the record
                rec.append_seq = self.fasta[chrom][rec.cyvcf2_record.end : self.last_end].seq.upper()

    def GetRawCalls(self):
        r"""
        Get string of inputs to use for debug info in VCF

        Returns
        -------
        out_dict : (dict of str: str)
           Key=sample, Value=comma-separated list of genotypes
        """
        out_dict = {}
        for sample in self.samples:
            samp_call_list = [ro.GetSampleString(sample) for ro in self.record_objs]
            out_dict[sample] = '|'.join(samp_call_list)
        return out_dict        

    def GetSampleCall(self, sample):
        r"""
        Get calls for an individual sample
        
        Parameters
        ----------
        sample : str

        Returns
        -------
        ret_dict : (dict of str: str)
           Key=VCF type, Value=sample call
        """
        ret_dict = {}
        for rec in self.record_objs:
            if rec.vcf_type in ret_dict:
                raise ValueError("Multiple records with same VCF type: " + str(rec.vcf_type))
            ret_dict[rec.vcf_type] = rec.GetROSampleCall(sample)
        return ret_dict

    def GetQualScore(self, sample):
        r"""
        Get quality scores for an individual sample

        Paramters
        ---------
        sample : str

        Returns
        -------
        ret_dict : (dict of str: str)
           Key=VCF type, Value=sample quality score
        """
        ret_dict = {}
        for rec in self.record_objs:
            if rec.vcf_type in ret_dict:
                raise ValueError("Multiple records with same VCF type: " + str(rec.vcf_type))
            score = rec.GetScore(sample)
            if ~np.isnan(score):
                ret_dict[rec.vcf_type] = score
        return ret_dict

class OverlappingRegion:
    """
    OverlappingRegion includes 1 or more record clusters.
    These RCs can have different motifs.
    """
    def __init__(self, rcs):
        self.RecordClusters = rcs

    def GetCanonicalMotifs(self):
        ret = []
        for rc in self.RecordClusters:
            ret.append(rc.canonical_motif)
        return ret

class Allele:
    """
    Object to store alleles (nodes)

    Parameters
    ----------
    ro : RecordObj
       record object from which the allele originates
    al_idx : int
       Allele index (0=ref, 1+ =alt alleles)

    Attributes
    ----------
    allele_sequence : str
       Full string of the allele
    allele_size : int
       Length difference from the reference genome (in bp)

    TODO: get rid of allele_ncopy, reference_ncopy
    """
    def __init__(self, ro, al_idx):
        self.record_object = ro
        self.al_idx = al_idx
        alleles = []
        if ro.hm_record.full_alleles is None:
            alleles = [ro.hm_record.ref_allele]+ro.hm_record.alt_alleles
        else:
            alleles = [ro.hm_record.full_alleles[0]] + ro.hm_record.full_alleles[1]
        allele_lengths = [ro.hm_record.ref_allele_length] + ro.hm_record.alt_allele_lengths
        self.reference_sequence = ro.prepend_seq + alleles[0] + ro.append_seq
        self.allele_sequence = ro.prepend_seq + alleles[al_idx] + ro.append_seq
        self.allele_size = len(self.allele_sequence) - len(self.reference_sequence)
        self.allele_ncopy = round(allele_lengths[al_idx],2)
        self.reference_ncopy = round(ro.hm_record.ref_allele_length,2)
        self.exp_flag = (len(ro.prepend_seq) > 0) or (len(ro.append_seq) > 0)


    def GetVCFType(self):
        return self.record_object.vcf_type

class PreAllele:
    def __init__(self, allele, callers):
        self.reference_sequence = allele.reference_sequence
        self.allele_sequence = allele.allele_sequence
        self.reference_ncopy = allele.reference_ncopy
        self.allele_ncopy = allele.allele_ncopy
        self.al_idx = allele.al_idx
        self.support = callers
        self.exp_flag = allele.exp_flag
    
    def add_support(self, callers):
        for caller in callers:
            if caller not in self.support:
                self.support.append(caller)

class ConnectedComponent:
    """
    Subgraph corresponding to nodes mapped to a single allele
    """
    def __init__(self, ccid, subgraph):
        self.cc_id = ccid
        self.subgraph = subgraph
        self.uniq_callers = self.GetUniqueCallers()
        self.caller_to_nodes = self.GetCallerToNodes()
        self.resolved_prealleles = self.GetResolvedPreAlleles()

    def GetUniqueCallers(self):
        uniq_callers = set()
        for node in self.subgraph.nodes():
            uniq_callers.add(node.GetVCFType())
        return uniq_callers

    def GetCallerToNodes(self):
        caller_to_nodes = {}
        for node in self.subgraph.nodes():
            if node.GetVCFType() not in caller_to_nodes:
                caller_to_nodes[node.GetVCFType()] = [node]
            else:
                caller_to_nodes[node.GetVCFType()].append(node)
        return caller_to_nodes

    def GetResolvedPreAlleles(self):
        resolved_prealleles = {}
        # If number of nodes == number of callers: 1-1-1
        if len(self.uniq_callers) == len(self.subgraph.nodes()):
            if trh.VcfTypes.hipstr in self.uniq_callers:
                tmp_node = self.caller_to_nodes[trh.VcfTypes.hipstr][0]
            else:
                tmp_node = list(self.subgraph.nodes())[0]
            pa = PreAllele(tmp_node, self.uniq_callers)
            resolved_prealleles['any'] = pa

        # If not 1-1-1    
        else:
            # This should only happen if hipSTR node exists
            # If that's not true, something bad happened...
            if trh.VcfTypes.hipstr not in self.uniq_callers:
                raise ValueError("HipSTR doesn't exist and we have a discrepancy!")

            # only one hipstr node
            if len(self.caller_to_nodes[trh.VcfTypes.hipstr]) == 1:
                tmp_node = self.caller_to_nodes[trh.VcfTypes.hipstr][0]
                pa = PreAllele(tmp_node, [trh.VcfTypes.hipstr])
                for node in self.subgraph:
                    if node != tmp_node and node.allele_sequence == tmp_node.allele_sequence:
                        pa.add_support([node.GetVCFType()])
                resolved_preallele['any'] = pa

            # More than one hipstr node: we need to assign 
            # different hipstr nodes for different allele idx
            else:
                tmp_node = None
                for node1 in self.subgraph:
                    if node1.GetVCFType() == trh.VcfTypes.hipstr:
                        if node1.al_idx not in resolved_prealleles:
                            tmp_node = node1
                            pa = PreAllele(tmp_node, [trh.VcfTypes.hipstr])
                            for node2 in self.subgraph:
                                if node2 != tmp_node and node2.allele_sequence == tmp_node.allele_sequence:
                                    pa.add_support([node2.GetVCFType()])
                            resolved_prealleles[node1.al_idx] = pa  
        return resolved_prealleles

class ClusterGraph:
    """
    Keeps track of graph of alleles called by each method

    Parameters
    ----------
    record_cluster : RecordCluster
       record cluster to build a graph out of 
    """

    def __init__(self, record_cluster):
        self.rclust = record_cluster
        self.graph = self.BuildGraph()

        # Get list of connected components
        # Sorted by the number of nodes in each
        sorted_ccs = sorted(nx.algorithms.components.connected_components(self.graph), \
                    key=len, reverse=True)
        self.connected_comps = [ConnectedComponent(CC_PREFIX+str(i), \
                                        self.graph.subgraph(sorted_ccs[i]).copy()) \
                                    for i in range(len(sorted_ccs))]

    def BuildGraph(self):
        r"""
        Build the allele graph
        """
        allele_list = self.GetAlleleList()
        graph = nx.Graph()

        # Add nodes - one for each allele
        for al in allele_list:
            graph.add_node(al)

        # Add edges between nodes of equivalent size
        # from different methods
        for nd1 in graph.nodes():
            for nd2 in graph.nodes():
                if nd1 == nd2:
                    continue
                if nd1.allele_size == nd2.allele_size \
                        and not graph.has_edge(nd1, nd2) \
                        and nd1.GetVCFType() != nd2.GetVCFType():
                    graph.add_edge(nd1, nd2)
        return graph

    def GetAlleleList(self):
        """
        Get list of Alleles called at least
        once across all RecordObj objects
        """
        alist = []
        for ro in self.rclust.record_objs:
            called_alleles = ro.GetCalledAlleles()
            for al_idx in called_alleles:
                alist.append(Allele(ro, al_idx))
        return alist

    def GetNodeObject(self, vcf_type, al_idx):
        for allele in self.graph.nodes:
            if allele.GetVCFType() == vcf_type and allele.al_idx == al_idx:
                return allele
        return None

    def GetSubgraphIndexForNode(self, node):
        if node is None:
            return None
        for i in range(len(self.connected_comps)):
            if node in self.connected_comps[i].subgraph:
                return i
        return None

    def GetSubgraphSize(self, ccid):
        if ccid < 0 or ccid >= len(self.connected_comps):
            return None
        return len(self.connected_comps[ccid].subgraph.nodes())
        


class RecordResolver:
    """
    Main class to resolve info for a record cluster

    Parameters
    ----------
    rc : RecordCluster
       the record cluster being resolved

    Attributes
    ----------
    rc_graph : ClusterGraph
       keeps track of alleles across records being merged
    resolved : bool
       Set to True once the record cluster has been resolved
    resolution_score : dict (str: float)
       Key=sample, Value=resolution score
       (i.e. all callers agreed)
    """
    def __init__(self, rc):
        self.record_cluster = rc
        self.rc_graph = ClusterGraph(rc)
        self.resolved = False

        # Get set after resolving
        self.resolved_prealleles = {}
        self.resolution_score = {}
        self.allele_support = {}
        self.resolution_method = {}
        self.empty_call = {}
        self.ref = None
        self.alts = []
        self.sample_to_info = {} # sample -> GT, NCOPY
        self.nocall = False

    def Resolve(self):
        resolved_prealleles = {}
        resolution_score = {}
        resolution_methods = {}
        allele_supports = {}
        for sample in self.record_cluster.samples:
            samp_call = self.record_cluster.GetSampleCall(sample)
            samp_qual_scores = self.record_cluster.GetQualScore(sample)
            resolved_connected_comp_ids, resolved_methods, score, allele_support = \
                self.GetConnectedCompForSingleCall(samp_call, samp_qual_scores)
            resolution_score[sample] = score
            resolution_methods[sample] = resolved_methods
            resolved_prealleles[sample] = self.ResolveSequenceForSingleCall(resolved_connected_comp_ids, samp_call)
            allele_supports[sample] = allele_support
        self.resolved_prealleles = resolved_prealleles
        self.resolution_score = resolution_score
        self.allele_support = allele_supports
        self.resolution_method = resolution_methods
        self.update()
        self.resolved = True
        return self.resolved
   
    def update(self):
        # First update alleles list
        for sample in self.resolved_prealleles:
            for pa in self.resolved_prealleles[sample]:
                if self.ref is None:
                    self.ref = pa.reference_sequence
                if pa.allele_sequence != self.ref and pa.allele_sequence != pa.reference_sequence:
                    if pa.allele_sequence not in self.alts:
                        if pa.allele_sequence != "":
                            self.alts.append(pa.allele_sequence)

        if self.ref is None:
            self.nocall = True 
        # Now update other info. need all alts for this
        for sample in self.resolved_prealleles:
            self.empty_call[sample] = False
            self.sample_to_info[sample] = {}
            GT_list = []
            GB_list = []
            NCOPY_list = []
            Expanded = []
            for pa in self.resolved_prealleles[sample]:
                if pa.exp_flag:
                    Expanded.append("1")
                else:
                    Expanded.append("0")
                if pa.al_idx != 0 and pa.allele_sequence != self.ref:
                    if pa.allele_sequence == "":
                        self.empty_call[sample] = True
                        break
                    GT_list.append(str(self.alts.index(pa.allele_sequence) + 1))
                    GB_list.append(str(len(pa.allele_sequence) - len(self.ref)))
                    NCOPY_list.append(str(pa.allele_ncopy))
                else:
                    GT_list.append('0')
                    GB_list.append('0')
                    NCOPY_list.append(str(pa.reference_ncopy))
            if len(GT_list) == 0 or self.empty_call[sample]:
                GT_list = ['.']
                NCOPY_list = ['.']
                GB_list = ['.']
                Expanded = ['.']
            self.sample_to_info[sample]["GT"] = '/'.join(GT_list)
            self.sample_to_info[sample]["GB"] = '/'.join(GB_list)
            self.sample_to_info[sample]["NCOPY"] = ','.join(NCOPY_list)
            self.sample_to_info[sample]["EXP"] = '/'.join(Expanded)

    def GetSampleScore(self, sample):
        if self.resolution_scoreGT[sample] == -1 or self.empty_call[sample]:
            return "."
        return str(self.resolution_score[sample])

    def GetSampleGTS(self, sample):
        if len(self.resolution_method[sample]) == 0 or self.empty_call[sample]:
            return "."
        return '|'.join([str(method) for method in self.resolution_method[sample]])

    def GetSampleALS(self, sample):
        if not self.allele_support[sample] or self.empty_call[sample]:
            return "."
        return ",".join([str(key) + "|" + str(val) for key,val in self.allele_support[sample].items()])

    def GetSampleGT(self, sample):
        return self.sample_to_info[sample]["GT"]

    def GetSampleGB(self, sample):
        return self.sample_to_info[sample]["GB"]

    def GetSampleNCOPY(self, sample):
        return self.sample_to_info[sample]["NCOPY"]

    def GetExpandedFlag(self, sample):
        return self.sample_to_info[sample]['EXP']

    def TestScore(self, score):
        if np.isnan(score):
            return False
        if score < 0 or score > 1:
            return False
        return True

    def GetConnectedCompForSingleCall(self, samp_call, samp_qual_scores):
        r"""

        Parameters
        ----------
        samp_call: dict(vcftype:allele)
                allele: [al1, al2, BOOL] or [-1, BOOL] for no calls

        Returns
        -------
        ret_cc_ids : list of int
           indices of CCs of resolved call
        score : float
           Score of caller agreement
        sup_method : list
           Supporting method for the ret_cc_ids
        allele_size_support: dict
            Key=bp diff of allele and the ref genome, Value=number of times we saw this allele.
        """
        method_cc = defaultdict(list)  # methods supporting each pair of CCID
        first_allele = []
        second_allele = []
        method_dict = {"advntr": [1, 0, 0, 0], "eh": [0, 1, 0, 0], "hipstr": [0, 0, 1, 0], "gangstr": [0, 0, 0, 1]}
        allele_size_support = {}
        for method in samp_call:
                # check for no calls
                if samp_call[method][0] == -1:
                        continue  # no call
                # Get the IDs of supported connected components
                ccids = []
                for i in [0, 1]:
                        node = self.rc_graph.GetNodeObject(method, samp_call[method][i])
                        ccid = self.rc_graph.GetSubgraphIndexForNode(node)
                        ccids.append(ccid)
                        allele_size_support[node.allele_size] = allele_size_support.get(node.allele_size, 0) + 1

                ccids.sort()
                ccids = (ccids[0],ccids[1])
                method_cc[ccids].append(method)
                if ccids[0] == ccids[1]:
                    first_allele.append((ccids[0], method))
                    second_allele.append((ccids[1],method))
                else:
                    first_allele.append((ccids[0], method))
                    first_allele.append((ccids[1], method))

        if len(method_cc) == 0:  # no call across all methods
            return [],[], -1, -1, {}
        ret_cc_ids, score = self.ResolveScore(samp_qual_scores, method_cc, first_allele, second_allele)
        assert(self.TestScore(score))

        ret_sup_methods = [method.value for method in method_cc[ret_cc_ids]]
        sup_method = [0, 0, 0, 0]
        for method in ret_sup_methods:
                sup_method = [sum(x) for x in zip(sup_method, method_dict[method])] 
        return list(ret_cc_ids), sup_method, round(score, 2), allele_size_support



    def ResolveScore(self, samp_qual_scores, method_cc, first_allele, second_allele):
        r"""

        Parameters
        ----------
        samp_qual_scores : (dict of str: str)
           Key=VCF type, Value=sample quality score

        method_cc : (dict of tuple)
           Key=pair of ccids, Value, list of supporting methods
        Returns
        -------
        ret_cc_ids : list of int
           indices of CCs of resolved call
        max_score : float
           score of ret_cc_ids
        """

        scores = {}
        sum_scores = 0
        max_seen_score = 0
        for pair in method_cc:
            score = 0
            for method in method_cc[pair]:
                score += samp_qual_scores[method]
                if samp_qual_scores[method] > max_seen_score:
                    max_seen_score = samp_qual_scores[method]
            sum_scores += score
            scores[pair] = score
        if sum_scores != 0:
            for pair in scores:
                scores[pair] = scores[pair]/sum_scores
        else:
            scores[pair] = 0

        max_score = max(scores.values())
        max_pair = []
        for pair in scores:
            if scores[pair] == max_score:
                max_pair.append(pair)

        if len(max_pair) == 1:
            pair = max_pair[0]
            return pair, max_score * max_seen_score

        else:  # ties
            for method in ['hipstr', 'gangstr', 'eh', 'advntr']:  # break ties with giving priority to methods
                for pair in max_pair:
                    for method_ in method_cc[pair]:
                        if method_.value == method:
                            return pair, max_score * max_seen_score



        
    def ResolveSequenceForSingleCall(self, ccid_list, samp_call):
        if len(ccid_list) == 0: return []
        pre_allele_list = []
        for cc_id in ccid_list:
            connected_comp = self.rc_graph.connected_comps[cc_id]
            resolved_prealleles = connected_comp.resolved_prealleles
            if "any" in resolved_prealleles:
                pre_allele_list.append(resolved_prealleles["any"])
                continue
            elif trh.VcfTypes.hipstr in connected_comp.uniq_callers and \
                trh.VcfTypes.hipstr in samp_call and samp_call[trh.VcfTypes.hipstr][0] != -1:
                # cc has hipstr nodes, and we have hipstr calls
                hip_call = samp_call[trh.VcfTypes.hipstr]
                for al_idx in hip_call[0:2]:
                    if al_idx in resolved_prealleles:
                        pre_allele_list.append(resolved_prealleles[al_idx])

            else:
                # TODO. For now just return a random allele if we don't have hipstr
                pre_allele_list.append(list(resolved_prealleles.values())[0])

        pre_allele_list = list(set(pre_allele_list))
        if len(pre_allele_list) == 1:
            pre_allele_list = [pre_allele_list[0], pre_allele_list[0]]
        if len(pre_allele_list) > 2:
            n_copies = []
            for i in range(len(pre_allele_list)):
                n_copies.append(pre_allele_list[i].allele_ncopy)
            n_copies.sort()
            if n_copies[0] == n_copies[1] or n_copies[1] == n_copies[2]:
                for pa in pre_allele_list:
                    if pa.allele_ncopy == n_copies[1]:
                        pre_allele_list.remove(pa)
                        return pre_allele_list

        return pre_allele_list    
