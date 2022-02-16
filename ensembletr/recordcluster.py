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
    def __init__(self, rec, vcf_type):
        self.cyvcf2_record = rec
        self.vcf_type = vcf_type
        self.hm_record = trh.HarmonizeRecord(vcf_type, rec)
        self.canonical_motif = GetCanonicalMotif(self.hm_record.motif)
        self.prepend_seq = ''
        self.append_seq = ''

    def GetSamples(self):
        return self.samples

    def GetCalledAlleles(self):
        al_idx = set()
        for call in self.hm_record.vcfrecord.genotypes:
            if call[0] != -1 and len(call) == 3:
                al_idx.add(call[0])
                al_idx.add(call[1])
        return al_idx

    def GetROSampleCall(self, samp_idx):
        return self.cyvcf2_record.genotypes[samp_idx]

    def GetSampleString(self, samp_idx):
        samp_call = self.GetROSampleCall(samp_idx)
        if samp_call is None or samp_call[0] == -1:
            sampdata = "."
        else:
            call_ncopy_list = []
            for idx in samp_call[0:2]:
                if idx == 0:
                    call_ncopy_list.append(str(self.hm_record.ref_allele_length))
                else:
                    call_ncopy_list.append(str(self.hm_record.alt_allele_lengths[idx - 1]))
            sampdata = ",".join(call_ncopy_list)
        callstr = "%s=%s"%(self.vcf_type.name, sampdata)
        return callstr

    def GetScores(self, samp_idx):
        vcf_idx = convert_type_to_idx[self.vcf_type]
        if vcf_idx == 0:
            return self.cyvcf2_record.format('ML')[samp_idx][0]
        if vcf_idx == 1:
            REPCI = self.cyvcf2_record.format('REPCI')[samp_idx]
            REPCN = self.cyvcf2_record.format('REPCN')[samp_idx]
            if REPCI == "." or REPCN == ".":
                return 0
            length = len(self.cyvcf2_record.INFO['RU'])
            return self.ehScore(REPCI, REPCN, length)
        if vcf_idx == 2 or vcf_idx == 3:
            return self.cyvcf2_record.format('Q')[samp_idx][0]

    def ehScore(self, conf_invs, CNs, length):
        conf_invs = conf_invs.split("/")
        CNs = CNs.split("/")
        CNs = [(int(CN) * length) for CN in CNs]
        score1 = self.CalcScore(conf_invs[0], CNs[0])
        score2 = self.CalcScore(conf_invs[1], CNs[1])
        return 0.8 * min(score1, score2) + 0.2 * max(score1, score2)

    def CalcScore(self, conf_inv, allele):
        conf_inv = conf_inv.split("-")
        dist = abs(int(conf_inv[0]) - int(conf_inv[1]))
        if dist > 100:
             return 0
        if allele == 0:
             return 1/math.exp(4 * (dist))
        return 1/math.exp(4 * (dist) / int(allele))
        

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
                rec.prepend_str = self.fasta[chrom][self.first_pos : rec.cyvcf2_record.POS].seq.upper()
            
            if rec.cyvcf2_record.end < self.last_end:
                # Found a record that ends before last end
                # Should append the record
                rec.append_str = self.fasta[chrom][rec.cyvcf2_record.end : self.last_end].seq.upper()  

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
            samp_call_list = [ro.GetSampleString(self.samples.index(sample)) for ro in self.record_objs]
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
            ret_dict[rec.vcf_type] = rec.GetROSampleCall(self.samples.index(sample))
        return ret_dict

    def GetSampleScore(self, sample):
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
            ret_dict[rec.vcf_type] = rec.GetScores(self.samples.index(sample))
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

        alleles = [ro.hm_record.ref_allele]+ro.hm_record.alt_alleles
        allele_lengths = [ro.hm_record.ref_allele_length] + ro.hm_record.alt_allele_lengths

        self.reference_sequence = ro.prepend_seq + alleles[0] + ro.append_seq
        self.allele_sequence = ro.prepend_seq + alleles[al_idx] + ro.append_seq
        self.allele_size = len(self.allele_sequence) - len(self.reference_sequence)
        self.allele_ncopy = allele_lengths[al_idx]
        self.reference_ncopy = ro.hm_record.ref_allele_length

    def GetVCFType(self):
        return self.record_object.vcf_type

class PreAllele:
    def __init__(self, allele, callers):
        self.reference_sequence = allele.reference_sequence
        self.allele_sequence = allele.allele_sequence
        self.reference_ncopy = allele.reference_ncopy
        self.allele_ncopy = allele.allele_ncopy
        self.support = callers
    
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
                self.resolved_preallele['any'] = pa

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
            samp_qual_scores = self.record_cluster.GetSampleScore(sample)
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
                if pa.allele_sequence != pa.reference_sequence:
                    if pa.allele_sequence not in self.alts:
                        self.alts.append(pa.allele_sequence)

        if self.ref is None:
            self.nocall = True 
        # Now update other info. need all alts for this
        for sample in self.resolved_prealleles:
            self.sample_to_info[sample] = {}
            GT_list = []
            NCOPY_list = []
            for pa in self.resolved_prealleles[sample]:
                if pa.allele_sequence != self.ref:
                    GT_list.append(str(self.alts.index(pa.allele_sequence) + 1))
                    NCOPY_list.append(str(pa.allele_ncopy))
                else:
                    GT_list.append('0')
                    NCOPY_list.append(str(pa.reference_ncopy))
            if len(GT_list) == 0:
                GT_list = ['.']
                NCOPY_list = ['.']
            self.sample_to_info[sample]["GT"] = '/'.join(GT_list)
            self.sample_to_info[sample]["NCOPY"] = ','.join(NCOPY_list)

    def GetSampleScore(self, sample):
        if self.resolution_score[sample] == -1:
            return "."
        return str(self.resolution_score[sample])

    def GetSampleGTS(self, sample):
        if len(self.resolution_method[sample]) == 0:
            return "."
        return '|'.join(self.resolution_method[sample])

    def GetSampleALS(self, sample):
        if not self.allele_support[sample]:
            return "."
        return ",".join([str(key) + "|" + str(val) for key,val in self.allele_support[sample].items()])

    def GetSampleGT(self, sample):
        return self.sample_to_info[sample]["GT"]

    def GetSampleNCOPY(self, sample):
        return self.sample_to_info[sample]["NCOPY"]

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

        TODO: work on logic
        Want to instead determine which caller most reliable
        Or alternatively vote
        Then return its call, plus support from other callers
        Instead of "certain", return a score
        """
        method_cc = defaultdict(list) # methods supporting each CCID
        allele_size_support = {}
        for method in samp_call:
                # check for no calls
                if samp_call[method][0] == -1:
                        continue # no call
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
        
        scores = {}
        sum_scores = 0
        for pair in method_cc:
                score = 0
                for method in method_cc[pair]:
                        score += samp_qual_scores[method]
                sum_scores += score
                scores[pair] = score
        if sum_scores != 0:
                for pair in scores:
                        scores[pair] = scores[pair]/sum_scores
        ret_cc_ids = max(scores, key=scores.get, default = -1) # get alleles with maximum score
        
        if ret_cc_ids == -1:
                return [],[], -1, {}

        ret_sup_methods = [method.value for method in method_cc[ret_cc_ids]]
        method_dict = {"advntr":[1,0,0,0],"eh":[0,1,0,0],"hipstr":[0,0,1,0],"gangstr":[0,0,0,1]}
        sup_method = [0,0,0,0]
        for method in ret_sup_methods:
                sup_method = [sum(x) for x in zip(sup_method, method_dict[method])] 
        sup_method = [str(method) for method in sup_method]
        return list(ret_cc_ids), sup_method, round(scores[ret_cc_ids],2), allele_size_support



        
    def ResolveSequenceForSingleCall(self, ccid_list, samp_call):
        if len(ccid_list) == 0: return []
        pre_allele_list = []
        for cc_id in ccid_list:
            connected_comp = self.rc_graph.connected_comps[cc_id]
            resolved_prealleles = connected_comp.resolved_prealleles
            if "any" in resolved_prealleles:
                pre_allele_list.append(resolved_prealleles["any"])
                continue
            if trh.VcfTypes.hipstr in connected_comp.uniq_callers and \
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
        return pre_allele_list    
