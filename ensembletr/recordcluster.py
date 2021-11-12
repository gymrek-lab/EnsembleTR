"""
RecordCluster object includes mergeable records.
Mergeable records are defined as records that overlap and share a motif
"""

import trtools.utils.tr_harmonizer as trh
from trtools.utils.utils import GetCanonicalMotif
from enum import Enum
import networkx as nx
import numpy as np

CC_PREFIX = 'cc'

class AlleleType(Enum):
    Reference = 0
    Alternate = 1

ColorDict = {trh.VcfTypes.advntr: "tab:blue",
             trh.VcfTypes.eh: "tab:orange",
             trh.VcfTypes.hipstr: "tab:red",
             trh.VcfTypes.gangstr: "tab:green",
             }
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
    vcf_type: trh.TRRecordHarmonizer.vcftype
       Type of the VCF file
    VCF : cyvcf2.VCF
       VCF reader from which the object came
    samples : list of str
       List of samples included
    rec : cyvcf2.VCF.vcfrecord
       VCF record for the object
    """
    def __init__(self, vcf_type, VCF, samples, rec):
        self.cyvcf2_record = rec
        self.VCF = VCF
        self.samples = samples
        self.vcf_type = vcf_type
        self.hm_record = trh.HarmonizeRecord(vcf_type, rec)
        self.ref = self.hm_record.ref_allele
        self.motif = self.hm_record.motif
        self.canonical_motif = GetCanonicalMotif(self.motif)
        self.prepend_seq = ''
        self.append_seq = ''

    def GetSamples(self):
        return self.samples

    def GetVcfRegion(self):
        return str(self.cyvcf2_record.CHROM) + ':' + str(self.cyvcf2_record.POS)
    
    def GetROSampleCall(self, sample):
        samp_idx = self.samples.index(sample)
        return self.cyvcf2_record.genotypes[samp_idx]


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

    def __str__(self):
        return "Ref: " + self.ref_seq + "\nSeq: " + self.seq + "\nSupp: " + str(self.support)
    def __repr__(self):
        return "Ref: " + self.ref_seq + "\nSeq: " + self.seq + "\nSupp: " + str(self.support)
    
class Allele:
    def __init__(self, ro, atype, allele_sequence, diff_from_ref, genotype_idx):
        if atype not in [AlleleType.Reference, AlleleType.Alternate]:
            raise ValueError('Unknown allele type: ' + atype)
        self.allele_type = atype
        self.record_object = ro
        self.original_allele_sequence = allele_sequence
        self.original_reference_sequnce = ro.hm_record.ref_allele
        self.allele_sequence = ro.prepend_seq + allele_sequence + ro.append_seq
        self.reference_sequence = ro.prepend_seq + ro.hm_record.ref_allele + ro.append_seq
        # self.allele_ncopy = ro.
        if genotype_idx == 0:
            self.allele_ncopy = ro.hm_record.ref_allele_length
        else:
            self.allele_ncopy = ro.hm_record.alt_allele_lengths[genotype_idx - 1]
        self.reference_ncopy = ro.hm_record.ref_allele_length
        self.allele_size = diff_from_ref
        self.genotype_idx = genotype_idx
        self.vcf_type = ro.vcf_type

    def GetLabel(self):
        if self.allele_type == AlleleType.Reference:
            return self.vcf_type.name + '_*'
        else:
            return self.vcf_type.name + '_' + str(self.allele_size)



class RecordCluster:
    def __init__(self, recobjs, ref_genome):
        self.motif = recobjs[0].motif
        self.canonical_motif = GetCanonicalMotif(self.motif)
        self.vcf_types = [False] * len(convert_type_to_idx.keys())
        self.samples = recobjs[0].samples
        self.first_pos = -1
        self.last_end = -1
        # References used for prepending and appending (can be different)
        self.fasta = ref_genome
        
        # First, find the first start and last end
        for rec in recobjs:
            if self.first_pos == -1:
                self.first_pos = rec.cyvcf2_record.POS
            if self.last_end == -1:
                self.last_end = rec.cyvcf2_record.end
            if rec.cyvcf2_record.POS < self.first_pos:
                self.first_pos = rec.cyvcf2_record.POS
            if rec.cyvcf2_record.end > self.last_end:
                self.last_end = rec.cyvcf2_record.end

            self.vcf_types[convert_type_to_idx[rec.vcf_type]] = True
            if rec.canonical_motif != self.canonical_motif:
                raise ValueError('RecordObjs in a cluster must share canonical motif.\n' +
                                 'Cannot add RO with ' + rec.canonical_motif + ' canonical motif to cluster' +
                                 'with ' + self.canonical_motif + 'canonical motif.')

        self.record_objs = recobjs
        for rec in self.record_objs:
            chrom = rec.cyvcf2_record.CHROM
            if rec.cyvcf2_record.POS > self.first_pos:
                # Found a record that starts after
                # Should prepend the record
                rec.prepend_str = self.fasta[chrom][self.first_pos : rec.cyvcf2_record.POS].seq.upper()
            
            if rec.cyvcf2_record.end < self.last_end:
                # Found a record that ends before last end
                # Should append the record
                rec.append_str = self.fasta[chrom][rec.cyvcf2_record.end : self.last_end].seq.upper()  

    def GetAllInputs(self):
        out_dict = {}
        for sample in self.samples:
            samp_call_list = []
            for ro in self.record_objs:
                temp_str = ro.vcf_type.name + '='
                
                samp_call = ro.GetROSampleCall(sample)
                if samp_call is None or samp_call[0] == -1:
                    temp_str = temp_str + '.'
                else:
                    call_idxs = samp_call[0:2]
                    call_ncopy_list = []
                    for idx in call_idxs:
                        if idx == 0:
                            call_ncopy_list.append(str(ro.hm_record.ref_allele_length))
                        else:
                            call_ncopy_list.append(str(ro.hm_record.alt_allele_lengths[idx - 1]))
                    temp_str = temp_str + ','.join(call_ncopy_list)
                samp_call_list.append(temp_str)
            out_dict[sample] = '|'.join(samp_call_list)
        return out_dict
        

    def AppendRecordObject(self, ro):
        if self.canonical_motif != ro.canonical_motif:
            raise ValueError('Canonical motif of appended record object mush match record cluster.')

        self.record_objs.append(ro)
        self.vcf_types[convert_type_to_idx[ro.vcf_type]] = True

        # Update first pos and last end and appropriate references for prepending and appending
        if ro.cyvcf2_record.POS < self.first_pos:
            self.first_pos = ro.cyvcf2_record.POS
        if ro.cyvcf2_record.end > self.last_end:
            self.last_end = ro.cyvcf2_record.end

        # Update prepend and append strings:
        for rec in self.record_objs:
            chrom = rec.cyvcf2_record.CHROM
            if rec.cyvcf2_record.POS > self.first_pos:
                # Found a record that starts after
                # Should prepend the record
                rec.prepend_seq = self.fasta[chrom][self.first_pos : rec.cyvcf2_record.POS].seq.upper()
            
            if rec.cyvcf2_record.end < self.last_end:
                # Found a record that ends before last end
                # Should append the record
                rec.append_seq = self.fasta[chrom][rec.cyvcf2_record.end : self.last_end].seq.upper()

    def GetVcfTypesTuple(self):
        return tuple(self.vcf_types)

    def GetAlleleList(self):
        alist = []
        for ro in self.record_objs:
            ref = ro.hm_record.ref_allele
            # find all genotype indexes (to avoid adding alleles that are not actually called)
            # Example: Hipstr has lines with REF ALT . . . . (all no calls)
            # This list of available genotypes ensures we only add alleles that are 
            # represented in the calls
            
            genot_set = set()
            for call in ro.hm_record.vcfrecord.genotypes:
                # check if no call
                if call[0] != -1 and len(call) == 3:
                    genot_set.add(call[0])
                    genot_set.add(call[1])
            if 0 in genot_set:
                alist.append(Allele(ro, AlleleType.Reference, ref, 0, 0))
            altnum = 1
            for alt in ro.hm_record.alt_alleles:
                if altnum in genot_set:
                    alist.append(Allele(ro, AlleleType.Alternate, alt, len(alt) - len(ref), altnum))
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
        self.ccid_to_connected_comp_dict = {}
        self.ccid_to_subgraph_dict = {}
        self.ccid_to_resolved_pre_allele = {}
        self.caller_to_nodes_dict = {}
        self.uniq_callers_dict = {}
        self.num_ccs = 0
        i = 0

        for cc in self.sorted_connected_comps:
            cc_id = CC_PREFIX + str(i)
            self.ccid_to_connected_comp_dict[cc_id] = cc
            self.ccid_to_subgraph_dict[cc_id] = self.graph.subgraph(cc).copy()
            # Update: For each cc_id, we have one resolved pre allele per genotype idx
            # If only one pre_allele, we have {'any': pa}
            self.ccid_to_resolved_pre_allele[cc_id] = {}

            uniq_callers = []
            caller_to_nodes = {}
            for node in self.ccid_to_subgraph_dict[cc_id].nodes():
                if node.vcf_type not in caller_to_nodes:
                    caller_to_nodes[node.vcf_type] = [node]
                else:
                    caller_to_nodes[node.vcf_type].append(node)
                
                if node.vcf_type not in uniq_callers:
                    uniq_callers.append(node.vcf_type)
            self.caller_to_nodes_dict[cc_id] = caller_to_nodes
            self.uniq_callers_dict[cc_id] = uniq_callers


            # If number of nodes == number of callers: 1-1-1
            if len(self.uniq_callers_dict[cc_id]) == len(self.ccid_to_subgraph_dict[cc_id].nodes()):
                if trh.VcfTypes.hipstr in self.uniq_callers_dict[cc_id]:
                    tmp_node = self.caller_to_nodes_dict[cc_id][trh.VcfTypes.hipstr][0]
                else:
                    tmp_node = list(self.ccid_to_subgraph_dict[cc_id].nodes())[0]
                pa = PreAllele(tmp_node, uniq_callers)
                self.ccid_to_resolved_pre_allele[cc_id]['any'] = pa
                
            # If not 1-1-1    
            else:
                # if hipstr exists (it should):
                if trh.VcfTypes.hipstr in self.uniq_callers_dict[cc_id]:
                    # only one hipstr node
                    if len(self.caller_to_nodes_dict[cc_id][trh.VcfTypes.hipstr]) == 1:
                        tmp_node = self.caller_to_nodes_dict[cc_id][trh.VcfTypes.hipstr][0]
                        pa = PreAllele(tmp_node, [trh.VcfTypes.hipstr])
                        for node in self.ccid_to_subgraph_dict[cc_id]:
                            if node != tmp_node and node.allele_sequence == tmp_node.allele_sequence:
                                pa.add_support([node.vcf_type])
                        self.ccid_to_resolved_pre_allele[cc_id]['any'] = pa
                        
                    else:
                        # More than one hipstr node: we need to assign different hipstr nodes for different genotype idx
                        tmp_node = None
                        for allele in self.ccid_to_subgraph_dict[cc_id]:
                            if allele.vcf_type == trh.VcfTypes.hipstr:
                                if allele.genotype_idx not in self.ccid_to_resolved_pre_allele[cc_id]:
                                    tmp_node = allele
                                    pa = PreAllele(tmp_node, [trh.VcfTypes.hipstr])
                                    for node in self.ccid_to_subgraph_dict[cc_id]:
                                        if node != tmp_node and node.allele_sequence == tmp_node.allele_sequence:
                                            pa.add_support([node.vcf_type])
                                    self.ccid_to_resolved_pre_allele[cc_id][allele.genotype_idx] = pa  
                else:
                    raise ValueError("HipSTR doesn't exist and we have a discrepancy!")

            i+=1
        self.subgraphs = list(self.ccid_to_subgraph_dict.values())

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
    def GetSubgraphIDForNode(self, node):
        if node is None:
            return None
        for cc_id in self.ccid_to_subgraph_dict:
            if node in self.ccid_to_subgraph_dict[cc_id].nodes():
                return cc_id
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
        # Resolve sequence for all connected components
        # for component in self.rc_graph.GetSortedConnectedComponents():
        #     num_unique_callers_in_component = 0
        #     callers_seen = []
        #     for allele in component:
        #         if allele.vcf_type not in callers_seen:
        #             callers_seen.append(allele.vcf_type)
        #             num_unique_callers_in_component += 1
        #     list_unique_caller_nodes_in_conn_comp.append(num_unique_callers_in_component)

    def GetConnectedCompForSingleCall(self, samp_call):
        r"""

        Parameters
        ----------
        samp_call: dict(vcftype:allele)
                allele: [al1, al2, BOOL] or [-1, BOOL] for no calls

        Returns
        -------
        list of connected component ids corresponding to resolved call
        """
        # Assign connected component ids to alleles
        conn_comp_cc_id_dict = {}
        cc_id_support = {}
        num_valid_methods = 0
        for method in samp_call:
            # check for no calls
            if samp_call[method][0] == -1:
                cc_id0 = None
                cc_id1 = None
            else:
                cc_id0 = self.rc_graph.GetSubgraphIDForNode(self.rc_graph.GetNodeObject(method, samp_call[method][0]))
                cc_id1 = self.rc_graph.GetSubgraphIDForNode(self.rc_graph.GetNodeObject(method, samp_call[method][1]))
                for cc_id in [cc_id0, cc_id1]:
                    if cc_id not in cc_id_support:
                        # idea: instead of 1, locus_reliability(method) * call_reliability
                        cc_id_support[cc_id] = 1
                    else:
                        cc_id_support[cc_id] += 1
                num_valid_methods += 1
            if cc_id0 is not None and cc_id1 is not None:
                conn_comp_cc_id_dict[method] = [cc_id0, cc_id1]

        sorted_ccid_support = dict(sorted(cc_id_support.items(), key=lambda item: item[1], reverse=True))
        ret_cc_ids = []
        certain = True
        for cc_id in sorted_ccid_support:
            # If all alleles in all methods point to one cc
            # for Hom require support to be perfect (2x num_valide) if not, report not certain
            # for Het do the same for support two ccs each with num_valid
            if sorted_ccid_support[cc_id] == num_valid_methods * 2: 
                return [cc_id, cc_id], certain
            elif sorted_ccid_support[cc_id] < num_valid_methods * 2 and sorted_ccid_support[cc_id] > num_valid_methods:
                certain = False
                return [cc_id, cc_id], certain
            elif sorted_ccid_support[cc_id] == num_valid_methods:
                if len(ret_cc_ids) < 2:
                    ret_cc_ids.append(cc_id)
                else:
                    certain = False
                    # We have already added 2 top supported ccs, but we still have another cc
                    print("WARNING: extra cc", cc_id)
                    
            else:
                certain = False
                print("WARNING: at least one CC is not fully supported by either one or both alleles", cc_id)
        # print (certain)
        return ret_cc_ids, certain
        


    def ResolveSequenceForSingleCall(self, ccid_list, samp_call):
        # Next update: Melissa's idea --> implemented!
        # Pre resolve possible sequences for each CC. Only check if samp_call contains the caller and genot_idx of the resolved seq
        # Add support based on overlap between samp_call and methods in cc
        pre_allele_list = []
        all_no_calls = True
        for vcftyp in samp_call:
            if samp_call[vcftyp][0] != -1:
                all_no_calls = False
        if all_no_calls:
            return pre_allele_list
        if len(ccid_list) == 0:
            return pre_allele_list
        for cc_id in ccid_list:
            resolved_ccid = False
            if cc_id in self.rc_graph.ccid_to_resolved_pre_allele:
                # Check if any caller points to same pa (1-1-1):
                if 'any' in self.rc_graph.ccid_to_resolved_pre_allele[cc_id]:
                    pre_allele_list.append(self.rc_graph.ccid_to_resolved_pre_allele[cc_id]['any'])
                else:
                    # We use hipster call to find the correct allele
                    # If hipstr calls exist in this cc_id, we can use them to find the correct pa
                    if trh.VcfTypes.hipstr in self.rc_graph.uniq_callers_dict[cc_id]:
                        # this cc has hipstr nodes
                        if trh.VcfTypes.hipstr in samp_call and samp_call[trh.VcfTypes.hipstr][0] != -1:
                            # cc has hipstr nodes, and we have hipstr calls
                            hip_call = samp_call[trh.VcfTypes.hipstr]
                            # find the matching genotype_idx
                            for genot_idx in hip_call[0:1]:
                                if genot_idx in self.rc_graph.ccid_to_resolved_pre_allele[cc_id]:
                                    pre_allele_list.append(self.rc_graph.ccid_to_resolved_pre_allele[cc_id][genot_idx])
                                    resolved_ccid = True
                                    break
                        else:
                            # samp_call doesn't have hipstr, but the resolved cc_id points to a cc that has hipstr calls
                            # (we can't resolve correctly)
                            print('Warning: samp_call does not have hipstr, but the resolved cc_id points to a cc that has hipstr calls ')
                    else:
                        # This cc doesn't have hipstr calls, we should resolve using whatever caller exists:
                        # TODO implement
                        raise ValueError("AAA lol")
            else:
                raise ValueError("No resolved sequence for cc_id: ", cc_id)
            

        if len(pre_allele_list) >= 2:
            return pre_allele_list[0:2]
        elif len(pre_allele_list) == 1:
            return [pre_allele_list[0], pre_allele_list[0]]
        else:
            # import networkx as nx
            # import matplotlib
            # matplotlib.use('TkAgg')
            # import matplotlib.pyplot as plt

            # pos = nx.spring_layout(self.rc_graph.graph, k=2 / np.sqrt(len(self.rc_graph.graph.nodes)))
            # nx.draw(self.rc_graph.graph, pos, node_color=self.rc_graph.colors)
            # nx.draw_networkx_labels(self.rc_graph.graph, pos, labels=self.rc_graph.labels)

            # plt.show()
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

