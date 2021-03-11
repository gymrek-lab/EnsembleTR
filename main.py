#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
python callsetmerger.py --vcfs hipstr.chr21.sorted.vcf.gz,advntr.chr21.sorted.vcf.gz,gangstr.chr21.sorted.vcf.gz
"""
import traceback
import argparse
from callsetmerger.callsetmerger import Readers, GetWriter
from callsetmerger.recordcluster import ClusterGraph
from callsetmerger.mergemodule import RecordClusterMerger
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcfs", help="Comma-separated list of VCFs to merge. Must be sorted/indexed", type=str,
                        required=True)
    parser.add_argument("--outvcftemplate", help="(TODO make general tmp or remove if pyvcf) Template to use for output VCF.", type=str,
                        required=True)
    parser.add_argument("--out", help="Output merged VCF file", type=str, required= True)

    args = parser.parse_args()

    readers = Readers(args.vcfs.split(","))
    out_path = args.out
    # CyVCF2 needs a VCF template, I'm using HipSTR for first iteration
    template_path = args.outvcftemplate # NOT USED IN PYVCF OUTPUT, REMOVE

    # Check samples same in each VCF
    # TODO


    # Create VCF writer for output
    outvcf = GetWriter(out_path, readers.samples)
    i = 0
    # Walk through sorted readers
    while not readers.done:
        # Get mergeable calls
        ov_region = readers.getMergableCalls()
        rc_list = ov_region.RecordClusters
        if rc_list is None:
            # Move on
            # TODO Write first variant to VCF
            readers.goToNext()
            continue

        for rc in rc_list:
            num_vcfs = len([i for i in rc.vcf_types if i == True])
            if num_vcfs >= 2:
                cg = ClusterGraph(rc)
                # nx.layout()
                pos = nx.spring_layout(cg.graph, k=2 / np.sqrt(len(cg.graph.nodes)))
                nx.draw(cg.graph, pos, node_color=cg.colors)
                nx.draw_networkx_labels(cg.graph, pos, labels=cg.labels)

                # plt.show()
                
                # Merge calls
                mo = RecordClusterMerger(rc, readers.samples)
                # rec = mo.GetPyVCFRecord()
                rec = mo.GetRawVCFRecord()
                print(rec)
                print("Writing:")
                try:
                    outvcf.write(rec)
                except:
                    traceback.print_exc()
                

                i = i + 1
                # input('Press any key to move on to next record.\n')

        # TODO will need to update is_min_pos_list if we merged something

        # Move on
        readers.goToNext()
        if i == 3:
            break
    outvcf.close()

if __name__ == "__main__":
    main()
