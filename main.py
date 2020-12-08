#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
python callsetmerger.py --vcfs hipstr.chr21.sorted.vcf.gz,advntr.chr21.sorted.vcf.gz,gangstr.chr21.sorted.vcf.gz
"""

import argparse
from callsetmerger.callsetmerger import Readers
from callsetmerger.recordcluster import ClusterGraph
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcfs", help="Comma-separated list of VCFs to merge. Must be sorted/indexed", type=str,
                        required=True)
    args = parser.parse_args()

    readers = Readers(args.vcfs.split(","))

    # Check samples same in each VCF
    # TODO

    # Walk through sorted readers
    while not readers.done:
        # Get mergeable calls
        ov_region = readers.getMergableCalls()
        rc_list = ov_region.RecordClusters
        if rc_list is None:
            # Move on
            readers.goToNext()
            continue

        for rc in rc_list:
            allele_list = rc.GetAlleleList()
            cg = ClusterGraph(allele_list)
            nx.layout
            pos = nx.spring_layout(cg.graph, k=2 / np.sqrt(len(cg.graph.nodes)))
            nx.draw(cg.graph, pos, node_color=cg.colors)
            # for p in pos:  # push text to right
            #     pos[p][0] += 0.2
            nx.draw_networkx_labels(cg.graph, pos, labels=cg.labels)
            # nx.draw_networkx(cg.graph, node_color=cg.colors, labels=cg.labels, seed=100)

            plt.show()
            pass

        # Merge calls
        # mo = MergeObject(allele_list, sample_calls)

        # TODO will need to update is_min_pos_list if we merged something

        # Move on
        readers.goToNext()
        input('Press any key to move on to next record.\n')


if __name__ == "__main__":
    main()
