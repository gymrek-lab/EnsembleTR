#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
python3 main.py --vcfs hipstr.chr21.sorted.vcf.gz,advntr.chr21.sorted.vcf.gz,gangstr.chr21.sorted.vcf.gz
"""

import argparse
from callsetmerger.callsetmerger import Readers
from callsetmerger.recordcluster import ClusterGraph
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcfs", help="Comma-separated list of VCFs to merge. Must be sorted/indexed", type=str,
                        required=True)
    parser.add_argument("--debug", help="Print helpful debug info like allele graphs", action="store_true")
    parser.add_argument("--verbose", help="Print helpful progress messages", action="store_true")
    args = parser.parse_args()

    readers = Readers(args.vcfs.split(","))

    # Check samples same in each VCF
    # TODO

    # Walk through sorted readers
    while not readers.done:
        if args.verbose:
            sys.stderr.write("Processing records in range %s:%s-%s\n"%(readers.cur_range_chrom, \
                                                                       readers.cur_range_start_pos, \
                                                                       readers.cur_range_end_pos))
        # Get mergeable calls
        rc_list = readers.getMergableCalls()
        if rc_list is None:
            # Move on
            readers.goToNext()
            continue

        # Try to merge calls in each record cluster
        for rc in rc_list:
            if args.verbose:
                sys.stderr.write("  Found record cluster with %s records\n"%len(rc.record_objs))
                for rec in rc.record_objs:
                    sys.stderr.write("    " + str(rec).strip()+"\n")
                
            allele_list = rc.GetAlleleList()
            cg = ClusterGraph(allele_list)

            ######## Print debug info #########
            if args.debug:
                nx.layout
                pos = nx.spring_layout(cg.graph, k=2 / np.sqrt(len(cg.graph.nodes)))
                nx.draw(cg.graph, pos, node_color=cg.colors)
                # for p in pos:  # push text to right
                #     pos[p][0] += 0.2
                nx.draw_networkx_labels(cg.graph, pos, labels=cg.labels)
                # nx.draw_networkx(cg.graph, node_color=cg.colors, labels=cg.labels, seed=100)

                plt.show()
            ###################################

        # Move on
        readers.goToNext()

        if args.debug:
            input('Press any key to move on to next record.\n')


if __name__ == "__main__":
    main()