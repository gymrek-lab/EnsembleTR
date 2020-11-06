#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
python callsetmerger.py --vcfs hipstr.chr21.sorted.vcf.gz,advntr.chr21.sorted.vcf.gz,gangstr.chr21.sorted.vcf.gz
"""

import argparse
from callsetmerger.callsetmerger import Readers


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
        res = readers.getMergableCalls()
        if res is not None:
            allele_list, sample_calls = res
        else:
            # Move on
            readers.goToNext()
            continue

        # Merge calls
        # mo = MergeObject(allele_list, sample_calls)

        # TODO will need to update is_min_pos_list if we merged something

        # Move on
        readers.goToNext()
        input('Press any key to move on to next record.\n')


if __name__ == "__main__":
    main()
