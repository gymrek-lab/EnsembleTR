#!/usr/bin/env python3

"""
Tool to merge STR calls across multiple tools
Work in progress

# Usage
EnsembleTR --out test.vcf --ref hg38.fa --vcfs ensembletr/ExampleData/advntr-chr20.vcf.gz,ensembletr/ExampleData/eh-chr20.vcf.gz,ensembletr/ExampleData/gangstr-chr20.vcf.gz,ensembletr/ExampleData/hipstr-chr20.vcf.gz
"""

import argparse
import networkx as nx
import numpy as np
import os
from pyfaidx import Fasta
import trtools.utils.utils as utils
import sys

from . import vcfio as vcfio
from . import recordcluster as recordcluster
from ensembletr import __version__

def main(args):
    if not os.path.exists(args.ref):
        common.WARNING("Error: %s does not exist"%args.ref)
        return 1
    for vcffile in args.vcfs.split(","):
        if not os.path.exists(vcffile):
            common.WARNING("Error: %s does not exist"%vcffile)
            return 1
    if not args.out.endswith("vcf"):
        common.WARNING("Error: --out must end with '.vcf'")
        return 1

    ref_genome = Fasta(args.ref)
    readers = vcfio.Readers(args.vcfs.split(","), ref_genome)
    writer = vcfio.Writer(args.out, readers.samples, " ".join(sys.argv))

    recnum = 0
    while not readers.done:
        rc_list = readers.getMergableCalls().RecordClusters
        for rc in rc_list:
            num_vcfs = len([i for i in rc.vcf_types if i == True])
            if not (num_vcfs == 1 and args.exclude_single):
                recresolver = recordcluster.RecordResolver(rc)
                if recresolver.Resolve(): 
                    writer.WriteRecord(recresolver)
            recnum += 1
        readers.goToNext()
        if args.end_after != -1 and recnum >= args.end_after:
            break
    writer.Close()

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcfs", help="Comma-separated list of VCFs to merge. Must be sorted/indexed", type=str,
                        required=True)
    inout_group.add_argument("--out", help="Output merged VCF file", type=str, required= True)
    inout_group.add_argument("--ref", help="Reference genome .fa file", type=str, required=True)
    filter_group = parser.add_argument_group("Filtering")
    debug_group = parser.add_argument_group("Debug")
    debug_group.add_argument("--end-after", help="Only process the first N records", type=int, default=-1)
    debug_group.add_argument("--exclude-single", help="Exclude TRs called by only one genotyper", default=False, action='store_true')
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()

