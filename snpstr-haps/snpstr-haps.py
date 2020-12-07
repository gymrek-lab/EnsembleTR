#!/usr/bin/env python3

"""
Example command:
./snpstr-haps.py --vcf /storage/s1saini/hipstr_allfilters/str_snp_parents/chr21.str.snp.feb18.vcf.gz --target 21:37064256 --window-kb 50 --min-maf 0.1
"""

import argparse
import pandas as pd
import sys
import vcf

def PruneSNPs(hapdata):
    """
    Return new data frame but remove pruned SNPs

    Idea: build graph where nodes=SNPs, edges connect SNPs with pairwise LD r2 > some threshold
         Then, only keep one SNP from each connected component

    Also look into: https://en.wikipedia.org/wiki/Four-gamete_test
         Four-gamete test detect recombination between pair of SNPs
         We might want to build up haplotypes starting at the STR and moving outwards
         Stop extending haplotypes if too many recombination events have occurred

    """
    return hapdata # TODO implement
    

def BuildTree(hapdata, outprefix):
    """
    Build phylogenetic tree of haplotypes. 
    Keep track of STR alleles/haplotypes at each leaf
    Output the tree to a file

    Notes:
      * See https://biopython.org/wiki/Phylo for working with trees in python
      * Can maybe use RaxML to build trees? Or Phylo might have its own tree building methods
      * Output tree in standard format like newick/phylip
      * See https://www.nature.com/articles/nature16549 Figure 2 for how we could visualize SNP/STR haplotypes
    """
    pass # TODO implement

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Path to SNP-STR VCF. Must be phased, gzipped, sorted, and indexed", type=str, required=True)
    parser.add_argument("--target", help="Target STR. chrom:start", type=str, required=True)
    parser.add_argument("--window-cm", help="Target STR window in cM", default=0.0001, type=float)
    parser.add_argument("--window-kb", help="Target STR window in kb", default=50, type=float)
    parser.add_argument("--gmap", help="Genetic map file. Needed to use --window-cM", type=str)
    parser.add_argument("--min-maf", help="Skip SNPs with MAF lower than this", default=0.0, type=float)
    parser.add_argument("--out", help="Output prefix", type=str)
    args = parser.parse_args()

    # Figure out target STR
    try:
        str_chrom, str_pos = args.target.split(":")
        str_pos = int(str_pos)
    except:
        sys.sterr.write("Incorrectly formatted target. must be chrom:pos of the target STR\n")
        sys.exit(1)

    # Figure out window
    window_start = -1
    window_end = -1
    if args.gmap is not None:
        # TODO implement using cM, rather than kb, to find the window
        sys.stderr.write("Using a genetic map for finding SNP window not yet implemented. Quitting. Use window-kb instead")
        sys.exit(1)
    else:
        window_start = str_pos-args.window_kb*1000
        window_end = str_pos+args.window_kb*1000
    if window_start == -1 or window_end == -1:
        sys.stderr.write("Couldn't set window. Quitting\n")
        sys.exit(1)

    # Set up VCF
    reader = vcf.Reader(open(args.vcf, "rb"))

    # Set up lists to read SNP haplotypes and corresponding STR alleles (2 for each sample, assuming diploid)
    snp_haps = []
    snp_labels = []
    str_alleles = []
    for sample in reader.samples:
        for hap in [1,2]:
            snp_haps.append([])
            str_alleles.append("NA")

    # Featch target region
    try:
        reader.fetch(str_chrom, start=window_start, end=window_end)
    except:
        sys.stderr.write("Error fetching %s:%s-%s from %s. Is VCF indexed and region set correctly?\n"%(str_chrom, window_start, window_end, args.vcf))
        sys.exit(1)
    
    # Get SNP haplotypes and corresponding STR alleles for each sample
    found_str = False # make sure we find the target STR
    for record in reader:
        is_str = False
        if (len(record.REF)!=1 or "STR" in record.ID) and (not record.POS==str_pos):
            continue # skip other STRs
        if record.POS == str_pos:
            found_str = True
            is_str = True
        if not is_str:
            alt_af = sum(record.aaf)
            if alt_af<args.min_maf or alt_af>(1-args.min_maf):
                #sys.stderr.write("Skipping... maf to low\n")
                continue
            else:
                snp_labels.append(record.ID)
        hap_ind = 0
        for sample in record:
            if not sample.called:
                if is_str:
                    str_alleles[hap_ind] = "NA"
                    str_alleles[hap_ind+1] = "NA"
                else:
                    snp_haps[hap_ind].append("NA")
                    snp_haps[hap_ind+1].append("NA")
                hap_ind += 2
            else:
                gts = sample.gt_alleles
                alleles = sample.gt_bases.split("|")
                assert(len(alleles)==2)
                assert(len(gts)==2)
                if is_str:
                    str_alleles[hap_ind] = alleles[0]
                    str_alleles[hap_ind+1] = alleles[1]
                else:
                    snp_haps[hap_ind].append(gts[0])
                    snp_haps[hap_ind+1].append(gts[1])
                hap_ind += 2
    if not found_str:
        sys.stderr.write("Could not find target str at %s\n"%args.target)
        sys.exit(1)

    # Load into pandas dataframe, one row per haplotype
    hapdata = pd.DataFrame({
        "str": str_alleles
        })
    for i in range(len(snp_haps[0])):
        hapdata[snp_labels[i]] = [hap[i] for hap in snp_haps]

    # Prune SNPs
    hapdata_pruned = PruneSNPs(hapdata)

    # Build tree
    BuildTree(hapdata_pruned, args.out)

    # Output haplotypes
    hapdata.sort_values("str").to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
