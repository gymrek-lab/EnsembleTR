#!/usr/bin/env python

"""
Usage:
./run_str_mend.py <vcffile> <famfile.ped>

Output: text file with one line per locus, with coluns:
chrom, start, repunit, total trios, num consistent, num informative and consistent, num inconsistnt,
mend inheritance rate, list of bad trios

where informative means not all samples in the trio are homozygous for the same allele

Example:
./run_str_mend.py /gymreklab-tscc/mousavi/results/1000genomes/outs/EUR/CEU/merged/CEU_filtered.vcf.gz 1000G.ped 
"""

import sys
import cyvcf2

try:
    VCFFILE = sys.argv[1]
    FAMFILE = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Get list of samples in the VCF
vcf = cyvcf2.VCF(VCFFILE)
vcfsamples = vcf.samples

# Infer trios from ped file
trio_list = [] # father, mother, child
trio_list_indices = [] # same as above but indices into the sample list
with open(FAMFILE, "r") as f:
    for line in f:
        items = line.strip().split()
        sample = items[1]
        father = items[2]
        mother = items[3]
        if father == "0" or mother == "0": continue
        if sample not in vcfsamples or mother not in vcfsamples or father not in vcfsamples: continue
        trio_list.append([father, mother, sample])
        trio_list_indices.append([vcfsamples.index(father), \
                                  vcfsamples.index(mother), \
                                  vcfsamples.index(sample)])

sys.stderr.write("### Found %s trios in the fam file overlapping the VCF\n"%len(trio_list))
if len(trio_list) == 0: sys.exit(0)

def IsHomRef(gt):
    return gt[0]==0 and gt[1]==0

def IsCalled(gt):
    return gt[0]>=0 and gt[1]>=0

def PrintGT(gt):
    return ",".join([str(gt[0]), str(gt[1])])

def GetTrioMetrics(gt_father, gt_mother, gt_child):
    metrics = {"informative": None, "consistent": None, "called": None, "GTs": None}
    # Check if all members called
    if not (IsCalled(gt_father) and IsCalled(gt_mother) and IsCalled(gt_child)):
        metrics["called"] = False
        return metrics
    else:
        metrics["called"] = True
    # Check if all homref
    if IsHomRef(gt_father) and IsHomRef(gt_mother) and IsHomRef(gt_child):
        metrics["informative"] = False
        metrics["consistent"] = True
        return metrics
    else:
        metrics["informative"] = True
    # Locus is called and not boring. Check for Mend. Inheritance
    if ((gt_child[0] in gt_father[0:2]) and (gt_child[1] in gt_mother[0:2])) or \
       ((gt_child[0] in gt_mother[0:2]) and (gt_child[1] in gt_father[0:2])):
        metrics["consistent"] = True
    else:
#        sys.stderr.write("Detected Mendelian error: father=%s mother=%s child=%s\n"%(PrintGT(gt_father), PrintGT(gt_mother), PrintGT(gt_child)))
        metrics["consistent"] = False
    metrics["GTs"]=";".join([PrintGT(gt_father),PrintGT(gt_mother),PrintGT(gt_child)])
    return metrics

# Go through each variant
for v in vcf:
    if v.FILTER is not None: continue
    # reset vals
    num_trios_consistent = 0
    num_trios_inconsistent = 0
    num_trios_inf_consistent = 0
    bad_trios = []
    total_trios = 0
    # Go through each trio to test    
    for i in range(len(trio_list_indices)):
        gt_father, gt_mother, gt_child = \
                [v.genotypes[trio_list_indices[i][j]] for j in range(3)]
        trio_metrics = GetTrioMetrics(gt_father, gt_mother, gt_child)
        if trio_metrics["called"]:
            total_trios += 1
            if trio_metrics["consistent"]:
                num_trios_consistent += 1
                if trio_metrics["informative"]:
                    num_trios_inf_consistent += 1
            else:
                num_trios_inconsistent += 1
                bad_trios.append(";".join(trio_list[i])+"/"+trio_metrics["GTs"])
    if total_trios == 0:
        rate= "NA"
    else: rate = num_trios_consistent/float(total_trios)
    outitems = [v.CHROM, v.POS, v.INFO["RU"], total_trios, num_trios_consistent, num_trios_inf_consistent, num_trios_inconsistent, rate, "|".join(bad_trios)] # Add trio metrics
    sys.stdout.write("\t".join([str(item) for item in outitems])+"\n")

