from cyvcf2 import VCF
import sys
import re


# Write base pair differences per position in a file
vcf = VCF(sys.argv[1], samples=sys.argv[2])
file = open(sys.argv[3],'a')
for v in vcf:
    pos = str(v.POS)
    refLen = len(v.REF)
    gt_bases = v.gt_bases[0]
    gt_bases = re.split('/|\|',gt_bases)
    if ('.' in gt_bases):
        gt_bases_len = ['NA']
    else:
        gt_bases_len = [str(len(i) - refLen) for i in gt_bases]
    final_list = [pos]+gt_bases_len+["\n"]
    file.write("\t".join(final_list))
file.close()
