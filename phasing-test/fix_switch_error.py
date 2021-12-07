#!/usr/bin/env python3
"""
Align phase between target VCF and reference VCF
Example:
./fix_switch_error.py \
  --phased-vcf /storage/s1saini/gl_template/beagle.str1.phased.reorder.vcf.gz \
  --ref-vcf /storage/s1saini/gl_template/shapeit.chr21.with.ref.reorder.vcf.gz \
  --switch-threshold 0.5 --min-maf 0.1 --check-snps 50
"""

import argparse
from cyvcf2 import VCF, Writer
import sys
import vcf

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--ref-vcf", help="Match to phase in this VCF", required=True, type=str)
    parser.add_argument("--phased-vcf", help="Target VCF", required=True, type=str)
    parser.add_argument("--switch-threshold", help="Switch error threshold", default=0.5, type=float)
    parser.add_argument("--new-vcf", help="Output VCF file. Default stdout", required=False, type=str)
    parser.add_argument("--check-snps", help="Only check this many SNPs", default=1000000, type=int)
    parser.add_argument("--min-maf", help="MAF threshold to check SNP", default=0.0, type=float)
    args = parser.parse_args()
    phased = args.phased_vcf
    reference = args.ref_vcf
    outvcf = args.new_vcf
    if outvcf is None: outvcf = "/dev/stdout"
    samples_to_switch = {} # sample -> [total_het, total_switch]
    
    target_reader = vcf.Reader(open(args.phased_vcf, "rb"))
    ref_reader = vcf.Reader(open(args.ref_vcf, "rb"))

    for sample in target_reader.samples:
        samples_to_switch[sample] = [0, 0]

    snp_counter = 0
    for target_record in target_reader:
        if min([target_record.aaf[0], 1-target_record.aaf[0]]) < args.min_maf: continue
        if snp_counter > args.check_snps: break
#        print target_record.POS, snp_counter, target_record.aaf
        snp_counter += 1
        # Fetch corresponding record in ref vcf
        records = ref_reader.fetch(target_record.CHROM, target_record.POS-1, target_record.POS)
        ref_record = None
        for r in records:
            if r.POS == target_record.POS:
                ref_record = r
                break
        if ref_record is None:
            sys.stderr.write("WARNING: COuld not find record for %s:%s\n"%(target_record.CHROM, target_record.POS))
            continue
        # Loop through samples
        for sample in target_record:
            if sample.is_het:
                if not ref_record.genotype(sample.sample).is_het:
                    sys.stderr.write("WARNING: Site %s:%s not het in both VCFs. Skipping\n"%(target_record.CHROM, target_record.POS))
                    continue
                samples_to_switch[sample.sample][0] += 1
                # Look for sample in ref
                if sample.gt_alleles[0] != ref_record.genotype(sample.sample).gt_alleles[0]:
                    samples_to_switch[sample.sample][1] += 1
                    
    # Get switch error rate per sample
    print("getting switch error per sample")
    sample_to_switch_rate = {}
    for sample in samples_to_switch:
        if samples_to_switch[sample][0] > 0:
            sample_to_switch_rate[sample] = samples_to_switch[sample][1]*1.0/samples_to_switch[sample][0]
        else: sample_to_switch_rate[sample] = 0

    print("open vcf")
    vcf1 = VCF(phased)
    samplesInvalid = [sample for sample in sample_to_switch_rate.keys() if sample_to_switch_rate[sample]>args.switch_threshold]
    sampinds = [vcf1.samples.index(sample) for sample in samplesInvalid]
    print(len(samplesInvalid))

    print("write vcf")
    # Switch phase for incorrect samples
    w = Writer(outvcf, vcf1)
    for v in vcf1:
        gtData = v.genotypes
        for sampind in sampinds:
            gtData[sampind][0],gtData[sampind][1] = gtData[sampind][1], gtData[sampind][0]
        v.genotypes = gtData
        w.write_record(v)  
    w.close()

if __name__ == "__main__":
    main()
