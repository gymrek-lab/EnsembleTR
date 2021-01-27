#!/bin/bash

VCFLIST=$1
OUTDIR=$2

mkdir -p ${OUTDIR}/tmp

# First, load capillary data from spreadsheet
# Get: chrom, pos, repunit, sample, gt
./load_cap.py > ${OUTDIR}/tmp/capillary.tab
cat ${OUTDIR}/tmp/capillary.tab | cut -f 5 | sort | uniq > ${OUTDIR}/tmp/capillary_samples.tab
cat ${OUTDIR}/tmp/capillary.tab | awk '{print $1 "\t" $2-10 "\t" $2+10}' | sort -k1,1 -k2,2n | uniq > ${OUTDIR}/tmp/capillary_regions.bed

# Second, make a smaller VCF file with only the target
# samples and loci needed
for vcf in $(echo $VCFLIST | sed 's/,/ /g')
do 
    vcfname=$(basename $vcf .vcf.gz)
    bcftools query -l $vcf > ${OUTDIR}/tmp/${vcfname}.samples
    grep -w -f ${OUTDIR}/tmp/capillary_samples.tab ${OUTDIR}/tmp/${vcfname}.samples > ${OUTDIR}/tmp/${vcfname}.samples.keep

    numsamp=$(wc -l ${OUTDIR}/tmp/${vcfname}.samples.keep | cut -f 1 -d' ')

    if [[ $numsamp -gt 0 ]]
    then
	intersectBed -header -a $vcf -b ${OUTDIR}/tmp/capillary_regions.bed | \
	    vcf-subset -c ${OUTDIR}/tmp/${vcfname}.samples.keep | \
	    bgzip -c > ${OUTDIR}/tmp/${vcfname}.vcf.gz
	tabix -p vcf -f ${OUTDIR}/tmp/${vcfname}.vcf.gz
    fi
done

# Merge to a single VCF
# TODO

# Third, convert capillary data to a VCF file in
# HipSTR-like format that can be read by compareSTR
# TODO

# Fourth, run compareSTR
# TODO


echo "working on it"
