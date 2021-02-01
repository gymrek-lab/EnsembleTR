#!/bin/bash

VCFLIST=$1
OUTDIR=$2

FAI=/storage/resources/dbase/human/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

mkdir -p ${OUTDIR}/tmp

# First, load capillary data from spreadsheet
# Get: chrom, pos, repunit, sample, gt
./load_cap.py > ${OUTDIR}/tmp/capillary.tab
cat ${OUTDIR}/tmp/capillary.tab | cut -f 5 | sort | uniq > ${OUTDIR}/tmp/capillary_samples.tab
cat ${OUTDIR}/tmp/capillary.tab | awk '{print $1 "\t" $2-10 "\t" $2+10}' | sort -k1,1 -k2,2n | uniq > ${OUTDIR}/tmp/capillary_regions.bed

# Second, make a smaller VCF file with only the target
# samples and loci needed
# only keep loci that are passing
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
	    awk '(($1~/^#/) || ($7=="PASS"))' \
	    > ${OUTDIR}/tmp/${vcfname}.vcf
	bgzip -f ${OUTDIR}/tmp/${vcfname}.vcf
	tabix -p vcf -f ${OUTDIR}/tmp/${vcfname}.vcf.gz
    fi
done

# Merge to a single VCF
mergevcfs=$(ls ${OUTDIR}/tmp/*_filtered.vcf.gz | awk '{print $0","}' | tr -d '\n' | sed 's/,$//')
echo "merging $mergevcfs"
mergeSTR --vcfs $mergevcfs --out ${OUTDIR}/tmp/merged
bgzip -f ${OUTDIR}/tmp/merged.vcf
tabix -p vcf -f ${OUTDIR}/tmp/merged.vcf.gz

# Pull out REPCN for cap calls for comparison later
bcftools query -f "[%CHROM\t%POS\t%RU\t%SAMPLE\t%REPCN\n]" ${OUTDIR}/tmp/merged.vcf.gz > ${OUTDIR}/tmp/wgscalls.tab

# Third, convert capillary data to a VCF file in
# HipSTR-like format that can be read by compareSTR
./convert_cap.py ${OUTDIR}/tmp/capillary.tab ${OUTDIR}/tmp/merged.vcf.gz > ${OUTDIR}/tmp/capillary.vcf
vcf-sort ${OUTDIR}/tmp/capillary.vcf | bgzip -c > ${OUTDIR}/tmp/capillary.sorted.vcf.gz
tabix -p vcf -f ${OUTDIR}/tmp/capillary.sorted.vcf.gz
bcftools reheader -f $FAI -o ${OUTDIR}/tmp/capillary.sorted.reheader.vcf.gz ${OUTDIR}/tmp/capillary.sorted.vcf.gz
tabix -p vcf -f  ${OUTDIR}/tmp/capillary.sorted.reheader.vcf.gz

# Fourth, run compareSTR
compareSTR \
    --vcf1 ${OUTDIR}/tmp/merged.vcf.gz \
    --vcf2 ${OUTDIR}/tmp/capillary.sorted.reheader.vcf.gz \
    --vcftype2 gangstr \
    --out ${OUTDIR}/comparestr
