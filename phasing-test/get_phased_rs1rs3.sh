#!/bin/bash

set -e

SUPERPOP=$1
POP=$2

VCF=/gymreklab-tscc/mousavi/results/1000genomes/hipstr_outs/${SUPERPOP}/${POP}/merged/${POP}_filtered.vcf.gz
bcftools query -l $VCF > ${POP}_samples.txt
cat ${POP}_samples.txt | grep -w -f snp_samples.txt > use_samples_${POP}.txt

for PREFIX in RS1 RS3
do
    if [ "x$PREFIX" = "xRS1" ]
    then
	REGION=12:63153304-63153304
    else
	REGION=12:63156351-63156351
    fi
    echo $SUPERPOP $POP $PREFIX $REGION

    # Get unphased STRs
    tabix --print-header $VCF chr$REGION | sed 's/^chr//' | sed '/^##/! s/|/\//g' | bgzip -c > 1000G_hg38_chr12subset_${PREFIX}_${POP}all.vcf.gz
    tabix -p vcf 1000G_hg38_chr12subset_${PREFIX}_${POP}all.vcf.gz

    # Subset samples
    bcftools view -S use_samples_${POP}.txt 1000G_hg38_chr12subset_snps.vcf.gz | bgzip -c > 1000G_hg38_chr12subset_snps_${POP}.vcf.gz
    bcftools annotate -x INFO/AC,AN 1000G_hg38_chr12subset_${PREFIX}_${POP}all.vcf.gz | bcftools view -S use_samples_${POP}.txt | \
	bgzip -c > 1000G_hg38_chr12subset_${PREFIX}_${POP}.vcf.gz
    bcftools index -f 1000G_hg38_chr12subset_snps_${POP}.vcf.gz
    bcftools index -f 1000G_hg38_chr12subset_${PREFIX}_${POP}.vcf.gz
    tabix -p vcf -f 1000G_hg38_chr12subset_snps_${POP}.vcf.gz

    # Concatenate STRs+SNPs
    bcftools concat -a 1000G_hg38_chr12subset_snps_${POP}.vcf.gz 1000G_hg38_chr12subset_${PREFIX}_${POP}.vcf.gz -O z -o ${PREFIX}.snpstr_${POP}.vcf.gz
    bcftools index -f ${PREFIX}.snpstr_${POP}.vcf.gz

    # Sort concatenated file
    bcftools sort ${PREFIX}.snpstr_${POP}.vcf.gz -O z -o ${PREFIX}.snpstr_${POP}.sorted.vcf.gz
    bcftools index -f ${PREFIX}.snpstr_${POP}.sorted.vcf.gz

    # Fix the missing marker
    bcftools view ${PREFIX}.snpstr_${POP}.sorted.vcf.gz | \
	sed 's/\.\:\./\.\/\.\:\./g' | \
	sed 's/\.\/\.\/\./\.\/\./g' | bgzip -c > ${PREFIX}.snpstr_${POP}.sorted.fixmissingvcf.gz
    bcftools index -t ${PREFIX}.snpstr_${POP}.sorted.fixmissingvcf.gz

    # Phase with Beagle
    java -Xmx5g -jar beagle.r1399.jar \
	gt=${PREFIX}.snpstr_${POP}.sorted.fixmissingvcf.gz \
	out=${PREFIX}.snpstr_${POP}.phased \
	usephase=true \
	nthreads=4
    tabix -p vcf -f ${PREFIX}.snpstr_${POP}.phased.vcf.gz

    # Make sure phase is consistent with original SNP file
    ./fix_switch_error.py \
	--phased-vcf ${PREFIX}.snpstr_${POP}.phased.vcf.gz \
	--ref-vcf 1000G_hg38_chr12subset_snps_${POP}.vcf.gz \
	--switch-threshold 0.5 --min-maf 0.1 --check-snps 100 \
	--new-vcf ${PREFIX}.snpstr_${POP}.phased.fixed.vcf.gz
done

bcftools concat RS1.snpstr_${POP}.phased.fixed.vcf.gz RS3.snpstr_${POP}.phased.fixed.vcf.gz | grep "^#\|Human" | bgzip -c > rs1rs3_${POP}.vcf.gz
bcftools index rs1rs3_${POP}.vcf.gz
bcftools query -f "[%TGT\t]\n" rs1rs3_${POP}.vcf.gz | \
    sed 's/|/\t/g' | datamash transpose | sort -k 1,1 -k2,2 | \
    uniq -c | awk '($1>=5)' | awk -v"pop=$POP" '{print pop "\t" $0}' > rs1rs3_haps_${POP}.tab



