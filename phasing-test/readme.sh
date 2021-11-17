#!/bin/bash

# Testing out phasing STR genotypes onto SNPs

# Get SNPs
tabix --print-header ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz 12:63117929-63176031 | bgzip -c > 1000G_hg38_chr12subset_snps.vcf.gz
bcftools query -l 1000G_hg38_chr12subset_snps.vcf.gz > snp_samples.txt

# Separate run for each population
./get_phased_rs1rs3.sh AFR YRI


#VCF=/gymreklab-tscc/mousavi/results/1000genomes/hipstr_outs/EUR/CEU/merged/CEU_filtered.vcf.gz
#bcftools query -l $VCF > ceu_samples.txt
#cat ceu_samples.txt | grep -w -f snp_samples.txt > use_samples.txt

# Get STRs - unphased
#PREFIX=RS1; REGION=12:63153304-63153304
#PREFIX=RS3; REGION=12:63156351-63156351
#tabix --print-header $VCF chr$REGION | sed 's/^chr//' | sed '/^##/! s/|/\//g' | bgzip -c > 1000G_hg38_chr12subset_${PREFIX}.vcf.gz
#tabix -p vcf 1000G_hg38_chr12subset_${PREFIX}.vcf.gz

# Subset samples
#bcftools view -S use_samples.txt 1000G_hg38_chr12subset_snps.vcf.gz | bgzip -c > 1000G_hg38_chr12subset_snps_ceu.vcf.gz
#bcftools annotate -x INFO/AC,AN 1000G_hg38_chr12subset_${PREFIX}.vcf.gz | bcftools view -S use_samples.txt | bgzip -c > 1000G_hg38_chr12subset_${PREFIX}_ceu.vcf.gz
#bcftools index -f 1000G_hg38_chr12subset_snps_ceu.vcf.gz
#bcftools index -f 1000G_hg38_chr12subset_${PREFIX}_ceu.vcf.gz
#tabix -p vcf -f 1000G_hg38_chr12subset_snps_ceu.vcf.gz

# Concatenate STRs+SNPs
#bcftools concat -a 1000G_hg38_chr12subset_snps_ceu.vcf.gz 1000G_hg38_chr12subset_${PREFIX}_ceu.vcf.gz -O z -o ${PREFIX}.snpstr_ceu.vcf.gz
#bcftools index -f ${PREFIX}.snpstr_ceu.vcf.gz

# Sort concatenated file
#bcftools sort ${PREFIX}.snpstr_ceu.vcf.gz -O z -o ${PREFIX}.snpstr_ceu.sorted.vcf.gz
#bcftools index -f ${PREFIX}.snpstr_ceu.sorted.vcf.gz

# Fix the missing marker
#bcftools view ${PREFIX}.snpstr_ceu.sorted.vcf.gz | \
#    sed 's/\.\:\./\.\/\.\:\./g' | \
#    sed 's/\.\/\.\/\./\.\/\./g' | bgzip -c > ${PREFIX}.snpstr_ceu.sorted.fixmissingvcf.gz
#bcftools index -t ${PREFIX}.snpstr_ceu.sorted.fixmissingvcf.gz

# Phase with Beagle
#java -Xmx5g -jar beagle.r1399.jar \
#    gt=${PREFIX}.snpstr_ceu.sorted.fixmissingvcf.gz \
#    out=${PREFIX}.snpstr_ceu.phased \
#    usephase=true \
#    nthreads=4
#tabix -p vcf -f ${PREFIX}.snpstr_ceu.phased.vcf.gz

# Make sure phase is consistent with original SNP file
#./fix_switch_error.py \
#    --phased-vcf ${PREFIX}.snpstr_ceu.phased.vcf.gz \
#    --ref-vcf 1000G_hg38_chr12subset_snps_ceu.vcf.gz \
#    --switch-threshold 0.5 --min-maf 0.1 --check-snps 100 \
#    --new-vcf ${PREFIX}.snpstr_ceu.phased.fixed.vcf.gz

# Extract and get counts
#bcftools concat RS1.snpstr_ceu.phased.fixed.vcf.gz RS3.snpstr_ceu.phased.fixed.vcf.gz | grep "^#\|Human" | bgzip -c > rs1rs3.vcf.gz
#bcftools index rs1rs3.vcf.gz
#bcftools query -f "[%TGT\t]\n" rs1rs3.vcf.gz | sed 's/|/\t/g' | datamash transpose | sort -k 1,1 -k2,2 | awk '{print length($1) "\t" length($2)}' | sort -k1,1 -k2,2 | uniq -c | awk '($1>=10)'
