#!/bin/bash

# Testing out phasing STR genotypes onto SNPs

# Get SNPs
tabix --print-header ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz 12:63117929-63176031 | bgzip -c > 1000G_hg38_chr12subset_snps.vcf.gz
bcftools query -l 1000G_hg38_chr12subset_snps.vcf.gz > snp_samples.txt

# Separate run for each population
./get_phased_rs1rs3.sh AFR YRI
./get_phased_rs1rs3.sh AFR LWK
./get_phased_rs1rs3.sh AFR ACB
./get_phased_rs1rs3.sh EUR CEU
./get_phased_rs1rs3.sh EUR GBR
./get_phased_rs1rs3.sh EAS CHB
./get_phased_rs1rs3.sh EAS JPT
./get_phased_rs1rs3.sh AMR PEL
./get_phased_rs1rs3.sh AMR MXL
./get_phased_rs1rs3.sh SAS GIH
./get_phased_rs1rs3.sh SAS BEB

