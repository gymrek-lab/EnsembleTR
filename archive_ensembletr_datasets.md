# EnsembleTR archived datasets

This page contains links to older versions of ensembleTR files

## Version III of reference SNP+TR haplotype panel for imputation of TR variants

### Dataset description
These files contain:
* [Phased SNP and indel variants](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) of 3,202 samples from the 1000 Genomes Project (1kGP).
* TRs phased/imputed from 3,202 1kGP samples based on EnsembleTR calls.

There are in total 1,070,685 TRs and 70,692,015 SNPs/indels.

All the coordinates are based on **hg38** human reference genome.

These files contain the same data as [Version II](archive_ensembletr_datasets.md), with the following updates to facilitate use in downstream imputation pipelines:

1. Remove TRs for which the REF allele does not match the expected sequence based on CHR:POS
2. For each TR, remove alelles with 0 count.
3. Remove TRs which have more than 100 alleles.
4. Remove TRs which have less than 2 alleles.
5. Remove the DS/GP fields which are large and not used by downstream steps.
6. Add unique IDs for each TR of the format EnsTR:CHROM:POS. For TRs with the same CHR:POS, add the duplicate number of the TR following format: EnsTR:CHROM:POS:Duplicate_num. Duplicated loci with identical alleles are removed.
7. Add VT field, set to VT=TR for TRs and VT=OTHER for other variant types
8. Add the bref format files which have the same information as the VCFs but can improve Beagle imputation performance.

### Availability

All file description and download links can be found [here](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_3_readme.txt). Data and links for each chromosome for the Verson IV panel are also provided below.

Chromosome 1 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr1.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr1.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr1.bref3)] SNPs/indels=5,759,060 TRs=92,372

Chromosome 2 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr2.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr2.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr2.bref3)] SNPs/indels=6,088,598 TRs=91,132

Chromosome 3 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr3.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr3.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr3.bref3)] SNPs/indels=4,983,185 TRs=75,233

Chromosome 4 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr4.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr4.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr4.bref3)] SNPs/indels=4,875,465 TRs=69,325

Chromosome 5 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr5.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr5.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr5.bref3)] SNPs/indels=4,536,819 TRs=66,492

Chromosome 6 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr6.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr6.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr6.bref3)] SNPs/indels=4,315,217 TRs=65,937

Chromosome 7 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr7.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr7.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr7.bref3)] SNPs/indels=4,137,254 TRs=59,409

Chromosome 8 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr8.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr8.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr8.bref3)] SNPs/indels=3,886,222 TRs=55,141

Chromosome 9 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr9.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr9.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr9.bref3)] SNPs/indels=3,165,513 TRs=44,188

Chromosome 10 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr10.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr10.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr10.bref3)] SNPs/indels=3,495,473 TRs=51,637

Chromosome 11 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr11.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr11.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr11.bref3)] SNPs/indels=3,423,341 TRs=49,598

Chromosome 12 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr12.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr12.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr12.bref3)] SNPs/indels=3,332,788 TRs=55,883

Chromosome 13 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr13.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr13.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr13.bref3)] SNPs/indels=2,509,179 TRs=35,720

Chromosome 14 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr14.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr14.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr14.bref3)] SNPs/indels=2,290,400 TRs=36,199

Chromosome 15 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr15.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr15.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr15.bref3)] SNPs/indels=2,109,285 TRs=32,336

Chromosome 16 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr16.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr16.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr16.bref3)] SNPs/indels=2,362,361 TRs=35,450

Chromosome 17 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr17.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr17.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr17.bref3)] SNPs/indels=2,073,624 TRs=38,377

Chromosome 18 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr18.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr18.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr18.bref3)] SNPs/indels=1,963,845 TRs=28,443

Chromosome 19 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr19.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr19.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr19.bref3)] SNPs/indels=1,670,692 TRs=33,535

Chromosome 20 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr20.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr20.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr20.bref3)] SNPs/indels=1,644,384 TRs=25,742

Chromosome 21 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr21.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr21.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr21.bref3)] SNPs/indels=1,002,753 TRs=12,893

Chromosome 22 [[VCF](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr22.vcf.gz)] [[tbi](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr22.vcf.gz.tbi)] [[bref](https://ensemble-tr.s3.us-east-2.amazonaws.com/ensembletr-refpanel-v3/ensembletr_refpanel_v3_chr22.bref3)] SNPs/indels=1,066,557 TRs=15,643

## Version II of reference SNP+TR haplotype panel for imputation of TR variants

### Dataset description

[Phased variants](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) of 3,202 samples from the 1000 Genomes Project (1kGP).

TRs imputed from 3,202 1kGP samples.

Total 70,692,015 variants + 1,091,550 TR markers.

All the coordinates are based on **hg38** human reference genome.

### Availability

Chromosome 1 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr1_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr1_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 2 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr2_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr2_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 3 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr3_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr3_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 4 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr4_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr4_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 5 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr5_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr5_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 6 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr6_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr6_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 7 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr7_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr7_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 8 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr8_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr8_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 9 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr9_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr9_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 10 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr10_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr10_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 11 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr11_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr11_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 12 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr12_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr12_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 13 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr13_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr13_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 14 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr14_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr14_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 15 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr15_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr15_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 16 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr16_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr16_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 17 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr17_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr17_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 18 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr18_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr18_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 19 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr19_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr19_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 20 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr20_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr20_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 21 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr21_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr21_final_SNP_merged_additional_TRs.vcf.gz.tbi)

Chromosome 22 [VCF file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr22_final_SNP_merged_additional_TRs.vcf.gz) and [tbi file](https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr22_final_SNP_merged_additional_TRs.vcf.gz.tbi)


## Version I

For version I of EnsembleTR calls, please use 
https://ensemble-tr.s3.us-east-2.amazonaws.com/split/ensemble_chr"$chr"_filtered.vcf.gz for VCF file and https://ensemble-tr.s3.us-east-2.amazonaws.com/split/ensemble_chr"$chr"_filtered.vcf.gz.tbi for tbi file.

For version I of phased panels, please use 
https://ensemble-tr.s3.us-east-2.amazonaws.com/phased-split/chr"$chr"_final_SNP_merged.vcf.gz for VCF file and https://ensemble-tr.s3.us-east-2.amazonaws.com/phased-split/chr"$chr"_final_SNP_merged.vcf.gz.csi for tbi file.
