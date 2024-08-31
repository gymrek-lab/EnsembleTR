# Clean up the EnsembleTR SNP-TR reference panel

## Step 1: Fix reference VCF files

This step:
* Removes records with incorrect reference allele
* Adds required header line for TRTools
* Adds informative variant IDs to STRs
* Remove loci with too many or too few alleles
* Remove loci with AF=0

```
# Test on chr11
chrom=11
./fix_ensembletr_snpstr_reference.py \
	--vcf chr${chrom}_final_SNP_merged_additional_TRs.vcf.gz \
	--ref Homo_sapiens_assembly38.fasta \
	--out ensembletr_refpanel_v3_chr${chrom} \
	--max-alleles 50 --min-alleles 2
```

## Step 2: Convert to bref

```
TODO
```