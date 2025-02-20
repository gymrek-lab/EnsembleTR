# Clean up the EnsembleTR SNP-TR reference panel

## Step 1: Fix reference VCF files

This step:
* Removes records with incorrect reference allele
* Adds required header line for TRTools
* Adds informative variant IDs to STRs
* Adds VT field
* Remove loci with too many or too few alleles
* Remove alleles with AF=0
* Remove loci that have the same chr:pos:ref:alt after the above steps
* Strip DS/GP fields

```
# Test on chr11
chrom=11
./fix_ensembletr_snpstr_reference.py \
	--vcf chr${chrom}_final_SNP_merged_additional_TRs.vcf.gz \
	--ref Homo_sapiens_assembly38.fasta \
	--max-alleles 100 --min-alleles 2 | bgzip -c > ensembletr_refpanel_v3_chr${chrom}.vcf.gz
tabix -p vcf ensembletr_refpanel_v3_chr${chrom}.vcf.gz
```

## Step 2: Convert to bref

```
wget https://faculty.washington.edu/browning/beagle/bref3.27May24.118.jar
chrom=11
zcat ensembletr_refpanel_v3_chr${chrom}.vcf.gz | \
    java -jar bref3.27May24.118.jar > ensembletr_refpanel_v3_chr${chrom}.bref3
```