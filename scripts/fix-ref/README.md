# Clean up the EnsembleTR SNP-TR reference panel

TODO:
* Remove duplicate records and those with incorrect reference allele
* Add required header line
* Add informative variant IDs to STRs
* Run script to add beagle info
* Release bref versions alongside the VCFs

## Step 1: Fix reference VCF files

```
chrom=11
./fix_ensembletr_snpstr_reference.py \
	--vcf chr${chrom}_final_SNP_merged_additional_TRs.vcf.gz \
	--ref Homo_sapiens_assembly38.fasta \
	--out ensembletr_refpanel_v3_chr${chrom}
```

## Step 2: Add beagle info

```
TODO
```

## Step 3: Convert to bref

```
TODO
```