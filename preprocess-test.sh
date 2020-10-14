#!/bin/bash

# Merge advntr calls
VCFS=$(ls /gymreklab-tscc/mousavi/analysis/callset_merge_project/data/advntr/*.vcf.gz  | awk '{print $0 ","}' | tr -d '\n' | sed 's/,$//')
mergeSTR --vcfs $VCFS --out advntr.chr21
cat advntr.chr21.vcf | vcf-sort | bgzip -c > advntr.chr21.sorted.vcf.gz
tabix -p vcf advntr.chr21.sorted.vcf.gz

# Merge gangstr calls
VCFS=$(ls /gymreklab-tscc/mousavi/analysis/callset_merge_project/data/gangstr/*.vcf.gz  | awk '{print $0 ","}' | tr -d '\n' | sed 's/,$//')
mergeSTR --vcfs $VCFS --out gangstr.chr21
cat gangstr.chr21.vcf | vcf-sort | bgzip -c > gangstr.chr21.sorted.vcf.gz
tabix -p vcf gangstr.chr21.sorted.vcf.gz

