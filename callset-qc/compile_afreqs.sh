#!/bin/bash

STATLIST=$1
PUBFILE=$2
OUTDIR=$3

mkdir -p $OUTDIR/tmp

# List of regions from compiled afreqs
cat $PUBFILE | grep -v chrom | awk '{print $1 "\t" $2-10 "\t" $2+10}' > ${OUTDIR}/tmp/freq_regions.bed

# Make smaller VCFs with target loci
for statfile in $(echo $STATLIST | sed 's/,/ /g')
do 
    statname=$(basename $statfile .statstr.tab)
    cat $statfile | grep -v chrom | \
	intersectBed -header -a $statfile -b ${OUTDIR}/tmp/freq_regions.bed | \
	awk -v "pop=$statname" '{print pop "\t" $0}' > ${OUTDIR}/tmp/$statname.freqs.bed
done

cat ${OUTDIR}/tmp/*.freqs.bed | grep -v "chrom" > ${OUTDIR}/tmp/allfreqs.txt

# Get one header
HEADER=$(cat ${OUTDIR}/AFR/ACB/ACB.statstr.tab | head -n 1 | sed 's/\t/,/g')

./combine_freqs.py ${OUTDIR}/tmp/allfreqs.txt $HEADER

exit 0
