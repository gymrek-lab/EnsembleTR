#!/bin/bash

STATLIST=$1
PUBFILE=$2
OUTDIR=$3

mkdir -p $OUTDIR/tmp

# List of regions from compiled afreqs
cat $PUBFILE | awk '{print $1 "\t" $2-10 "\t" $2+10}' > ${OUTDIR}/tmp/freq_regions.bed

# Make smaller VCFs with target loci
for statfile in $(echo $STATLIST | sed 's/,/ /g')
do 
    statname=$(basename $statfile .statstr.tab)
    cat $statfile | grep -v chrom | \
	intersectBed -header -a $statfile -b ${OUTDIR}/tmp/freq_regions.bed -sorted | \
	awk -v "pop=$statname" '{print pop "\t" $0}' > ${OUTDIR}/tmp/$statname.freqs.bed
done

cat ${OUTDIR}/tmp/*.freqs.bed | grep -v "chrom" > ${OUTDIR}/tmp/allfreqs.txt

./combine_freqs.py ${OUTDIR}/tmp/allfreqs.txt 

exit 0
