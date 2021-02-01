#!/bin/bash

# Output a file with:
# chrom, pos, motif, pathogenic_threshold, freq_EUR, freq_AFR, freq_EAS, freq_AMR, freq_SAS

echo "chrom,pos,motif,thresh,freq_EUR,freq_AFR,freq_EAS,freq_AMR,freq_SAS" | sed 's/,/\t/g'

####### HTT ########
# Taken from dbgap
chrom="chr4"
pos=3074877
motif=CAG
thresh=40
freqEUR=$(cat /storage/s1saini/manuscript_strsnp/fig3/tredparse_calls/htt_control_freqs.txt | sort | uniq -c | awk '{print $2":"$1/222 "," }' | tr -d '\n' | sed 's/,$//')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freq_SAS | sed 's/;/\t/g'

exit 0
