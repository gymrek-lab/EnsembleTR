#!/bin/bash

# Output a file with:
# chrom, pos, motif, pathogenic_threshold, freq_EUR, freq_AFR, freq_EAS, freq_AMR, freq_SAS

echo "chrom,pos,motif,thresh,freq_EUR,freq_AFR,freq_EAS,freq_AMR,freq_SAS" | sed 's/,/\t/g'

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

####### HTT ########
# Taken from dbgap
freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr4"
pos=3074877
motif=CAG
thresh=40
freqEUR=$(cat /storage/s1saini/manuscript_strsnp/fig3/tredparse_calls/htt_control_freqs.txt | sort | uniq -c | awk '{print $2":"$1/222 "," }' | tr -d '\n' | sed 's/,$//')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'


####### PPP2R2B (SCZ12) ########
# https://www.sciencedirect.com/science/article/pii/S0304394007010518
freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr5"
pos=146878728
motif="GCT"
thresh=43 # https://pubmed.ncbi.nlm.nih.gov/27864267/
freqEUR="9:0.02,10:0.59,11:0.01,13:0.11,14:0.09,15:0.15,16:0.01,18:0.02"
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'

####### DMPK ###########
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5387946/ CHINESE
freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr19"
pos=45770205
motif="CAG"
thresh=-1
freqEAS="5:0.31,11:0.14,12:0.30,13:0.14,14:0.03,15:0.05,16:0.01,24:0.01,27:0.01"
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'
exit 0
