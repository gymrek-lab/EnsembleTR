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

####### D5S818 ###########

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr5"
pos=123775556
motif="AGAT"
thresh=-1
freqEUR=$(cat CODIS_afreq.csv |grep D5S818 | grep EUR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAFR=$(cat CODIS_afreq.csv |grep D5S818 | grep AFR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqEAS=$(cat CODIS_afreq.csv |grep D5S818 | grep EAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAMR=$(cat CODIS_afreq.csv |grep D5S818 | grep AMR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqSAS=$(cat CODIS_afreq.csv |grep D5S818 | grep SAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'

####### TPOX ###########

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr2"
pos=1489653
motif="AATG"
thresh=-1
freqEUR=$(cat CODIS_afreq.csv |grep TPOX | grep EUR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAFR=$(cat CODIS_afreq.csv |grep TPOX | grep AFR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqEAS=$(cat CODIS_afreq.csv |grep TPOX | grep EAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAMR=$(cat CODIS_afreq.csv |grep TPOX | grep AMR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqSAS=$(cat CODIS_afreq.csv |grep TPOX | grep SAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'

####### CSF1PO ###########

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr5"
pos=150076324
motif="AGAT"
thresh=-1
freqEUR=$(cat CODIS_afreq.csv |grep CSF1PO | grep EUR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAFR=$(cat CODIS_afreq.csv |grep CSF1PO | grep AFR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqEAS=$(cat CODIS_afreq.csv |grep CSF1PO | grep EAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAMR=$(cat CODIS_afreq.csv |grep CSF1PO | grep AMR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqSAS=$(cat CODIS_afreq.csv |grep CSF1PO | grep SAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'

####### D7S820 ###########

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr7"
pos=84160226
motif="GATA"
thresh=-1
freqEUR=$(cat CODIS_afreq.csv |grep D7S820 | grep EUR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAFR=$(cat CODIS_afreq.csv |grep D7S820 | grep AFR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqEAS=$(cat CODIS_afreq.csv |grep D7S820 | grep EAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAMR=$(cat CODIS_afreq.csv |grep D7S820 | grep AMR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqSAS=$(cat CODIS_afreq.csv |grep D7S820 | grep SAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'

####### D8S1179 ###########

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr8"
pos=124894873
motif="TATC"
thresh=-1
freqEUR=$(cat CODIS_afreq.csv |grep D8S1179 | grep EUR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAFR=$(cat CODIS_afreq.csv |grep D8S1179 | grep AFR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqEAS=$(cat CODIS_afreq.csv |grep D8S1179 | grep EAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAMR=$(cat CODIS_afreq.csv |grep D8S1179 | grep AMR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqSAS=$(cat CODIS_afreq.csv |grep D8S1179 | grep SAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'

####### TH01 ###########

freqEUR=""
freqAFR=""
freqEAS=""
freqAMR=""
freqSAS=""

chrom="chr11"
pos=2171088
motif="TCAT"
thresh=-1
freqEUR=$(cat CODIS_afreq.csv |grep TH01 | grep EUR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAFR=$(cat CODIS_afreq.csv |grep TH01 | grep AFR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqEAS=$(cat CODIS_afreq.csv |grep TH01 | grep EAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqAMR=$(cat CODIS_afreq.csv |grep TH01 | grep AMR | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
freqSAS=$(cat CODIS_afreq.csv |grep TH01 | grep SAS | sed 's/, /|/g' | cut -d',' -f5 | sed 's/|/,/g' | sed 's/"/ /g' | sed 's/ //g')
echo $chrom";"$pos";"$motif";"$thresh";"$freqEUR";"$freqAFR";"$freqEAS";"$freqAMR";"$freqSAS | sed 's/;/\t/g'


exit 0

