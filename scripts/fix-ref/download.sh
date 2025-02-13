#!/bin/bash

#download and index hg38 ref panel
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta

#download bref3.jar file
wget  https://faculty.washington.edu/browning/beagle/bref3.27May24.118.jar 