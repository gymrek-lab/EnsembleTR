#!/bin/bash

ref=$1
bref="bref3.27May24.118.jar"



#bgzip and index ref
echo "bgzip and index ref"
bgzip ${ref}.vcf
tabix -p vcf ${ref}.vcf.gz

#dowaload bref3.jar
# Check if the file exists locally
if [ ! -f "$bref" ]; then
    echo "File does not exist. Downloading..."
    wget https://faculty.washington.edu/browning/beagle/bref3.27May24.118.jar 
else
    echo "File already exists. Continuing..."
fi

#convert vcf to bref3
echo "converting vcf to bref3 format"
zcat ${ref}.vcf.gz | java -jar $bref > ${ref}.bref3
echo "Done converting"

 