while getopts s:c:b:p: option
do
        case "${option}"
        in
                s) sss=${OPTARG};;
                c) chr=$OPTARG;;
                b) beagle=$OPTARG;;
                p) pop=$OPTARG;;
        esac
done

str=/projects/ps-gymreklab/helia/ensembl/ensemble_out/chrs/merged_chr"$chr"_sorted_ver2.vcf.gz
procStr=/projects/ps-gymreklab/helia/ensembl/biallelic/snpstrs/snpstrs_"$chr".vcf.gz
snp=/projects/ps-gymreklab/helia/ensembl/biallelic/snpstrs/snps_"$chr".vcf.gz

SAMPLEID=$sss
echo $SAMPLEID
echo $chr

#extract snp IDs
bcftools query -f '%ID\n' $snp > ID.${SAMPLEID}.${chr}.txt

#remove our target sample from snpstrs
bcftools view $procStr --samples ^$SAMPLEID --no-update --output-type z --output-file ref.${SAMPLEID}.${chr}.vcf.gz --force-samples

#keep only our target sample from snps
bcftools view $snp --samples $SAMPLEID --no-update --output-type z --output-file exclude.${SAMPLEID}.${chr}.vcf.gz --force-samples
bcftools index -f ref.${SAMPLEID}.${chr}.vcf.gz

#impute strs for target sample
java -Xmx4g -jar $beagle gt=exclude.${SAMPLEID}.${chr}.vcf.gz ref=ref.${SAMPLEID}.${chr}.vcf.gz out=imputed.${SAMPLEID}.${chr}
bcftools index imputed.${SAMPLEID}.${chr}.vcf.gz

#remove snps from imputed file
bcftools view imputed.${SAMPLEID}.${chr}.vcf.gz --exclude ID=@ID.${SAMPLEID}.${chr}.txt -O z -o imputed_${pop}/${chr}/imputed.str.${SAMPLEID}.${chr}.vcf.gz

#write gts for both imputed and original strs
python write_bases.py imputed_${pop}/${chr}/imputed.str.${SAMPLEID}.${chr}.vcf.gz $SAMPLEID $SAMPLEID.${chr}.imputeResult.txt
python write_bases.py $str $SAMPLEID $SAMPLEID.${chr}.groundTruth.txt

#sort and join 
sort -f ${SAMPLEID}.${chr}.groundTruth.txt > ${SAMPLEID}.${chr}.ground.sorted.txt
sort -f ${SAMPLEID}.${chr}.imputeResult.txt > ${SAMPLEID}.${chr}.imputed.sorted.txt
join -1 1 -2 1 ${SAMPLEID}.${chr}.ground.sorted.txt ${SAMPLEID}.${chr}.imputed.sorted.txt | awk 'NF==5{print}' > diff_${pop}/${chr}/${SAMPLEID}.${chr}.diff.txt
rm ${SAMPLEID}.${chr}.groundTruth.txt ${SAMPLEID}.${chr}.imputeResult.txt ${SAMPLEID}.${chr}.ground.sorted.txt ${SAMPLEID}.${chr}.imputed.sorted.txt

echo "done"

rm exclude.${SAMPLEID}.${chr}.vcf.gz ref.${SAMPLEID}.${chr}.vcf.gz ref.${SAMPLEID}.${chr}.vcf.gz.csi imputed.${SAMPLEID}.${chr}.vcf.gz imputed.${SAMPLEID}.${chr}.vcf.gz.csi
rm ID.${SAMPLEID}.${chr}.txt imputed.${SAMPLEID}.${chr}.log
