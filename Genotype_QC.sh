#!/bin/bash
#$ -M lindsay.rizzardi@jhmi.edu
#$ -m e

INPUT_DIR=$1
VCF_NAME=$2

set -e
cd ${INPUT_DIR}

echo -e "Making bam list..."
#ls -1 ${INPUT_DIR}/HiSeq*/GTEX*/*subset.sorted.bam > ${INPUT_DIR}/List_Sorted_BAM.txt


echo -e "Running mpileup..."
/jhpce/shared/community/core/samtools/1.6/bin/samtools mpileup -v -t DP,AD -f /dcl01/FB2/data/personal/gtex/hg38/hg38_noALT_with_lambda.fasta --positions /dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/EPIC_QC_SNP_CHR_positions.bed -R -b ${INPUT_DIR}/List_Sorted_BAM.txt | bcftools call -m > ${VCF_NAME}.vcf

wait
echo -e "Indexing vcf..."
bgzip ${INPUT_DIR}/${VCF_NAME}.vcf
wait
bcftools index ${INPUT_DIR}/${VCF_NAME}.vcf.gz

wait


for s in $(cat ${INPUT_DIR}/List_Sorted_BAM.txt)
do
sample=`basename ${s} | cut -d "." -f1`	
REF_ID=`basename ${s} | awk -F "-" '{FS = "-"} {print "GTEX-"$2}'`
echo -e "Running gtcheck for $s against $REF_ID"
cd ${INPUT_DIR}
bcftools gtcheck -a -p ${sample}_output -S ${REF_ID} -s ${s} -R /dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/EPIC_QC_SNP_CHR_positions.bed -g /dcl01/FB2/data/personal/gtex/dbgap/dbGaP-7882/59386/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz ${INPUT_DIR}/${VCF_NAME}.vcf.gz 

tail -n 1 ${INPUT_DIR}/${sample}_output.tab | awk '{print $5}' >> ${INPUT_DIR}/matching_GTEX_IND.txt
done
echo -e "COMPLETE!"


