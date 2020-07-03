#Use 50-65 SNPs from 450K that are used to assess genotype concordance with that platform
snps=read.csv("/Volumes/Transcend/eGTEx/EPIC_QC_SNPs.csv",stringsAsFactors=F,header=T) # SuppTable5 from PMID:26673039
library(biomaRt) #

snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

snp_ids = snps$Name
snp_attributes = c("refsnp_id", "chr_name", "chrom_start","")

snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", 
                      values=snp_ids, mart=snp_mart)
snp_locations$chr_name=paste0("chr",snp_locations$chr_name)
colnames(snp_locations)=c("refsnp_id","chr","start")
snp_locations$end=snp_locations$start
write.csv(snp_locations,file="/Volumes/Transcend/eGTEx/EPIC_QC_SNP_CHR_positions.csv") 
#moved to cluster as well /dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/EPIC_QC_SNP_CHR_positions.csv

# Load in snps and convert to granges
snps=read.csv("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/EPIC_QC_SNP_CHR_positions.csv",header=T,stringsAsFactors=F)
snps$X=NULL
snps=GRanges(snps)
export(snps,con="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/EPIC_QC_SNP_CHR_positions.bed",format="bed")


