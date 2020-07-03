##### Justifying removal of outliers from analysis #####
library(bsseq)
library(tidyverse)
library(circlize)
library(EnrichedHeatmap)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/DMR_OR_Enrich_results.RData")

 setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis")
BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
BS <- keepSeqlevels(
  x = BS,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/QC_failures.rda")
bsseq_filtered <- BS[, !colnames(BS) %in%
                                   genotype_failures]

bsseq <- bsseq_filtered

bsseq$group <- factor(
  dplyr::case_when(
    bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~ "Cortical",
    bsseq_filtered$Tissue ==
      "Brain - Anterior cingulate cortex (BA24)" ~ "Cortical",
    bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "Hippocampus",
    bsseq_filtered$Tissue == "Brain - Amygdala" ~ "Amygdala",
    bsseq_filtered$Tissue == "Brain - Hypothalamus" ~ "Hypothalamus",
    bsseq_filtered$Tissue ==
      "Brain - Nucleus accumbens (basal ganglia)" ~ "Basal_ganglia",
    bsseq_filtered$Tissue ==
      "Brain - Putamen (basal ganglia)" ~ "Basal_ganglia",
    bsseq_filtered$Tissue ==
      "Brain - Caudate (basal ganglia)" ~ "Basal_ganglia"),
  levels = c(
    "Cortical", "Hippocampus", "Amygdala", "Hypothalamus", "Basal_ganglia")
)

mismatch=c("GTEX-13OVH-0011-R8b-SM-AFUN2","GTEX-13NYS-0011-R6b-SM-AFUKH","GTEX-139TS-0011-R6a-SM-AFUJM","GTEX-1399T-0011-R6b-SM-AFUJG","GTEX-139TS-0011-R1a-SM-AFUJL")
bsseq2=bsseq[,!colnames(bsseq) %in% c("GTEX-13OVH-0011-R8b-SM-AFUN2","GTEX-13NYS-0011-R6b-SM-AFUKH","GTEX-139TS-0011-R6a-SM-AFUJM","GTEX-1399T-0011-R6b-SM-AFUJG","GTEX-139TS-0011-R1a-SM-AFUJL")]
bsseq2=bsseq2[,!colnames(bsseq2) %in% colnames(bsseq2)[which(pData(bsseq2)$Tissue %in% c("Brain - Frontal Cortex (BA9)","Brain - Caudate (basal ganglia)","Brain - Putamen (basal ganglia)","Brain - Anterior cingulate cortex (BA24)"))]]

M=as.matrix(getMeth(bsseq2,type="smooth"))

Amy_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq2)[which(pData(bsseq2)$Tissue=="Brain - Amygdala")])])
NAcc_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq2)[which(pData(bsseq2)$Tissue=="Brain - Nucleus accumbens (basal ganglia)")])])
HC_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq2)[which(pData(bsseq2)$Tissue=="Brain - Hippocampus")])])
Hyp_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq2)[which(pData(bsseq2)$Tissue=="Brain - Hypothalamus")])])
BA9_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq)[which(pData(bsseq)$Tissue=="Brain - Frontal Cortex (BA9)")])])
BA24_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq)[which(pData(bsseq)$Tissue=="Brain - Anterior cingulate cortex (BA24)")])])
CAU_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq)[which(pData(bsseq)$Tissue=="Brain - Caudate (basal ganglia)")])])
PUT_M=rowMeans2(M[,which(colnames(M) %in% sampleNames(bsseq)[which(pData(bsseq)$Tissue=="Brain - Putamen (basal ganglia)")])])

new_M=cbind(Amy_M,NAcc_M,HC_M,Hyp_M,BA9_M,BA24_M,CAU_M,PUT_M)

#mismatch_M=as.matrix(getMeth(bsseq[,mismatch],type="smooth"))
all_M=cbind(new_M,mismatch_M)
rm(mismatch_M,Amy_M,NAcc_M,Hyp_M,HC_M,new_M,BA9_M,BA24_M,PUT_M,CAU_M)
gr=granges(bsseq2)
mcols(gr)=all_M
dmrs=GRanges(list_of_all_DMRs[["general_CG-DMRs"]])
gr=subsetByOverlaps(gr,dmrs)
dmrs_list=split(dmrs,as.factor(dmrs))
a=findOverlaps(gr,dmrs_list)
gr$group=subjectHits(a)
test=as_tibble(mcols(gr))

res <- test %>%
group_by(group) %>%
summarize_all(funs(mean))

res$group=NULL
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

Sample=c("Amy_M","NAcc_M","HC_M","Hyp_M","BA9_M","BA24_M","CAU_M","PUT_M","GTEX.13OVH.0011.R8b.SM.AFUN2","GTEX.13NYS.0011.R6b.SM.AFUKH","GTEX.139TS.0011.R6a.SM.AFUJM","GTEX.1399T.0011.R6b.SM.AFUJG","GTEX.139TS.0011.R1a.SM.AFUJL")

Tissue=c("Amygdala","Nucleus accumbens","Hippocampus","Hypothalamus","BA9","BA24","Caudate","Putamen","Hypothalamus","Nucleus accumbens","Nucleus accumbens","Nucleus accumbens","Hippocampus")
class=c("correct","correct","correct","correct","correct","correct","correct","correct","wrong","wrong","wrong","wrong","wrong")
sample_info=cbind(Sample,Tissue,class)
ha = HeatmapAnnotation(type = Tissue, class=class, 
    col = list(type = c("Amygdala"="#e7298a","Hippocampus"="#1b9e77","Hypothalamus"="#66a61e","Nucleus accumbens"="#d95f02","Caudate" = "#666666","Putamen"="#7570b3","BA9"="#e6ab02","BA24"="#a6761d"),
  class=c("correct"="blue","wrong"="red")))


set.seed(12345)
random_20k=res[sample(1:nrow(res),20000,replace=FALSE),]

pdf("Outlier_test_heatmap_dendrogram2.pdf")
Heatmap(random_20k, name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = FALSE, show_column_names = FALSE)

dev.off()

