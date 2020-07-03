
library(bsseq)
library(HDF5Array)
library(DelayedMatrixStats)
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs_wSDS.rda")
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_non-brain_samples.mCG.small_smooth")

BS <- keepSeqlevels(
  x = BS,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")
BS=sort(sortSeqlevels(BS))

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/QC_failures.rda")
bsseq_filtered <- BS[, !colnames(BS) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]
BS <- bsseq_filtered
rm(bsseq_filtered)

# 15 lung and 18 thyroid samples remain that passed QC

M_Lung=as.matrix(getMeth(BS[,BS$Tissue == "Lung"],regions=list_of_filtered_autosomal_vmrs_wSDS[[9]],type="smooth",what="perRegion"))
M_Thyroid=as.matrix(getMeth(BS[,BS$Tissue == "Thyroid"],regions=list_of_filtered_autosomal_vmrs_wSDS[[10]],type="smooth",what="perRegion"))
M_nonbrain_list=list(M_Thyroid,M_Lung)
names(M_nonbrain_list)=c("Thyroid","Lung")
save(M_nonbrain_list,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/Lung_Thyroid_M_list.rda")

####################### VMR effect size analysis
library(GenomicRanges)
library(matrixStats)
### VMR coordinates
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs_wSDS.rda")
## HC1 is HC(DG)
## HC2 is HC
## HC is both regions

### data points
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/VMR_all_brain_list_w_Meth.rda")
## This is missing lung and thyroid
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/Lung_Thyroid_M_list.rda")

iqrs <- lapply(M_list, rowIQRs)
sds <- lapply(M_list, rowSds)
maxs <- lapply(M_list, function(xx) { rr <- rowRanges(xx); rr[,2] - rr[,1] })

pdf(file = "tmp.pdf")
boxplot(iqrs, main = "IQRs")
boxplot(sds, main = "SDs")
boxplot(maxs, main = "Max - min")
# plot(sds[[1]], iqrs[[1]])
dev.off()


names(maxs) = c("HC (DG)", "HC", "Amygdala", "Hypo.", "NAcc", "Putamen", "BA9", "BA24", "Caudate")

pdf(file = "VMR_effect_size.pdf", height = 2, width = 2, pointsize = 6)
par(bty = "n")
boxplot(maxs, ylab = "Max - min (methylation)", outline = FALSE, xaxt = "n", yaxt = "n")
axis(side = 2, at = c(0.2, 0.35, 0.5))
axis(side = 1, las = 2, at = 1:9, labels = names(maxs))
dev.off()

