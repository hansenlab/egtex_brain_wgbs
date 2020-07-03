
library(SummarizedExperiment)
library(bsseq)
library(podkat)
library(HDF5Array)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

####### Annotate ALL F-stat DMRs #########

dmrs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/general_CG-DMRs.rds")
BS<-HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth/eGTEx.Phase2_brain_samples.mCG.small_smooth.general_DMRs")
BS <- keepSeqlevels(
  x = BS,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/QC_failures.rda")
bsseq_filtered <- BS[, !colnames(BS) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]
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

design.no.intercept <- model.matrix(~0 + bsseq$group)

colnames(design.no.intercept) <- gsub("bsseq\\$group", "",
                                      colnames(design.no.intercept))

contrasts.no.intercept <- limma::makeContrasts(
CortvsHC  = Cortical - Hippocampus,
CortvsBG  = Cortical - Basal_ganglia,
CortvsAMY = Cortical - Amygdala,
CortvsHYP = Cortical - Hypothalamus,
HCvsAMY = Hippocampus - Amygdala,
HCvsBG = Hippocampus - Basal_ganglia,
HCvsHYP = Hippocampus - Hypothalamus,
BGvsHYP = Basal_ganglia - Hypothalamus,
BGvsAMY = Basal_ganglia - Amygdala,
HYPvsAMY = Hypothalamus - Amygdala,
 levels = design.no.intercept)

fac <- bsseq$group
fstat_comparisons_res <-
  bsseq:::fstat.comparisons.pipeline(BSseq = bsseq,
                                     design = design.no.intercept,
                                     contrasts = contrasts.no.intercept,
                                     cutoff = 4.6,
                                     fac = fac,
                                     verbose = TRUE)

tstat_dmrs <- fstat_comparisons_res[[2]]


annotation.pipeline <- function(fstat_dmrs, fstat_comparisons_res, overlap = 0.5) {
  if (overlap > 1 || overlap <= 0) {
    stop("require 0 < 'overlap' <= 1")
  }
  # For each f-stat comparisons DMR data.frame, annotate F-stat DMRs by overlap
  dmr_labels <- sapply(fstat_comparisons_res$dmrs, function(x) {
    x <- makeGRangesFromDataFrame(x)
    h <- findOverlaps(fstat_dmrs, x)
    w <- width(pintersect(fstat_dmrs[queryHits(h)], x[subjectHits(h)]))
    ok <- w > (overlap * (width(fstat_dmrs[queryHits(h)])))
    h <- h[ok]
    countQueryHits(h) > 0
  })

  dmr_labels
}

# Add labels to dmrs
dmr_labels <- annotation.pipeline(fstat_dmrs = dmrs,fstat_comparisons_res = fstat_comparisons_res, overlap = 0.5)

mcols(dmrs) <- cbind(mcols(dmrs), as(dmr_labels, "DataFrame"))

# Labels of F-stat DMRs by number of comparisons with t-stat DMRs

cs <- sapply(split(dmr_labels, rowSums(dmr_labels)), function(x) {
  colSums(matrix(x, ncol = ncol(dmr_labels), dimnames = dimnames(dmr_labels)))
})

cs <- cbind(cs, matrix(rowSums(cs), dimnames = list(NULL, "Total")))

save(fstat_comparisons_res,dmr_labels,cs,dmrs,tstat_dmrs,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/Fstat_annotation_generalCG-DMRs.RData",compress=TRUE)



### Annotate the OVERALL DMRs (between all 8 tissues not just groups)
dmrs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/overall_CG-DMRs.rds")
dmrs=GRanges(dmrs) #174482 total
bsseq$group <- factor(
  dplyr::case_when(
    bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~ "BA9",
    bsseq_filtered$Tissue ==
      "Brain - Anterior cingulate cortex (BA24)" ~ "BA24",
    bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC",
    bsseq_filtered$Tissue == "Brain - Amygdala" ~ "AMY",
    bsseq_filtered$Tissue == "Brain - Hypothalamus" ~ "HYPO",
    bsseq_filtered$Tissue ==
      "Brain - Nucleus accumbens (basal ganglia)" ~ "NACC",
    bsseq_filtered$Tissue ==
      "Brain - Putamen (basal ganglia)" ~ "PUT",
    bsseq_filtered$Tissue ==
      "Brain - Caudate (basal ganglia)" ~ "CAU"),
  levels = c(
    "BA9", "BA24","HC", "AMY", "HYPO", "PUT","NACC","CAU")
)

design.no.intercept <- model.matrix(~0 + bsseq$group)

colnames(design.no.intercept) <- gsub("bsseq\\$group", "",
                                      colnames(design.no.intercept))

contrasts.no.intercept <- limma::makeContrasts(
BA9vsHC  = BA9 - HC,
BA9vsBA24  = BA9 - BA24,
BA9vsAMY = BA9 - AMY,
BA9vsHYP = BA9 - HYPO,
BA9vsNACC  = BA9 - NACC,
BA9vsPUT = BA9 - PUT,
BA9vsCAU = BA9 - CAU,
HCvsAMY = HC - AMY,
HCvsBA24 = HC- BA24,
HCvsHYP = HC - HYPO,
HCvsCAU = HC - CAU,
HCvsNACC = HC- NACC,
HCvsPUT = HC - PUT,
HYPvsBA24= HYPO - BA24,
HYPvsAMY= HYPO - AMY,
HYPvsPUT= HYPO - PUT,
HYPvsCAU= HYPO - CAU,
HYPvsNACC= HYPO - NACC,
AMYvsBA24 = AMY - BA24,
AMYvsPUT = AMY - PUT,
AMYvsCAU = AMY - CAU,
AMYvsNACC = AMY - NACC,
CAUvsBA24 = CAU - BA24,
CAUvsPUT = CAU - PUT,
CAUvsNACC = CAU - NACC,
PUTvsNACC = PUT - NACC,
PUTvsBA24 = PUT - BA24,
NACCvsBA24 = NACC - BA24,
 levels = design.no.intercept)

fac <- bsseq$group
fstat_comparisons_res <-
  bsseq:::fstat.comparisons.pipeline(BSseq = bsseq,
                                     design = design.no.intercept,
                                     contrasts = contrasts.no.intercept,
                                     cutoff = 4.6,
                                     fac = fac,
                                     verbose = TRUE)

tstat_dmrs <- fstat_comparisons_res[[2]]


# Add labels to dmrs
dmr_labels <- annotation.pipeline(fstat_dmrs = dmrs,fstat_comparisons_res = fstat_comparisons_res, overlap = 0.5)

mcols(dmrs) <- cbind(mcols(dmrs), as(dmr_labels, "DataFrame"))

cs <- sapply(split(dmr_labels, rowSums(dmr_labels)), function(x) {
  colSums(matrix(x, ncol = ncol(dmr_labels), dimnames = dimnames(dmr_labels)))
})

cs <- cbind(cs, matrix(rowSums(cs), dimnames = list(NULL, "Total")))

save(fstat_comparisons_res,dmr_labels,cs,dmrs,tstat_dmrs,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/Fstat_annotation_OVERALL_CG-DMRs.RData",compress=TRUE)


