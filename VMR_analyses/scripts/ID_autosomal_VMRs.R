##### Identify VMRs #####
library(bsseq)
library(DelayedMatrixStats)
library(scales)
library(HDF5Array)
# NOTE: Do not use /dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain/HDF5/small_autosomes
#       without filtering low coverage sites; these low coverage regions are
#       exactly the regions that otherwise come up as 'VMRs'
gtex <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
gtex2 <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_non-brain_samples.mCG.small_smooth")

gtex=combineList(gtex,gtex2)
rm(gtex2)
options("mc.cores" = 12)

# Exclude several samples that appear to be outliers in the PCA
## SM-AFUOS: Nominally BA9 but clusters with one of the two HC clusters ( failed GT QC )
## SM-AFUJL: Nominally HC but clusters with hypothalamus
## SM-AFUKH: Nominally NAcc but clusters with amygdala
## SM-AFUJM: Nominally NAcc but (loosely) clusters with amygdala
## SM-AFUJG: Nominally NAcc but clusters with hypothalamus
## SM-AFUJX: Nominally NAcc but clusters with hypothalamus (failed GT QC)
## SM-AFUN2: Nominally hypothalamus but clusters with HC
outliers <- c("SM-AFUJL", "SM-AFUKH", "SM-AFUJM", "SM-AFUJG", "SM-AFUN2")

gtex <- gtex[, !gtex$Sample_ID %in% outliers]
QC_fail=c("GTEX-13G51-0011-R6b-SM-AFUJX","GTEX-13G51-0011-R10b-SM-AFUJU","GTEX-13N1W-0826-SM-ATE4J","GTEX-13O3Q-0011-R3a-SM-AFUKX",
"GTEX-T2IS-0011-R10a-SM-AFUOS","GTEX-13N1W-0011-R5b-SM-AFUQW","GTEX-13N1W-0011-R1b-SM-AFUQX","GTEX-13N1W-0011-R7b-SM-AFUQT",
"GTEX-13O3Q-0011-R6a-SM-AFUKS","GTEX-WZTO-0426-SM-ATE6F","GTEX-13N1W-0726-SM-ATE66","GTEX-13FHO-1026-SM-ATE65","GTEX-13N1W-0011-R4a-SM-AFUQV",
"GTEX-13NYB-0011-R1b-SM-9VIIM","GTEX-11GSP-0011-R7b-SM-9VIIR","GTEX-XV7Q-2526-SM-9VXYG","GTEX-13FTW-2626-SM-9VXZI","GTEX-13OW8-1426-SM-9VXZL",
"GTEX-U3ZN-2226-SM-9VXZT","GTEX-1211K-0826-SM-9VXYJ","GTEX-XV7Q-0326-SM-9VXYI","GTEX-1122O-2426-SM-9VXZJ","GTEX-1399R-0726-SM-9VXZ3",
"GTEX-U3ZN-0326-SM-9VXYQ","GTEX-YEC4-0526-SM-9VXYP","GTEX-U3ZN-2226-SM-9VXZT","GTEX-1211K-2026-SM-9VXZN","GTEX-13FTW-0626-SM-9VXYR",
"GTEX-13OW8-1426-SM-9VXZL","GTEX-11DXX-0326-SM-9VXZM")
gtex=gtex[,!sampleNames(gtex) %in% QC_fail]
tissues <- setNames(unique(gtex$Tissue), unique(gtex$Tissue))
gtex.cov=getCoverage(gtex)
colnames(gtex.cov)=colnames(gtex)
all<- which(rowSums2(gtex.cov >= 5) >= 100) #24630044 sites


# ------------------------------------------------------------------------------
# Compute means and sds
#
M=as.matrix(getMeth(gtex))
M2=M
gtex=gtex[all,]
M=M[all,]

list_of_sds <- mclapply(tissues, function(tissue) {
  rowSds(M[, gtex$Tissue == tissue])
})
list_of_means <- mclapply(tissues, function(tissue) {
  rowMeans2(M[, gtex$Tissue == tissue])
})
save(list_of_means, list_of_sds,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/vmr_stats.rda")

# HC1, HC2 'empirical' clusters (based on PCA)
HC1 <- c("SM-AFUP5", "SM-AFUL4", "SM-AFUNM", "SM-AFUOK", "SM-AFUKC",
         "SM-AFUKP", "SM-AFUNG", "SM-AFUOR", "SM-AFUJO", "SM-AFUJF",
         "SM-AFUJZ", "SM-AFUMU")
HC2 <- c("SM-AFUMH", "SM-AFUKV", "SM-AFUKJ", "SM-AFUK5",
         "SM-AFUMM", "SM-AFUNS")
clusters <- list(HC1 = HC1,
                 HC2 = HC2)
list_of_cluster_means <- mclapply(clusters, function(cluster) {
  rowMeans2(M[, gtex$Sample_ID %in% cluster])
})
list_of_cluster_sds <- mclapply(clusters, function(cluster) {
  rowSds(M[, gtex$Sample_ID %in% cluster])
})
save(list_of_cluster_means, list_of_cluster_sds,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/cluster_vmr_stats.rda")

# Combine 'tissue' and 'cluster' data
list_of_autosomal_means <- mclapply(
  c(list_of_means, list_of_cluster_means), function(x) {
    x[which(seqnames(gtex) %in% paste0("chr", 1:22))]
  })
list_of_autosomal_sds <- mclapply(
  c(list_of_sds, list_of_cluster_sds), function(x) {
    x[which(seqnames(gtex) %in% paste0("chr", 1:22))]
  })
save(list_of_autosomal_means, list_of_autosomal_sds,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/autosomal_vmr_stats.rda")

# ------------------------------------------------------------------------------
# Plot mean-variance relationship
#

# NOTE: Downsample autosomal data
groups <- setNames(names(list_of_autosomal_means),
                   names(list_of_autosomal_means))
group_colors <- sapply(groups, function(group) {
  if (group %in% gtex$Tissue) {
    unique(gtex$Brain_color[gtex$Tissue == group])
  } else {
    if (grepl("amygdala", group)) {
      unique(gtex$Brain_color[gtex$Tissue == "Brain - Amygdala"])
    } else if (grepl("HC", group)) {
      unique(gtex$Brain_color[gtex$Tissue == "Brain - Hippocampus"])
    } else {
      stop("Unknown group")
    }
  }
})
group_lty <- sapply(groups, function(group) {
  if (group %in% gtex$Tissue) {
    return(1L)
  } else {
    if (grepl("1", group)) {
      return(2L)
    } else if (grepl("2", group)) {
      return(3L)
    } else {
      stop("Unknown group")
    }
  }
})

autosomal_loess_fits <- mclapply(groups, function(group) {
    message(group)
    means <- list_of_autosomal_means[[group]]
    sds <- list_of_autosomal_sds[[group]]
    idx <- sample(seq_along(means), replace = FALSE, size = 10 ^ 6)
    fit <- lowess(means[idx], sds[idx], f = 1/10)
  })

pdf("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/mean-variance_relationship.loess_fits.autosomal.pdf")
plot(autosomal_loess_fits[[1]],
     col = group_colors[names(autosomal_loess_fits)[1]],
     type = "l",
     ylim = c(0, max(sapply(autosomal_loess_fits, "[[", 2))),
     lty = group_lty[names(autosomal_loess_fits)[1]])
lapply(names(autosomal_loess_fits[-1]), function(group) {
  lines(autosomal_loess_fits[[group]],
        col = group_colors[group],
        lty = group_lty[group])
})
legend("bottom",
       legend = groups,
       col = group_colors,
       lty = group_lty,
       bty = "n",
       cex = 0.9)
dev.off()

# ------------------------------------------------------------------------------
# Call VMRs
# First determine the 99% SD cutoff in each tissue 
vmrs2 <- function(gr, means, sds,
                  sdCut = quantile(sds, prob = 0.99, na.rm = TRUE)) {
  ## VMRs by thresholding sd
  isHigh <- rep(0, length(means))
  isHigh[sds > sdCut] <- 1
  vmrs <- bsseq:::regionFinder3(isHigh, as.character(seqnames(gr)),
                                start(gr), maxGap = 1000)$up
  vmrs
}

## What are the sd cutoffs used for each tissue?
list_of_sdCutoffs <- lapply(list_of_autosomal_sds,function(x){quantile(x, prob = 0.99, na.rm = TRUE)})
sd_cut=as.data.frame(unlist(list_of_sdCutoffs))

#### The lowest value was selected so to properly reflect overlaps among variable regions

## ID VMRs using 0.095 cutoff for all tissues so edit vmrs2 as below to reflect
vmrs2 <- function(gr, means, sds, sdCut = 0.095) {
  ## VMRs by thresholding sd
  isHigh <- rep(0, length(means))
  isHigh[sds > sdCut] <- 1
  vmrs <- bsseq:::regionFinder3(isHigh, as.character(seqnames(gr)),
                                start(gr), maxGap = 1000)$up
  vmrs
}

list_of_autosomal_vmrs.nq <- mcmapply(function(means, sds, gr.wgbs) {
    vmrs2(gr.wgbs, means, sds)
}, means = list_of_autosomal_means,
sds = list_of_autosomal_sds,
MoreArgs = list(
  gr.wgbs = keepSeqlevels(rowRanges(gtex), paste0("chr", 1:22), "coarse")),
SIMPLIFY = FALSE)
save(list_of_autosomal_vmrs.nq,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_autosomal_vmrs.nq.rda")


# ------------------------------------------------------------------------------
# Filter out VMRs if they are driven by outliers
#

filterVMRs <- function(vmrs_df, BSseq, cutoff = 0.7) {
  # NOTE: The maximum value of cooks.distance() is 1 when y = [0, 1] (based on
  #       empirical observations)
  stopifnot(cutoff < 1 & cutoff > 0)
  message("Computing perRegion meth")
  meth <- as.matrix(getMeth(BSseq, vmrs_df, what = "perRegion"))
  message("Identifying outlier samples")
  outliers <- lapply(seq_len(nrow(meth)), function(i) {
    idx <- cooks.distance(lm(meth[i, ] ~ 1)) > cutoff
    colnames(meth)[idx]
  })
  list(has_outliers = vmrs_df[lengths(outliers) > 0, ],
       no_outliers = vmrs_df[lengths(outliers) == 0, ])
}

list_of_filtered_autosomal_vmrs.nq <- mclapply(groups, function(group) {
  vmrs_df <- list_of_autosomal_vmrs.nq[[group]]
  vmrs_df <- vmrs_df[vmrs_df$n > 10, ]  # filter out vmrs with 10 CpGs or less
  if (group %in% gtex$Tissue) {
    BSseq <- gtex[, gtex$Tissue == group]
  } else if (group == "HC1") {
    BSseq <- gtex[, gtex$Sample_ID %in% HC1]
  } else if (group == "HC2") {
    BSseq <- gtex[, gtex$Sample_ID %in% HC2]
  } else {
    stop("Unknown group")
  }
  filterVMRs(vmrs_df, BSseq, cutoff = 0.7)
})

save(list_of_filtered_autosomal_vmrs.nq,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs.rda")

#### annotate with sds over each vmr

results=map2(list_of_filtered_autosomal_vmrs.nq,list_of_autosomal_sds,function(X,Y)
  sapply(1:nrow(X$no_outliers), function(i) mean(Y[c(X$no_outliers$idxStart[i]:X$no_outliers$idxEnd[i])])))
# now append these numbers onto the vmr list
list_of_filtered_autosomal_vmrs_wSDS=map2(list_of_filtered_autosomal_vmrs.nq,results,function(x,meanSDS) cbind(x$no_outliers,meanSDS))
save(list_of_filtered_autosomal_vmrs_wSDS,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs_wSDS.rda")

