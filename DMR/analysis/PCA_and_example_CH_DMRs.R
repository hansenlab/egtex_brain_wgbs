# Plot some example CH-DMRs and generate PCA plot
# Peter Hickey
# 2017-12-19
# Edited by Lindsay Rizzardi
# 6/26/19


library(bsseq)
library(matrixStats)
library(limma)
library(scales)
library(dplyr)
library(bsseq)
library(HDF5Array)
library(DelayedMatrixStats)
library(dplyr)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

options("mc.cores" = 15)
#options("DelayedArray.block.size" = DelayedArray:::DEFAULT_BLOCK_SIZE * 100L)

extdir <- "/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq"

### ----------------------------------------------------------------------------
### Load data
###

CA_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth"))
CA_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "eGTEx.Phase2_brain_samples.mCA_pos.small_smooth"))
CT_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "eGTEx.Phase2_brain_samples.mCT_pos.small_smooth"))
CT_neg_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth"))

strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)
pretty_names <- paste0("m", contexts, " (",
                       ifelse(strands == "pos", "+", "-"), ")")

# ------------------------------------------------------------------------------
# mCH
#

CH_BSseq_names <- c("pos_CA", "neg_CA", "pos_CT", "neg_CT")
names(CH_BSseq_names) <-  paste0("m", contexts, " (",
                                 ifelse(strands == "pos", "+", "-"), ")")
list_of_CH_BSseq <- list(CA_pos_BSseq,CA_neg_BSseq,CT_pos_BSseq,CT_neg_BSseq)
names(list_of_CH_BSseq)=names(CH_BSseq_names)

# ------------------------------------------------------------------------------
# mCG
#

CG_small_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "eGTEx.Phase2_brain_samples.mCG.small_smooth"))

CG_large_BSseq <- loadHDF5SummarizedExperiment(
  dir = file.path(extdir, "eGTEx.Phase2_brain_samples.mCG.large_smooth"))

list_of_CG_BSseq <- list("mCG (S)" = CG_small_BSseq,
                         "mCG (L)" = CG_large_BSseq)

# ------------------------------------------------------------------------------
# mC
#

list_of_BSseq <- c(list_of_CG_BSseq, list_of_CH_BSseq)
# NOTE: Fix up colours
c("Brain - Hippocampus"="#1b9e77","Brain - Amygdala" ="#e7298a","Brain - Putamen (basal ganglia)" = "#7570b3",
  "Brain - Hypothalamus" = "#66a61e", "Brain - Nucleus accumbens (basal ganglia)" = "#d95f02", 
  "Brain - Anterior cingulate cortex (BA24)" = "#a6761d","Brain - Frontal Cortex (BA9)" = "#e6ab02", "Brain - Caudate (basal ganglia)" = "#666666"))


#"Cortical"="#a6761d"
#"Basal_ganglia"="#7570b3"

list_of_BSseq <- lapply(list_of_BSseq, function(BSseq) {
  BSseq$col <- BSseq$Brain_color
  BSseq
})
col=list_of_BSseq[[1]]$Brain_color

col2=list_of_BSseq[[1]]$Tissue
col2[which(col2=="Brain - Anterior cingulate cortex (BA24)")]="#a6761d"
col2[which(col2=="Brain - Frontal Cortex (BA9)")]="#a6761d"
col2[which(col2=="Brain - Hippocampus")]="#1b9e77"
col2[which(col2=="Brain - Amygdala")]="#e7298a"
col2[which(col2=="Brain - Putamen (basal ganglia)")]="#7570b3"
col2[which(col2=="Brain - Caudate (basal ganglia)")]="#7570b3"
col2[which(col2=="Brain - Nucleus accumbens (basal ganglia)")]="#7570b3"
col2[which(col2=="Brain - Hypothalamus")]="#66a61e"

# ------------------------------------------------------------------------------
# DMRs and blocks
#
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/Reduced_CH-DMRs.rda") #general_CH_DMRs,hippocampus_CH_DMRs,basal_ganglia_CH_DMRs
CT_neg_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CT_neg-DMRs.rds")
CA_neg_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CA_neg-DMRs.rds")
CT_pos_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CT_pos-DMRs.rds")
CA_pos_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CA_pos-DMRs.rds")
gen_CG_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/general_CG-DMRs.rds")
gen_CG_blocks=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/general_CG-blocks.rds")

list_of_CH_DMRs=list(CT_neg_DMRs,CA_neg_DMRs,CT_pos_DMRs,CA_pos_DMRs)
names(list_of_CH_DMRs)=c("mCT (-)","mCA (-)","mCT (+)","mCA (+)")
list_of_CH_DMRs=lapply(list_of_CH_DMRs,function(x) GRanges(x))

DMRs_and_blocks <- c(list("mCG (S)" = makeGRangesFromDataFrame(gen_CG_DMRs),
                          "mCG (L)" = makeGRangesFromDataFrame(gen_CG_blocks)),
                     list_of_CH_DMRs)

# ------------------------------------------------------------------------------
# Gene models for plots
#

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/Gencodev26_exons_plotting.rda") #exon

### ============================================================================
### Plots
###

source("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/EDA_functions.R")

# ------------------------------------------------------------------------------
# Marker-gene-focused 
# RASGRF2, RGS4,NRGN,PAX6,CACNA1A,FEZF2,GAD65,CAMK2A,HTR2A and HTR4
#
library(biomaRt)
mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")
attributes <- listAttributes(mart.hs)
x=getBM(attributes=c('hgnc_symbol','chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=c("RASGRF2","RGS4","NRGN","PAX6","CACNA1A","SP8","PROX1","FEZF2","SLC17A6","CAMK2A","HTR2A","HTR4"),
      mart=mart.hs)
colnames(x)=c("hgnc_symbol","chr","start","end","strand")
x$chr=paste0("chr",x$chr)
x$strand[which(x$strand==-1)]="-"
x$strand[which(x$strand==1)]="+"

x_gr=GRanges(x)
# Find DEGs where at least 50% of gene is covered by a CA_pos-DMR

regions <- subsetByOverlaps(x_gr,DMRs_and_blocks[["mCA (+)"]])
T=lapply(list_of_CH_BSseq,function(x) subsetByOverlaps(x,regions))

ca=subsetByOverlaps(CA_pos_BSseq,regions)

test2=lapply(list_of_CH_BSseq,function(x) subsetByOverlaps(x_gr,x))
test=lapply(list_of_CH_BSseq,function(x) subsetByOverlaps(x,x_gr))
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/")

pdf("MarkerGenes_covered_by_CA_pos-DMRs5.pdf")
plotManyRegions2(
  list_of_CG_BSseq = list_of_CG_BSseq,
  list_of_CH_BSseq = list_of_CH_BSseq,
  #list_of_CA_BSseq = list_of_CH_BSseq[c(1,3)],
  #list_of_CT_BSseq = list_of_CH_BSseq[c(2,4)],
  col=col2,
  regions = regions,
  list_of_CG_addRegions = DMRs_and_blocks[c("mCG (S)", "mCG (L)")],
  list_of_CH_addRegions = DMRs_and_blocks[names(list_of_CH_DMRs)],
  #list_of_CT_addRegions = DMRs_and_blocks[c("mCT (-)","mCT (+)")],
  #list_of_CA_addRegions = DMRs_and_blocks[c("mCA (-)","mCA (+)")],
  ATAC_SE = NULL,
  geneTrack = exon,
  extend = 10000)
dev.off()


##########  PCA OF mCH  #########
library(bsseq)
library(matrixStats)
library(viridis)
library(gplots)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(parallel)
options("scipen" = 20)


#extdir <- "../extdata"
strands <- rep(c("pos", "neg"), each = 2)
contexts <- rep(c("CA", "CT"), times = 2)

list_of_CH_BSseq <- lapply(list_of_CH_BSseq, function(bsseq) {
  bsseq$col <- col2
  bsseq
})
### ============================================================================
### Per-sample average methylation in bins
###

w <- 1000
bins <- tileGenome(
  seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr", 1:22)],
  tilewidth = w,
  cut.last.tile.in.chrom = TRUE)
# Add dinucleotide counts of each bin
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, bins)
dns <- c("AA", "AC", "AG", "AT",
         "CA", "CC", "CG", "CT",
         "GA", "GC", "GG", "GT",
         "TA", "TC", "TG", "TT")
names(dns) <- dns
dn_counts <- lapply(X = dns,
                    FUN = vcountPattern,
                    subject = seq)
mcols(bins) <- as(dn_counts, "DataFrame")

# Save bins
saveRDS(bins, paste0("bins.", w, "kb.rds"))

list_of_meanMeth <- mclapply(list_of_CH_BSseq, function(BSseq) {
  getMeth(BSseq = BSseq, regions = bins, type = "smooth",
          what = "perRegion")
},mc.cores=10)

bsseq=list_of_CH_BSseq[[1]]
j=1
meth <- do.call(
  cbind, 
  lapply(seq_len(ncol(bsseq)), function(j) {
    message(j)
    # Subset to just the j-th sample
    tmp <- bsseq[, j]
    # Load the data for the j-th sample into memory
    assays(tmp, withDimnames = FALSE) <- endoapply(
      assays(tmp, withDimnames = FALSE),
      as.matrix)
    getMeth(
      tmp,
      regions = bins,
      type = "smooth",
      what = "perRegion",
      withDimnames = FALSE)
  }))
colnames(meth) <- colnames(bsseq)
# meth should be a matrix of meth values in 1kb bins (rows) for each sample 
# (columns)

saveRDS(meth,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/mCAplus_meanMeth.1kb_bins.rds")

### TRY FOR  A SINGLE strand/CONTEXT
meanMeth <- meth[!rowAlls(meth, value = NA), ]
BiocManager:::install("Morpho")
library(Morpho)
BiocManager:::install("caret")
library(caret)
library(foreach)
zero=nearZeroVar(meanMeth,allowParallel=TRUE) #yeilded nothing so no zeroVar regions

pca=prcompfast(t(meanMeth))
percent_variance <- setNames(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2),
                               colnames(pca$x))
df <- data.frame(x = pca$x[, "PC1"],y = pca$x[, "PC2"])
df2=df*-1
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PCA")
pdf("eGTEx_Sorted_PCA_mCAplus_1kb_2.pdf",width=12,height=10)

theme_set(theme_bw(base_size=16))
g <- ggplot(data = df2, aes(x = x, y = y, colour = pData(bsseq)$Tissue)) + coord_fixed(1)+
      geom_point(size=3.5)+scale_colour_manual(values = c("Brain - Anterior cingulate cortex (BA24)"="#a6761d","Brain - Frontal Cortex (BA9)"="#a6761d","Brain - Hippocampus"="#1b9e77","Brain - Amygdala"="#e7298a","Brain - Putamen (basal ganglia)"="#7570b3","Brain - Caudate (basal ganglia)"="#7570b3","Brain - Nucleus accumbens (basal ganglia)"="#7570b3","Brain - Hypothalamus"="#66a61e"))+
      labs(colour="Tissue")
g + xlab(paste0("PC1", " (", round(percent_variance["PC1"], 1), " %)")) +
    ylab(paste0("PC2", " (", round(percent_variance["PC2"], 1), " %)")) 

h <- ggplot(data = df2, aes(x = x, y = y, colour = pData(bsseq)$Tissue)) +coord_fixed(1)+
    geom_point(size=3.5)+scale_colour_manual(values = c("Brain - Hippocampus"="#1b9e77","Brain - Amygdala" ="#e7298a","Brain - Putamen (basal ganglia)" = "#7570b3",
  "Brain - Hypothalamus" = "#66a61e", "Brain - Nucleus accumbens (basal ganglia)" = "#d95f02", 
  "Brain - Anterior cingulate cortex (BA24)" = "#a6761d","Brain - Frontal Cortex (BA9)" = "#e6ab02", "Brain - Caudate (basal ganglia)" = "#666666"))+
    labs(colour="Tissue")
h + xlab(paste0("PC1", " (", round(percent_variance["PC1"], 1), " %)")) +
    ylab(paste0("PC2", " (", round(percent_variance["PC2"], 1), " %)")) 

dev.off()

