### Plot heatmap of CH-DMRs with greatest maxDiff (>0.1)
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/")
library(SummarizedExperiment)
library(dplyr)
library(scales)
library(gplots)
library(matrixStats)
library(limma)
library(dplyr)
library(bsseq)
library(HDF5Array)
library(DelayedMatrixStats)
library(dplyr)
library(circlize)
library(EnrichedHeatmap)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)
## Make list of CH-DMRs
# General CH-DMRs
CAneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CA_neg-DMRs.rds"))
CApos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CA_pos-DMRs.rds"))
CTneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CT_neg-DMRs.rds"))
CTpos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CT_pos-DMRs.rds"))
CCneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CC_neg-DMRs.rds"))
CCpos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CC_pos-DMRs.rds"))
# BG CH-DMRs
BG_CAneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/basal_ganglia_CA_neg-DMRs.rds"))
BG_CApos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/basal_ganglia_CA_pos-DMRs.rds"))
BG_CTneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/basal_ganglia_CT_neg-DMRs.rds"))
BG_CTpos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/basal_ganglia_CT_pos-DMRs.rds"))
BG_CCneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/basal_ganglia_CC_neg-DMRs.rds"))
BG_CCpos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/basal_ganglia_CC_pos-DMRs.rds"))

#HC CH-DMRs
HC_CAneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/hippocampus_CA_neg-DMRs.rds"))
HC_CApos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/hippocampus_CA_pos-DMRs.rds"))
HC_CTneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/hippocampus_CT_neg-DMRs.rds"))
HC_CTpos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/hippocampus_CT_pos-DMRs.rds"))
HC_CCneg=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/hippocampus_CC_neg-DMRs.rds"))
HC_CCpos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/hippocampus_CC_pos-DMRs.rds"))

list_of_CH_DMRs<-GenomicRangesList(CAneg,CApos,CTneg,CTpos,CCneg,CCpos,BG_CAneg,BG_CApos,BG_CTneg,BG_CTpos,BG_CCneg,BG_CCpos,HC_CAneg,HC_CApos,HC_CTneg,HC_CTpos,HC_CCneg,HC_CCpos)
names(list_of_CH_DMRs)=c("general mCA (-)","general mCA (+)","general mCT (-)","general mCT (+)","general mCC (-)","general mCC (+)",
  "basal ganglia mCA (-)","basal ganglia mCA (+)","basal ganglia mCT (-)","basal ganglia mCT (+)","basal ganglia mCC (-)","basal ganglia mCC (+)",
  "hippocampus mCA (-)","hippocampus mCA (+)","hippocampus mCT (-)","hippocampus mCT (+)","hippocampus mCC (-)","hippocampus mCC (+)")


### Make general CH-DMRs (exactly like was done with SLDSR script)
general_CH_DMRs=reduce(Reduce(c,list_of_CH_DMRs[1:6]),ignore.strand=T)
basal_ganglia_CH_DMRs=reduce(Reduce(c,list_of_CH_DMRs[7:12]),ignore.strand=T)
hippocampus_CH_DMRs=reduce(Reduce(c,list_of_CH_DMRs[13:18]),ignore.strand=T)
save(basal_ganglia_CH_DMRs,general_CH_DMRs,hippocampus_CH_DMRs,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/Reduced_CH-DMRs.rda")


CA_pos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CA_pos-DMRs.rds"))
CA_pos_max=CA_pos[which((CA_pos$maxDiff) > 0.1),] #3195

plot_these=subsetByOverlaps(general_CH_DMRs,CA_pos_max) #3060



extdir <- "/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq"

### ----------------------------------------------------------------------------
### Load data
###

CA_pos_BSseq <- loadHDF5SummarizedExperiment(
  file.path(extdir, "eGTEx.Phase2_brain_samples.mCA_pos.small_smooth"))
CA_pos_BSseq$group <- factor(
  dplyr::case_when(
    CA_pos_BSseq$Tissue == "Brain - Frontal Cortex (BA9)" ~ "Cortical",
    CA_pos_BSseq$Tissue ==
      "Brain - Anterior cingulate cortex (BA24)" ~ "Cortical",
    CA_pos_BSseq$Tissue == "Brain - Hippocampus" ~ "Hippocampus",
    CA_pos_BSseq$Tissue == "Brain - Amygdala" ~ "Amygdala",
    CA_pos_BSseq$Tissue == "Brain - Hypothalamus" ~ "Hypothalamus",
    CA_pos_BSseq$Tissue ==
      "Brain - Nucleus accumbens (basal ganglia)" ~ "Basal_ganglia",
    CA_pos_BSseq$Tissue ==
      "Brain - Putamen (basal ganglia)" ~ "Basal_ganglia",
    CA_pos_BSseq$Tissue ==
      "Brain - Caudate (basal ganglia)" ~ "Basal_ganglia"),
  levels = c(
    "Cortical", "Hippocampus", "Amygdala", "Hypothalamus", "Basal_ganglia")
)

### now make heatmap of these

ha = HeatmapAnnotation(type = pData(CA_pos_BSseq)$Tissue, group=pData(CA_pos_BSseq)$group, 
    col = list(type = c("Brain - Amygdala"="#e7298a","Brain - Anterior cingulate cortex (BA24)"="#a6761d","Brain - Caudate (basal ganglia)"="#666666",
"Brain - Frontal Cortex (BA9)"="#e6ab02","Brain - Hippocampus"="#1b9e77","Brain - Hypothalamus"="#66a61e","Brain - Nucleus accumbens (basal ganglia)"="#d95f02","Brain - Putamen (basal ganglia)"="#7570b3"),
  group=c("Amygdala"="#e7298a","Cortical"="#a6761d","Hippocampus"="#1b9e77","Basal_ganglia"="#7570b3","Hypothalamus"="#66a61e")))

test_meth=as.matrix(getMeth(CA_pos_BSseq,regions=plot_these,type="smooth",what="perRegion"))
# quantile(test_meth)
#           0%          25%          50%          75%         100% 
# 3.206627e-10 1.059533e-01 1.516009e-01 1.955165e-01 9.693957e-01 
## Make the scale fit the data
test_meth2[test_meth2>0.2]=0.2
meth_col_fun = colorRamp2(c(0, 0.2), c("white", "red"))
meth_col_fun = colorRamp2(c(0, 0.2), c("blue","white"))

pdf("Top_CA_DMRs_heatmap2.pdf")
Heatmap(as.matrix(test_meth2), name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = FALSE, show_column_names = FALSE) 
dev.off()
