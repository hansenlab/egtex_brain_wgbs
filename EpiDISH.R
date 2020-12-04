##################
# EpiDISH Analysis
# Perform analysis of cell proportions using our previous NeuN sorted data from Kozlenkov et al 2014
##################

library(bsseq)
library(HDF5Array)
library(DelayedMatrixStats)
library(matrixStats)
library(DelayedArray) 
library(rtracklayer)
library(EpiDISH)
library(ggplot2)
library(EnrichedHeatmap)
library(openxlsx)
library(dplyr)
library(purrr)
library(tidyr)
library(gplots)
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis")

Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)
BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
BS=keepStandardChromosomes(BS,pruning.mode="coarse")
BS=dropSeqlevels(BS,c("chrX","chrY","chrM"),pruning.mode="coarse")
BS=sort(sortSeqlevels(BS))
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/QC_failures.rda")
bsseq_filtered <- BS[, !colnames(BS) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]
BS <- bsseq_filtered
rm(bsseq_filtered)
# 131 samples


### Stella's 2014 data
glu1=read.csv("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/KozlenkovNAR2014_HM450K_glialDMsites.csv",header=T)
glu2=read.csv("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/KozlenkovNAR2014_HM450K_neuronalDMsites.csv",header=T)
glu1=rbind(glu1,glu2)
glu=glu1[which(abs(glu1$DeltaBeta.G.N.) > 0.7),] # 426 sites

glu=glu[,c(2,3,5,6)]
glu$CHR=paste0("chr",glu$CHR)
colnames(glu)=c("chr","start","Glia","Neuron")
glu$end=glu$start
glu=GRanges(glu)
### liftOver to hg38
ch = import.chain("/amber3/feinbergLab/personal/lrizzard/genomes/hg19ToHg38.over.chain")
glu_hg38=unlist(liftOver(glu,ch))
Stell_posvneg=glu_hg38


#### Find Overlaps with our 8-group DMRs and HC DMRs
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/Fstat_annotation_OVERALL_CG-DMRs.RData")
egtex_dmrs=dmrs   
hc_dmrs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/hippocampus_CG-DMRs.rds")
hc_dmrs=GRanges(hc_dmrs)

x2=subsetByOverlaps(Stell_posvneg,egtex_dmrs,invert=T) 
x2=subsetByOverlaps(x2,hc_dmrs,invert=T) 
Stella.ref.M=cbind(mcols(x2)[,1],mcols(x2)[,2])
colnames(Stella.ref.M)=c("NeuN_neg","NeuN_pos")

gtex.M.Stella=getMeth(BS,regions=x2,type="smooth",what="perRegion",withDimnames=FALSE)
colnames(gtex.M.Stella)=sampleNames(BS)

rownames(x2)=NULL
ref_df=as.data.frame(x2,row.names=NULL)
ref_names=paste(ref_df$seqnames,ref_df$start,sep="_")

rownames(gtex.M.Stella)=ref_names
rownames(Stella.ref.M)=ref_names


# Run EpiDISH
out.l <- epidish(gtex.M.Stella,Stella.ref.M, method = "RPC")
out.l$estF  #estimated cell fraction matrix
est=out.l$estF 
e=as.data.frame(est)
t=pData(BS)$Tissue[which(pData(BS)$Sample_ID_for_Data_Sharing_and_Public_Release %in% rownames(e))]
e$Tissue=t


pdf("EpiDISH.pdf")
ggplot(e, aes(x=Tissue, y=NeuN_pos,fill=Tissue)) + theme_classic() + ylim(0,1)+
  geom_boxplot() + geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




