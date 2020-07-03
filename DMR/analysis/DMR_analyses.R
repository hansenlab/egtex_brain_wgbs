
#### DMR HEATMAPS ####

library(EnrichedHeatmap)
library(SummarizedExperiment)
library(circlize)
library(viridis)
library(bsseq)
library(podkat)
library(HDF5Array)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
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

rm(BS,bsseq_filtered)

BG_samples=sampleNames(bsseq)[grep("basal ganglia",pData(bsseq)$Tissue)]
HC_samples=sampleNames(bsseq)[grep("Hippocampus",pData(bsseq)$Tissue)]

meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))


### NOTE FOR ALL HEATMAP FIGURES, run this on command line so can edit in illustrator
# convert -density 300 Random20k_general_DMRs_heatmap_dendrogram.pdf -quality 80 testing.png
## will only want to use the actual heatmap part and use the pdf version of the rest.

x=getMeth(bsseq,regions=list_of_all_DMRs$`basal_ganglia_CG-DMRs`,type="smooth",what="perRegion")

x=x[,BG_samples]
BG=bsseq[,BG_samples]
sex=pData(BG)$Self_Reported_Sex
tissue=pData(BG)$Tissue

ha = HeatmapAnnotation(type = tissue, sex=sex, 
    col = list(type = c("Brain - Caudate (basal ganglia)"="#666666",
"Brain - Nucleus accumbens (basal ganglia)"="#d95f02","Brain - Putamen (basal ganglia)"="#7570b3"),
    sex=c("Male"="blue","Female"="red")))

pdf("BG_DMRs_heatmap_dendrogram2.pdf")
Heatmap(as.matrix(x), name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = FALSE, show_column_names = T) 
dev.off()

bsseq_HC=bsseq[,HC_samples]
pData(bsseq_HC)$group <- factor(
  dplyr::case_when(
    pData(bsseq_HC)$Sample_ID %in% c("SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM",
                           "SM-AFUK5", "SM-AFUKJ", "SM-AFUJL") ~ "HC1",
    TRUE ~ "HC2"),
  levels = c("HC1", "HC2")
)

sex=pData(bsseq_HC)$Self_Reported_Sex
group=as.character(pData(bsseq_HC)$group)
ha = HeatmapAnnotation(type = group,sex=sex, 
    col = list(type = c("HC1"="#d95f02","HC2"="#7570b3"),sex=c("Male"="blue","Female"="red")))

HC=getMeth(bsseq_HC,regions=list_of_all_DMRs$`hippocampus_CG-DMRs`,type="smooth",what="perRegion")
HC2=HC[,HC_samples]
pdf("HC_DMRs_heatmap_dendrogram.pdf")
Heatmap(as.matrix(HC2), name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = FALSE, show_column_names = FALSE) 
dev.off()

#### Figure 1B
## Identify most discriminatory CG-DMRs among 5 tissue groups

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/Fstat_annotation_generalCG-DMRs.RData")
bg=dmrs[which(dmrs$BGvsAMY == TRUE & dmrs$CortvsBG ==TRUE & dmrs$HCvsBG ==TRUE & dmrs$BGvsHYP==TRUE & dmrs$CortvsHYP ==FALSE & dmrs$CortvsAMY ==FALSE & dmrs$HYPvsAMY==FALSE & dmrs$HCvsAMY ==FALSE & dmrs$HCvsHYP==FALSE&  dmrs$CortvsHC==FALSE),]  #25509
cort=dmrs[which(dmrs$CortvsHC == TRUE & dmrs$CortvsBG ==TRUE & dmrs$CortvsAMY ==TRUE & dmrs$CortvsHYP==TRUE & dmrs$BGvsHYP ==FALSE & dmrs$BGvsAMY ==FALSE & dmrs$HYPvsAMY==FALSE & dmrs$HCvsAMY ==FALSE & dmrs$HCvsHYP==FALSE&  dmrs$HCvsBG==FALSE),] #1604
amy=dmrs[which(dmrs$CortvsHC == FALSE & dmrs$CortvsBG ==FALSE & dmrs$CortvsAMY ==TRUE & dmrs$CortvsHYP==FALSE & dmrs$BGvsHYP ==FALSE & dmrs$BGvsAMY ==TRUE & dmrs$HYPvsAMY==TRUE & dmrs$HCvsAMY ==TRUE & dmrs$HCvsHYP==FALSE&  dmrs$HCvsBG==FALSE),] #14
hc=dmrs[which(dmrs$CortvsHC == TRUE & dmrs$CortvsBG ==FALSE & dmrs$CortvsAMY ==FALSE & dmrs$CortvsHYP==FALSE & dmrs$BGvsHYP ==FALSE & dmrs$BGvsAMY ==FALSE & dmrs$HYPvsAMY==FALSE & dmrs$HCvsAMY ==TRUE & dmrs$HCvsHYP==TRUE&  dmrs$HCvsBG==TRUE),] #14158
hyp=dmrs[which(dmrs$CortvsHC == FALSE & dmrs$CortvsBG ==FALSE & dmrs$CortvsAMY ==FALSE & dmrs$CortvsHYP==TRUE & dmrs$BGvsHYP ==TRUE & dmrs$BGvsAMY ==FALSE & dmrs$HYPvsAMY==TRUE & dmrs$HCvsAMY ==FALSE & dmrs$HCvsHYP==TRUE&  dmrs$HCvsBG==FALSE),] #2006

#since bg so many....sort by maxStat and choose top 10k
bg2=sort(bg,by=~maxStat,decreasing=TRUE)
test=c(bg2[1:10000,],cort,amy,hc,hyp)

### now make heatmap of these

ha = HeatmapAnnotation(type = pData(bsseq)$Tissue, group=pData(bsseq)$group, 
    col = list(type = c("Brain - Amygdala"="#e7298a","Brain - Anterior cingulate cortex (BA24)"="#a6761d","Brain - Caudate (basal ganglia)"="#666666",
"Brain - Frontal Cortex (BA9)"="#e6ab02","Brain - Hippocampus"="#1b9e77","Brain - Hypothalamus"="#66a61e","Brain - Nucleus accumbens (basal ganglia)"="#d95f02","Brain - Putamen (basal ganglia)"="#7570b3"),
  group=c("Amygdala"="#e7298a","Cortical"="#a6761d","Hippocampus"="#1b9e77","Basal_ganglia"="#7570b3","Hypothalamus"="#66a61e")))

test_meth=as.matrix(getMeth(bsseq,regions=test,type="smooth",what="perRegion"))

pdf("SpecificTissue_DMRs_heatmap_dendrogram2.pdf") # used in paper
Heatmap(as.matrix(test_meth), name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = FALSE, show_column_names = FALSE) 
dev.off()

### outside of R:
convert -density 300 Outlier_test_heatmap_dendrogram2_edit.pdf -quality 80 Outlier_test_heatmap_dendrogram2_edit.png
####


#### Figure 1E
## Plot example DMR
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis")

pData(bsseq)$group_color=as.character(pData(bsseq)$group)

bsseq$group_color[which(bsseq$group_color=="Cortical")]="#a6761d"
bsseq$group_color[which(bsseq$group_color=="Hippocampus")]="#1b9e77"
bsseq$group_color[which(bsseq$group_color=="Amygdala")]="#e7298a"
bsseq$group_color[which(bsseq$group_color=="Hypothalamus")]="#66a61e"
bsseq$group_color[which(bsseq$group_color=="Basal_ganglia")]="#7570b3"

#NPTXR gene chr22:38,822,609-38,843,858
NPTXR=GRanges(seqnames = "chr22", strand = "-",
              ranges = IRanges(start = 38822609,end=38843858))
NPTXRext=NPTXR+5000
x=subsetByOverlaps(dmrs,NPTXRext)
pdf(file = "NPTXR_General_CG-DMRs_plot.pdf")
plotManyRegions(BSseq = bsseq,
                regions = NPTXRext,col=pData(bsseq)$group_color,
                stat = "stat", addRegions = dmrs, geneTrack=exon,extend = 5000,
                stat.ylim = c(-10, 10))
dev.off()

### get expression of NPTXR from GTEx for Fig 1F

es=readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/GTEx_v8_HC_RNAseq_tpm.rds")
brain_tissues=c("Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)","Brain - Caudate (basal ganglia)",
  "Brain - Hippocampus","Brain - Hypothalamus","Brain - Amygdala","Brain - Frontal Cortex (BA9)","Brain - Anterior cingulate cortex (BA24)")
names=pData(es)$SAMPID[which(pData(es)$SMTSD %in% brain_tissues)]
es_brain<-es[,which(colnames(es) %in% names)]

pData(es_brain)$group="other"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Hippocampus"]="Hippocampus"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Hypothalamus"]="Hypothalamus"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Amygdala"]="Amygdala"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Frontal Cortex (BA9)"]="Cortical"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Anterior cingulate cortex (BA24)"]="Cortical"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Nucleus accumbens (basal ganglia)"]="Basal_ganglia"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Putamen (basal ganglia)"]="Basal_ganglia"
pData(es_brain)$group[pData(es_brain)$SMTSD=="Brain - Caudate (basal ganglia)"]="Basal_ganglia"

expression=exprs(es)


NPTXR=expression[grep("ENSG00000221890",rownames(expression)),]
NPTXR2=as.data.frame(NPTXR)
### Now just get the one gene NPTXR
names=pData(es_brain)$SAMPID
NPTXR_3=as.data.frame(NPTXR2[rownames(NPTXR2) %in% names,])
NPTXR_3$sample=rownames(NPTXR2)[rownames(NPTXR2) %in% names]
NPTXR_3$group=pData(es_brain)$group[match(pData(es_brain)$SAMPID ,NPTXR_3$sample)]
len=gen_len[match(rownames(test), rownames(gen_len)), ]  # returns the length value of genes in the right order

colnames(NPTXR_3)=c("value","sample")
pdf("NPTXR_expression.pdf")  # used in Figure 1
theme_bw()
ggplot(NPTXR_3, aes(x=group, y=value,fill=group)) + 
  geom_violin(trim=F)+geom_boxplot(width=0.1)

ggplot(NPTXR_3, aes(x=group, y=value, fill=group)) + 
  geom_boxplot() 
dev.off()

#### Figure 3E (Venn diagrams)

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/Reduced_CH-DMRs.rda") ##general_CH_DMRs,hippocampus_CH_DMRs,basal_ganglia_CH_DMRs
gen_CG_DMRs=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/general_CG-DMRs.rds"))
BG_CG_DMRs=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/basal_ganglia_CG-DMRs.rds"))
HC_CG_DMRs=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/hippocampus_CG-DMRs.rds"))

library(Vennerable)
library(ChIPpeakAnno)
all_overlaps2=findOverlapsOfPeaks(general_CH_DMRs,hippocampus_CH_DMRs,basal_ganglia_CH_DMRs,maxgap=10,connectedPeaks="merge")
all_overlaps=findOverlapsOfPeaks(gen_CG_DMRs,BG_CG_DMRs,HC_CG_DMRs,maxgap=10,connectedPeaks="merge")
all_overlaps3=findOverlapsOfPeaks(general_CH_DMRs,gen_CG_DMRs,maxgap=10,connectedPeaks="merge")



pdf("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/DMR_Venn_Diagrams2.pdf")
makeVennDiagram(all_overlaps2, NameOfPeaks=c("All","Hippocampus","Basal ganglia"), totalTest=500000,
                main="CH-DMRs", main.fontfamily="sans", 
                sub.fontfamily="sans", cat.fontfamily="sans", col = "transparent",  fontfamily = "sans",
                fill = c("darkorange1", "seagreen3", "orchid3"))
makeVennDiagram(all_overlaps, NameOfPeaks=c("All","Basal ganglia","Hippocampus"), totalTest=500000,
                main="CG-DMRs", main.fontfamily="sans", 
                sub.fontfamily="sans", cat.fontfamily="sans", col = "transparent",  fontfamily = "sans",
                fill = c("darkorange1", "seagreen3", "orchid3"))

makeVennDiagram(all_overlaps3, NameOfPeaks=c("CH-DMRs","CG-DMRs"), totalTest=500000,
                main="CG-DMRs that overlap CH-DMRs", main.fontfamily="sans", 
                sub.fontfamily="sans", cat.fontfamily="sans", col = "transparent",  fontfamily = "sans",
                fill = c("darkorange1", "seagreen3"))

dev.off()
### because of the way peaks are combined and counted manually adjusted the CH-DMR vs CG-DMR venn diagram numbers
a=intersect(gen_CG,general_CH_DMRs)
b=subsetByOverlaps(gen_CG,general_CH_DMRs,invert=TRUE)
c=subsetByOverlaps(general_CH_DMRs,gen_CG,invert=TRUE)





