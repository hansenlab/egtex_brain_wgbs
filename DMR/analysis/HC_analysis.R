###### Understanding 2 HC groups ###########
library(bsseq)
BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
HC=BS[,which(pData(BS)$Tissue=="Brain - Hippocampus")]
QC_fail=c("GTEX-13G51-0011-R6b-SM-AFUJX","GTEX-13G51-0011-R10b-SM-AFUJU","GTEX-13N1W-0826-SM-ATE4J","GTEX-13O3Q-0011-R3a-SM-AFUKX",
"GTEX-T2IS-0011-R10a-SM-AFUOS","GTEX-13N1W-0011-R5b-SM-AFUQW","GTEX-13N1W-0011-R1b-SM-AFUQX","GTEX-13N1W-0011-R7b-SM-AFUQT",
"GTEX-13O3Q-0011-R6a-SM-AFUKS","GTEX-WZTO-0426-SM-ATE6F","GTEX-13N1W-0726-SM-ATE66","GTEX-13FHO-1026-SM-ATE65","GTEX-13N1W-0011-R4a-SM-AFUQV",
"GTEX-13NYB-0011-R1b-SM-9VIIM","GTEX-11GSP-0011-R7b-SM-9VIIR","GTEX-XV7Q-2526-SM-9VXYG","GTEX-13FTW-2626-SM-9VXZI","GTEX-13OW8-1426-SM-9VXZL",
"GTEX-U3ZN-2226-SM-9VXZT","GTEX-1211K-0826-SM-9VXYJ","GTEX-XV7Q-0326-SM-9VXYI","GTEX-1122O-2426-SM-9VXZJ","GTEX-1399R-0726-SM-9VXZ3",
"GTEX-U3ZN-0326-SM-9VXYQ","GTEX-YEC4-0526-SM-9VXYP","GTEX-U3ZN-2226-SM-9VXZT","GTEX-1211K-2026-SM-9VXZN","GTEX-13FTW-0626-SM-9VXYR",
"GTEX-13OW8-1426-SM-9VXZL","GTEX-11DXX-0326-SM-9VXZM")
x=sampleNames(HC) %in% QC_fail
HC=HC[,which(x==FALSE)]
pData(HC)$HC_group="HC2"
pData(HC)$HC_group <- factor(
  dplyr::case_when(
    pData(HC)$Sample_ID %in% c("SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM",
                           "SM-AFUK5", "SM-AFUKJ", "SM-AFUJL") ~ "HC2",
    TRUE ~ "HC1"),
  levels = c("HC2", "HC1")
)

### Plot DMRs in PROX1
PROX1=GRanges(seqnames="chr1",IRanges(start=213984935, end=214043502))
dmrs=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/hippocampus_CG-DMRs.rds"))
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/Gencodev26_exons_plotting.rda")

pdf("Prox1_HC_DMRs.pdf")
plotRegion(HC,region=PROX1,col=HC$HC_group,addRegions=dmrs,geneTrack=exon,extend=5000)
dev.off()


### 66 Granule cell Marker genes from Figure2A heatmap
genes=read.table("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/elife-14997-supp1-v1.txt",sep="\t",header=T)

genes2=genes[which(genes$enriched %in% c("dg_d-dg_v","ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v","ca3_d-ca3_v-ca2-ca1_d-ca1_v") & genes$depleted %in% c("dg_d-dg_v","ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v","ca3_d-ca3_v-ca2-ca1_d-ca1_v")),]
musGenes=genes2$gene_short_name
#16 more granule cell genes from Table_S1 PMID: 29204516
# Add more from PMID: 32203495 "Kcnk1","Camk1","Gabrd"
musGenes=unique(c(as.character(genes2$gene_short_name),"Kcnk1","Camk1","Gabrd","Cdhr1","Dsg2","Dsp","Ehd1","Neurod1","Npnt","Npy5r","Ogn","Spint2","St8sia4","Stx3","Tdo2","Trpc6","Eml5","Pter")) 


converted_genes <- convertMouseGeneList(musGenes)
converted_genes$Strand[which(converted_genes$Strand<0)]="-"
converted_genes$Strand[which(converted_genes$Strand>0)]="+"
colnames(converted_genes)=c("MGI_symbol","chr","start","end","strand","HGNC_symbol","Gene_ID")
converted_genes$chr=paste0("chr",converted_genes$chr)
converted_genes2=GRanges(converted_genes)
## add upstream 2000bp for promoter
converted_genes3=resize(converted_genes2,width(converted_genes2)+2000,fix="end")

### Do the same thing but with just the genes overlapping the HC DMRs
dg_dmrs=subsetByOverlaps(dmrs,converted_genes3) #116

### Save data as table
library(openxlsx)
library(readr)

wb=createWorkbook()
addWorksheet(wb,"Granule_marker_gene_list")
addWorksheet(wb,"CG-DMRs_Granule_marker_genes")

writeDataTable(wb,"Granule_marker_gene_list",x=as.data.frame(converted_genes3,stringsAsFactors=F))
writeDataTable(wb,"CG-DMRs_Granule_marker_genes",x=as.data.frame(dg_dmrs,stringsAsFactors=F))

saveWorkbook(wb, "SuppTableX–DG Marker Genes-DMRs.xlsx", overwrite = TRUE)
library(EnrichedHeatmap)


x=findOverlaps(dmrs,converted_genes3)
ol_gene=converted_genes3$HGNC_symbol[subjectHits(x)]

meth_all_dg=as.matrix(getMeth(bsseq,region=dg_dmrs,type="smooth",what="perRegion"))
rownames(meth_all_dg)=ol_gene
new_meth=meth_all_dg[order(rownames(meth_all_dg)), ] 


pData(bsseq)$HC_group="none"
pData(bsseq)$HC_group <- factor(
  dplyr::case_when(
    pData(bsseq)$Sample_ID %in% c("SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM",
                           "SM-AFUK5", "SM-AFUKJ", "SM-AFUJL") ~ "none",
    TRUE ~ "HC1"),
  levels = c("none", "HC1")
)
ha = HeatmapAnnotation(type = Tissue, class=pData(bsseq)$HC_group,
    col = list(class=c("none"="grey","HC1"="lightblue"),type = c("Brain - Amygdala"="#e7298a","Brain - Anterior cingulate cortex (BA24)"="#a6761d","Brain - Caudate (basal ganglia)"="#666666",
"Brain - Frontal Cortex (BA9)"="#e6ab02","Brain - Hippocampus"="#1b9e77","Brain - Hypothalamus"="#66a61e","Brain - Nucleus accumbens (basal ganglia)"="#d95f02","Brain - Putamen (basal ganglia)"="#7570b3")))

pdf("HC_DMRs_DG_MarkerGenes_ALL_TISSUES_heatmap.pdf")
Heatmap(as.matrix(meth_all_dg), name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),show_column_names = F) 
dev.off()

pdf("HC_DMRs_DG_MarkerGenes_ALL_TISSUES_heatmap2.pdf")
Heatmap(as.matrix(new_meth), name = "methylation", col = meth_col_fun, top_annotation = ha, top_annotation_height = unit(4, "mm"),  
    show_row_names = TRUE, cluster_rows = FALSE,row_names_gp = gpar(fontsize = 6),show_column_names = F) 
dev.off()



