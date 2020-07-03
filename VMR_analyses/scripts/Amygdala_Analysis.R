## Amygdala Variability analysis

library(bsseq)
library(HDF5Array)
library(DelayedMatrixStats)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(dplyr)
library(stringr)
library(tidyr)
library(gplots)
library(EnrichedHeatmap)
library(SummarizedExperiment)
library(circlize)
library(viridis)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)
# read in TableS2 from PMID 29024657
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs_wSDS.rda")
library(openxlsx)
gene_list=read.xlsx("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/TableS2_PMID_29024657.xlsx",2,startRow=5)
convertMouseGeneList <- function(x){
 
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("chromosome_name","start_position","end_position","strand","hgnc_symbol","ensembl_gene_id"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
 
# Print the first 6 genes found to the screen
print(head(genesV2))
return(genesV2)
}
genes=as.factor(gene_list[,1])

converted_genes <- convertMouseGeneList(genes)
converted_genes$Chromosome.scaffold.name=paste0("chr",converted_genes$Chromosome.scaffold.name)
colnames(converted_genes)=c("mgi_symbol","chr","start","end","strand","symbol","gene_id")
converted_genes$strand="*"
converted_genes=GRanges(converted_genes)

setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/Gencodev26_exons_plotting.rda")

BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
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

amy=BS[,which(pData(BS)$Tissue=="Brain - Amygdala")]
amy_VMRs=GRanges(list_of_filtered_autosomal_vmrs_wSDS[[6]])
converted_genes_TSS=resize(converted_genes,width=10000,fix="start")
pdf("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/AMY_VMRs_differential_genes.pdf")
plotManyRegions(amy,
                regions=converted_genes_TSS,
                extend = 1000,geneTrack=exon,
                addRegions = amy_VMRs)
dev.off()


#### Let's cluster samples based on VMRs within 1kb of these genes (remove chrX genes)
converted_genes2=converted_genes+1000
interest=subsetByOverlaps(amy_VMRs,converted_genes2)
ol=findOverlaps(amy_VMRs,converted_genes2)

meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
Individual=pData(amy)$Participant_ID
sex=pData(amy)$Self_Reported_Sex

ha = HeatmapAnnotation(type=Individual,sex=sex, 
    col = list(Individual = c("PT-1399T"="#e7298a","PT-13N2G"="#a6761d","PT-13O3Q"="#666666",
"PT-13O3O"="#e6ab02","PT-13OVH"="#1b9e77","PT-13NYS"="#66a61e","PT-13OW7"="#d95f02","PT-P44H"="#7570b3","PT-QDT8"="green","PT-X4XX"="black","PT-WZTO"="yellow"),
  sex=c("Male"="blue","Female"="red")))

gen_DMR_meth2=as.matrix(getMeth(amy,regions=amy_VMRs[queryHits(ol),],type="smooth",what="perRegion"))

pdf("AmyVMRs_olapDiffGenes_heatmap_dendrogram2.pdf")
Heatmap(as.matrix(gen_DMR_meth), name = "methylation", col = meth_col_fun, top_annotation = ha, 
    show_row_names = FALSE, km=4,show_column_names = FALSE,use_raster=TRUE,raster_device="png") 
dev.off()


clusters=hclust(dist(gen_DMR_meth2))

clusterCut=cutree(clusters,4)
table(clusterCut)

clusterCut2=as.data.frame(clusterCut) 

converted_genes2_ol=converted_genes2[subjectHits(ol),]$symbol

v=amy_VMRs[queryHits(ol),]
v$cluster=clusterCut2[,1]
v$gene=converted_genes2_ol

x=unique(v[which(v$cluster==2),]$gene)
y=unique(v[which(v$cluster==1),]$gene)
z=unique(v[which(v$cluster==3),]$gene)
w=unique(v[which(v$cluster==4),]$gene)


hmm=v[which(v$gene %in% rownames(t(t(t[which(t[,1]==1),])))),]
subsetByOverlaps(Amy_VMRs,v)
##Plot examples
pData(amy)$Participant_ID
col=c("#e7298a","#e7298a","#e6ab02","#a6761d","#a6761d","#a6761d","#a6761d","#e7298a","#e7298a","#e6ab02","#e6ab02")
#Group 2="#e7298a"
#Group 1="#a6761d"
#Group 3="#e6ab02"
CART=v[grep("CART",v$gene),]
GAB=v[grep("GAB",v$gene),]


pdf("AMY_VMRs_SLC17A7_CART_GAB_examples.pdf")
plotManyRegions(amy,col=col,
                regions=c(CART,SLC,GAB),
                extend = 5000,geneTrack=exon,
                addRegions = amy_VMRs)

dev.off()



