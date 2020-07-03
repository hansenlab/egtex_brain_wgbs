
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
library(readr)
library(purrr)
library(UpSetR)

setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs_wSDS.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/Gencodev26_exons_plotting.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/snp_matrix.rda")
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

snp=snp_matrix$map$snp.names
a=strsplit(snp,"[_]")
snpspos=as.data.frame(cbind(snp,sapply(a,"[[",1),as.numeric(sapply(a,"[[",2))),stringsAsFactors=FALSE)
colnames(snpspos)=c("snp","chr","start")
snpspos$end=snpspos$start
snp=GRanges(snpspos)

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
HC1 <- c("SM-AFUP5", "SM-AFUL4", "SM-AFUNM", "SM-AFUOK", "SM-AFUKC",
         "SM-AFUKP", "SM-AFUNG", "SM-AFUOR", "SM-AFUJO", "SM-AFUJF",
         "SM-AFUJZ", "SM-AFUMU")
HC2 <- c("SM-AFUMH", "SM-AFUKV", "SM-AFUKJ", "SM-AFUK5",
         "SM-AFUMM", "SM-AFUNS")

VMRs=list_of_filtered_autosomal_vmrs_wSDS # w/o outliers already

BA9=GRanges(VMRs[[1]])
Hypo=GRanges(VMRs[[2]])
Put=GRanges(VMRs[[3]])
Cau=GRanges(VMRs[[4]])
NAcc=GRanges(VMRs[[5]])
Amy=GRanges(VMRs[[6]])
BA24=GRanges(VMRs[[7]])
HC=GRanges(VMRs[[8]])
Lung=GRanges(VMRs[[9]])
Thy=GRanges(VMRs[[10]])
HC1=GRanges(VMRs[[11]])
HC2=GRanges(VMRs[[12]])

VMR_list=GRangesList(BA9,BA24,Amy,Put,Cau,NAcc,Hypo,HC1,HC2,Lung,Thy)
names(VMR_list)=c("BA9","BA24","Amy","Put","Cau","NAcc","Hypo","HC1","HC2","Lung","Thy")

#Need lists...use dmr_labels objects from Annotation pipeline
union=GenomicRanges::reduce(c(Thy,Lung,NAcc,Put,Cau,BA24,Hypo,BA9,Amy,HC1,HC2)) 

ol=findOverlaps(HC,union)
union$Hippocampus=FALSE
union$Hippocampus[subjectHits(ol)] = TRUE

ol=findOverlaps(HC1,union)
union$Hippocampus_1=FALSE
union$Hippocampus_1[subjectHits(ol)] = TRUE

ol=findOverlaps(HC2,union)
union$Hippocampus_2=FALSE
union$Hippocampus_2[subjectHits(ol)] = TRUE

ol=findOverlaps(Amy,union)
union$Amygdala=FALSE
union$Amygdala[subjectHits(ol)] = TRUE 

ol=findOverlaps(Hypo,union)
union$Hypothalamus=FALSE
union$Hypothalamus[subjectHits(ol)] = TRUE 

ol=findOverlaps(NAcc,union)
union$`Nucleus accumbens`=FALSE
union$`Nucleus accumbens`[subjectHits(ol)] = TRUE 


ol=findOverlaps(Put,union)
union$Putamen=FALSE
union$Putamen[subjectHits(ol)] = TRUE 

ol=findOverlaps(BA9,union)
union$`Frontal cortex (BA9)`=FALSE
union$`Frontal cortex (BA9)`[subjectHits(ol)] = TRUE 

ol=findOverlaps(BA24,union)
union$`Anterior cingulate cortex (BA24)`=FALSE
union$`Anterior cingulate cortex (BA24)`[subjectHits(ol)] = TRUE 

ol=findOverlaps(Cau,union)
union$Caudate=FALSE
union$Caudate[subjectHits(ol)] = TRUE 

ol=findOverlaps(Thy,union)
union$Thyroid=FALSE
union$Thyroid[subjectHits(ol)] = TRUE 
ol=findOverlaps(Lung,union)
union$Lung=FALSE
union$Lung[subjectHits(ol)] = TRUE 
test=as.data.frame(mcols(union))

# First convert T/F matrix to 1/0 matrix
test=test + 0

test=test[,c(7,8,3,6,9,5,4,1,2,10,11)]
sets=colnames(test)

x=list(list("Amygdala"),list("Hippocampus_2"),list("Hippocampus_1"), list("Hypothalamus"),list("Nucleus.accumbens"),list("Putamen"),list("Caudate"),
 list("Frontal.cortex..BA9."),list("Anterior.cingulate.cortex..BA24."),list("Lung"),list("Thyroid"))
x2=list(list("Amygdala","Hippocampus_2")
,list("Amygdala","Hippocampus_1")
,list("Amygdala","Nucleus.accumbens")
,list("Amygdala","Hypothalamus")
,list("Amygdala","Hippocampus_1","Hippocampus_2")
,list("Putamen","Caudate")
,list("Hippocampus_1","Hippocampus_2")
,list("Amygdala","Hippocampus_2","Hypothalamus")
,list("Hippocampus_2","Hypothalamus")
,list("Nucleus.accumbens","Putamen","Caudate")
,list("Nucleus.accumbens","Amygdala","Hippocampus_2")
,list("Nucleus.accumbens","Hippocampus_2")
,list("Nucleus.accumbens","Hypothalamus")
,list("Nucleus.accumbens","Putamen")
,list("Nucleus.accumbens","Caudate")
,list("Nucleus.accumbens","Hypothalamus","Amygdala")
,list("Amygdala","Frontal.cortex..BA9.","Hippocampus_2")
,list("Amygdala","Frontal.cortex..BA9.")
,list("Hippocampus_1","Hippocampus_2","Amygdala","Hypothalamus","Nucleus.accumbens","Putamen","Frontal.cortex..BA9.","Anterior.cingulate.cortex..BA24.","Caudate","Thyroid","Lung")
,list("Hippocampus_1","Hippocampus_2","Amygdala","Hypothalamus","Nucleus.accumbens","Putamen","Frontal.cortex..BA9.","Anterior.cingulate.cortex..BA24.","Caudate")
,list("Thyroid","Lung")
,list("Amygdala","Anterior.cingulate.cortex..BA24.")
,list("Hippocampus_1","Hypothalamus")
,list("Amygdala","Putamen")
,list("Amygdala","Anterior.cingulate.cortex..BA24.","Hippocampus_2")
,list("Nucleus.accumbens","Hypothalamus","Amygdala","Hippocampus_2")
,list("Amygdala","Anterior.cingulate.cortex..BA24.","Hippocampus_2","Frontal.cortex..BA9.")
,list("Nucleus.accumbens","Hippocampus_1")
,list("Amygdala","Anterior.cingulate.cortex..BA24.","Hippocampus_2","Hippocampus_1","Frontal.cortex..BA9.")
,list("Nucleus.accumbens","Hippocampus_1","Amygdala")
,list("Hippocampus_2","Frontal.cortex..BA9.")
,list("Hippocampus_2","Putamen")
,list("Amygdala","Hypothalamus","Hippocampus_1"))


pdf("Full_VMR_upSet.pdf")
upset(test,number.angles = 90, sets.bar.color="deepskyblue",intersections = c(x,x2),order.by= "freq",keep.order=T,empty.intersection="off",sets.x.label = "VMRs")
dev.off()

## To get which VMRs in common between all groups just do this....

get_intersect_members <- function (x, comp){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(comp)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='numeric'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='numeric'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}


test3=test[,2:12]
y=get_intersect_members(test3, comp=colnames(test3))

# how many VMRs are common among just brain regions?
#common_VMRs=union[as.numeric(rownames(x)),] #333  

## how many VMRs are common to all tissues?
common_VMRs=union[as.numeric(rownames(y)),] #202

#tissue-specific VMRs

x=get_intersect_members(test2, comp="Hippocampus_1")
HC1_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Hippocampus_2")
HC2_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Amygdala")
Amy_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Hypothalamus")
Hypo_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Nucleus.accumbens")
NACC_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Putamen")
PUT_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Frontal.cortex..BA9.")
BA9_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Anterior.cingulate.cortex..BA24.")
BA24_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Caudate")
CAU_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Thyroid")
Thy_VMRs=union[as.numeric(rownames(x)),]
x=get_intersect_members(test2, comp="Lung")
Lung_VMRs=union[as.numeric(rownames(x)),]


## All VMRs for each tissue

NAcc_all=union[as.numeric(rownames(test2[which(test2$Nucleus.accumbens == 1),]))]
BA24_all=union[as.numeric(rownames(x[which(x$Anterior.cingulate.cortex..BA24.== 1),]))]
HC1_all=union[as.numeric(rownames(x[which(x$Hippocampus_1 == 1),]))]
HC2_all=union[as.numeric(rownames(x[which(x$Hippocampus_2 == 1),]))]
AMY_all=union[as.numeric(rownames(x[which(x$Amygdala == 1),]))]
HYP_all=union[as.numeric(rownames(x[which(x$Hypothalamus == 1),]))]
BA9_all=union[as.numeric(rownames(x[which(x$Frontal.cortex..BA9. == 1),]))]
CAU_all=union[as.numeric(rownames(x[which(x$Caudate == 1),]))]
PUT_all=union[as.numeric(rownames(x[which(x$Putamen == 1),]))]
Lung_all=union[as.numeric(rownames(x[which(x$Lung == 1),]))]
Thyroid_all=union[as.numeric(rownames(x[which(x$Thyroid == 1),]))]

# Now can plot these VMRs in each tissue to confirm shared ones and unique ones
pdf("Common_VMRs.pdf")
plotManyRegions(BS[,which(pData(BS)$Tissue=="Brain - Hippocampus")],
                common_VMRs[1:5, ],col="#a6cee3",
                extend = 5000,
                addRegions = union)

plotManyRegions(BS[,which(pData(BS)$Tissue=="Brain - Amygdala")],
                common_VMRs[1:5, ],col="#d95f02",
                extend = 5000,
                addRegions = union)

plotManyRegions(BS[,which(pData(BS)$Tissue=="Thyroid")],
                common_VMRs[1:5, ],col="#377eb8",
                extend = 5000,
                addRegions = union)

dev.off()

## Plot tissue-specific VMRs 
new_VMR_list=list(
"Lung_VMR"=Lung_VMRs,
"NAcc_all"=NAcc_all,
"BA24_all"=BA24_all,
"HC1_all"=HC1_all,
"HC2_all"=HC2_all,
"AMY_all"=AMY_all,
"HYP_all"=HYP_all,
"BA9_all"=BA9_all,
"CAU_all"=CAU_all,
"PUT_all"=PUT_all,
"Lung_all"=Lung_all,
"Thyroid_all"=Thyroid_all,
"Common_VMRs"=common_VMRs,
"HC1_only"=HC1_VMRs,
"HC2_only"=HC2_VMRs,
"Amy_only"=Amy_VMRs,
"Hypo_only"=Hypo_VMRs,
"NACC_only"=NACC_VMRs,
"PUT_only"=PUT_VMRs,
"BA9_only"=BA9_VMRs,
"BA24_only"=BA24_VMRs,
"CAU_only"=CAU_VMRs,
"Thy_only"=Thy_VMRs,
"Lung_only"=Lung_VMRs)

## Just do the ONLY category for now
 ba9=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[1]]),new_VMR_list[[20]])
ba24=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[7]]),new_VMR_list[[21]])
 cau=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[4]]),new_VMR_list[[22]])
 put=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[3]]),new_VMR_list[[19]])
 nac=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[5]]),new_VMR_list[[18]])
 hyp=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[2]]),new_VMR_list[[17]])
 hc1=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[11]]),new_VMR_list[[14]])
 hc2=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[12]]),new_VMR_list[[15]])
 amy=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[6]]),new_VMR_list[[16]])
 thy=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[10]]),new_VMR_list[[23]])
lung=subsetByOverlaps(GRanges(list_of_filtered_autosomal_vmrs_wSDS[[9]]),new_VMR_list[[24]])


new_VMR_list=list(amy,ba9,ba24,cau,put,nac,hyp,hc1,hc2,thy,lung)
names(new_VMR_list)=c("AMY","BA9","BA24","CAU","PUT","NAcc","HYP","HC1","HC2","THY","Lung")
new_VMR_list=lapply(new_VMR_list,function(x) x[order(-x$meanSDS),])

## Tissue-specific VMRs in Fig4B were assembled from these
pdf("Example_specific_VMRs_for_figs.pdf")
plotManyRegions(bsseq[,which(pData(bsseq)$Tissue=="Brain - Amygdala")],
                c(new_VMR_list[[1]][41,],new_VMR_list[[7]][41,],new_VMR_list[[4]][17,]),
                  col="#e7298a",main="Amy only VMRs",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = new_VMR_list[[1]])

plotManyRegions(bsseq[,which(pData(bsseq)$Tissue=="Brain - Hypothalamus")],
                c(new_VMR_list[[1]][41,],new_VMR_list[[7]][41,],new_VMR_list[[4]][17,]),
                col="#66a61e",main="HYP only VMRs",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = new_VMR_list[[7]])

plotManyRegions(bsseq[,which(pData(bsseq)$Tissue=="Brain - Caudate (basal ganglia)")],
                c(new_VMR_list[[1]][41,],new_VMR_list[[7]][41,],new_VMR_list[[4]][17,]),
                col="#666666",main="CAU only VMRs",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = new_VMR_list[[4]])
dev.off()

####
## Ubiquitous VMR Figure 4C
gtex =bsseq
gtex2 <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_non-brain_samples.mCG.small_smooth")


use=c("GTEX-12WSC","GTEX-12ZZY","GTEX-13FHP","GTEX-13FXS","GTEX-RU72","GTEX-X4XX","GTEX-ZE7O","GTEX-ZVT3","GTEX-13OW7","GTEX-13O3Q","GTEX-139TS","GTEX-13NZA")

same=pData(BS)$Collaborator_Participant_ID
gtex3=gtex2[,pData(gtex2)$Collaborator_Participant_ID %in% use]
gtex3=gtex3[,pData(gtex3)$Tissue =="Thyroid"]

gtex3 <- keepSeqlevels(
  x = gtex3,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")
gtex3=sort(sortSeqlevels(gtex3))


gtex3 <- gtex3[, !colnames(gtex3) %in%
                                   genotype_failures]

gtex=BS[,pData(BS)$Tissue %in% c("Brain - Frontal Cortex (BA9)","Brain - Caudate (basal ganglia)")]
gtex=gtex[,pData(gtex)$Collaborator_Participant_ID %in% use]

a=common_VMRs[seqnames(common_VMRs)=="chr19",]

pdf("Common_VMRs_for_figure4C.pdf")
plotManyRegions(gtex[,which(pData(gtex)$Tissue=="Brain - Frontal Cortex (BA9)")],
                a[3,],
                  col="#e6ab02",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = common_VMRs)

plotManyRegions(gtex3[,which(pData(gtex3)$Tissue=="Thyroid")],
                a[3,],
                col="#377eb8",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = common_VMRs)

plotManyRegions(gtex[,which(pData(gtex)$Tissue=="Brain - Caudate (basal ganglia)")],
                a[3,],
                col="#666666",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = common_VMRs)
dev.off()


### Are there VMRs in the highly variable MHC region? yes
MHC_1=GRanges(seqnames = "chr6",IRanges(start=28477797,end=33448354))
start(MHC_1)=start(MHC_1) - 100000
end(MHC_1)=end(MHC_1) + 100000

findOverlaps(union,MHC_1)  # 12 VMRs in MHC1 region



#### Get VMRs common to Hypo and AMY for Shared VMR fig 4B
y=get_intersect_members(test[,2:12], comp=c("Amygdala","Hypothalamus"))  # 949
amy_hypo=union[as.numeric(rownames(y)),]
library(rGREAT)  
job = submitGreatJob(amy_hypo,species="hg38",adv_span=100,rule="basalPlusExt")
tb = getEnrichmentTables(job)
names(tb)
GO_BP=tb[[2]]
GO_BP=GO_BP[order(GO_BP$Hyper_Adjp_BH),]
GO_BP_sig=GO_BP[which(GO_BP$Hyper_Adjp_BH<0.01),]
GO_BP_sig=GO_BP_sig[which(GO_BP_sig$Binom_Adjp_BH<0.01),]

GO_BP_sig2=GO_BP_sig[order(GO_BP_sig$Hyper_Fold_Enrichment,decreasing=T),]
res = plotRegionGeneAssociationGraphs(job, ontology = "GO Biological Process",
                                      termID = "GO:0006836")
#SLC32A1 and SLC6A1
slc2= amy_hypo[seqnames(amy_hypo)=="chr20"][19,]
slc1=amy_hypo[seqnames(amy_hypo)=="chr3"][5,]
slc=c(slc1,slc2)
pdf("SLC_VMRs.pdf")
plotManyRegions(BS[,which(pData(BS)$Tissue=="Brain - Caudate (basal ganglia)")],
                slc,col="#666666",
                extend = 5000,geneTrack=exon,annoTrack=list("SNP"=snp),
                addRegions = union)

plotManyRegions(BS[,which(pData(BS)$Tissue=="Brain - Amygdala")],
                slc,col="#d95f02",geneTrack=exon,annoTrack=list("SNP"=snp),
                extend = 5000,
                addRegions = union)

plotManyRegions(BS[,which(pData(BS)$Tissue=="Brain - Hypothalamus")],
                slc,col="#66a61e",geneTrack=exon,annoTrack=list("SNP"=snp),
                extend = 5000,
                addRegions = union)

dev.off()

### Get methylation for each tissue/individual over regions


#First strip list_of_VMRs of metadata to make results cleaner
test=list_of_VMRs
test2=unlist(test)
mcols(test2)=NULL
test3=relist(test2,test)
list_of_VMRs2=relist(test2,test)
list_of_VMRs2=list_of_VMRs2[grep("_all",names(list_of_VMRs2))]
list_of_VMRs2=list_of_VMRs2[1:9]
## Get Meth of all VMRs
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
HC1 <- c("SM-AFUP5", "SM-AFUL4", "SM-AFUNM", "SM-AFUOK", "SM-AFUKC",
         "SM-AFUKP", "SM-AFUNG", "SM-AFUOR", "SM-AFUJO", "SM-AFUJF",
         "SM-AFUJZ", "SM-AFUMU")
HC2 <- c("SM-AFUMH", "SM-AFUKV", "SM-AFUKJ", "SM-AFUK5",
         "SM-AFUMM", "SM-AFUNS")
## get meth only for the "all" categories....can subset out the only if needed
M_HC1=as.matrix(getMeth(BS[,which(pData(BS)$Sample_ID %in% HC1)],regions=list_of_VMRs2[[1]],type="smooth",what="perRegion"))
M_HC2=as.matrix(getMeth(BS[,which(pData(BS)$Sample_ID %in% HC2)],regions=list_of_VMRs2[[2]],type="smooth",what="perRegion"))
M_AMY=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Amygdala"],regions=list_of_VMRs2[[3]],type="smooth",what="perRegion"))
M_HYP=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Hypothalamus"],regions=list_of_VMRs2[[4]],type="smooth",what="perRegion"))
M_NACC=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Nucleus accumbens (basal ganglia)"],regions=list_of_VMRs2[[5]],type="smooth",what="perRegion"))
M_PUT=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Putamen (basal ganglia)"],regions=list_of_VMRs2[[6]],type="smooth",what="perRegion"))
M_BA9=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Frontal Cortex (BA9)"],regions=list_of_VMRs2[[7]],type="smooth",what="perRegion"))
M_BA24=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Anterior cingulate cortex (BA24)"],regions=list_of_VMRs2[[8]],type="smooth",what="perRegion"))
M_CAU=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Caudate (basal ganglia)"],regions=list_of_VMRs2[[9]],type="smooth",what="perRegion"))

M_BA24=M_BA24[,-grep("GTEX-1399T",colnames(M_BA24))]
M_BA24=M_BA24[,-grep("GTEX-13OVJ",colnames(M_BA24))]

M_BA9=M_BA9[,-grep("GTEX-1399T",colnames(M_BA9))]
M_BA9=M_BA9[,-grep("GTEX-13OVH",colnames(M_BA9))]
M_BA9=M_BA9[,-grep("GTEX-13FHP",colnames(M_BA9))]
M_BA9=M_BA9[,-grep("GTEX-RU72",colnames(M_BA9))]

M_list=list(M_HC1,M_HC2,M_AMY,M_HYP,M_NACC,M_PUT,M_BA9,M_BA24,M_CAU)
names(M_list)=c("HC1","HC2","AMY","HYP","NAcc","PUT","BA9","BA24","CAU")

M_list=lapply(M_list,function(x) {
  a=paste0("GTEX-",sapply(strsplit(as.character(colnames(x)),"[-]"),"[[",2))
  colnames(x)=a
  x})

save(M_list,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/VMR_all_brain_list_w_Meth.rda")

#### make karyotype plot in R
BiocManager::install("karyoploteR")
library(karyoploteR)

library(karyoploteR)
pdf("VMR_karyoplot.pdf")
kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, list_of_VMRs[[1]])
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.5,
                 minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")

dev.off()

pdf("CommonVMR_karyoplot.pdf")
kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, list_of_VMRs[[2]])
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.5,
                 minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")

dev.off()


#### VMR DMR overlap
hits <- findOverlaps(list_of_VMRs[[1]],dmrs,type="within")
overlaps <- pintersect(list_of_VMRs[[1]][queryHits(hits)], dmrs[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(list_of_VMRs[[1]][queryHits(hits)])
hits <- hits[percentOverlap > 0.5]


