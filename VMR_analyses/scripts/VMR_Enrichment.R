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
library(circlize)
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/")

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs.rda")

VMRs=list_of_filtered_autosomal_vmrs.nq

BA9=GRanges(VMRs[[1]]$no_outliers)
Hypo=GRanges(VMRs[[2]]$no_outliers)
Put=GRanges(VMRs[[3]]$no_outliers)
Cau=GRanges(VMRs[[4]]$no_outliers)
NAcc=GRanges(VMRs[[5]]$no_outliers)
Amy=GRanges(VMRs[[6]]$no_outliers)
BA24=GRanges(VMRs[[7]]$no_outliers)
Lung=GRanges(VMRs[[9]]$no_outliers)
Thy=GRanges(VMRs[[10]]$no_outliers)
HC1=GRanges(VMRs[[11]]$no_outliers)
HC2=GRanges(VMRs[[12]]$no_outliers)

union=GenomicRanges::reduce(c(Thy,Lung,NAcc,Put,Cau,BA24,Hypo,BA9,Amy,HC1,HC2)) 

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
test=test + 0
sets=colnames(test)
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


x=get_intersect_members(test, comp=colnames(test))
common_VMRs_all=union[as.numeric(rownames(x)),] #202
#unique VMRs

x=get_intersect_members(test, comp="Hippocampus_1")
HC1_VMRs=subsetByOverlaps(HC1,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Hippocampus_2")
HC2_VMRs=subsetByOverlaps(HC2,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Amygdala")
Amy_VMRs=subsetByOverlaps(Amy,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Hypothalamus")
Hypo_VMRs=subsetByOverlaps(Hypo,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Nucleus.accumbens")
NACC_VMRs=subsetByOverlaps(NAcc,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Putamen")
PUT_VMRs=subsetByOverlaps(Put,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Frontal.cortex..BA9.")
BA9_VMRs=subsetByOverlaps(BA9,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Anterior.cingulate.cortex..BA24.")
BA24_VMRs=subsetByOverlaps(BA24,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Caudate")
CAU_VMRs=subsetByOverlaps(Cau,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Thyroid")
Thy_VMRs=subsetByOverlaps(Thy,union[as.numeric(rownames(x)),])
x=get_intersect_members(test, comp="Lung")
Lung_VMRs=subsetByOverlaps(Lung,union[as.numeric(rownames(x)),])

# get VMRs with and without SNPS
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/snp_matrix.rda")
x=snp_matrix$map
a=strsplit(x[,1],"[_]")
snpspos_gr=as.data.frame(cbind(x[,1],sapply(a,"[[",1),as.numeric(sapply(a,"[[",2))),stringsAsFactors=FALSE)
colnames(snpspos_gr)=c("snp","chr","start")
snpspos_gr$end=snpspos_gr$start
snpspos_gr=GRanges(snpspos_gr)

mcols(common_VMRs_all)=NULL
mcols(union)=NULL

list_of_VMRs<-GRangesList("All_VMRs"=union,"common_VMRs"=common_VMRs_all,"HC1_all"=HC1,"HC1_only"=HC1_VMRs,"HC2_only"=HC2_VMRs,"HC2_all"=HC2,"AMY_only"=Amy_VMRs,"AMY_all"=Amy,
	"HYP_only"=Hypo_VMRs,"HYP_all"=Hypo,"NAcc_only"=NACC_VMRs,"NAcc_all"=NAcc,"PUT_only"=PUT_VMRs,"PUT_all"=Put,"BA9_only"=BA9_VMRs,"BA9_all"=BA9,
"BA24_only"=BA24_VMRs,"BA24_all"=BA24,"CAU_all"=Cau,"CAU_only"=CAU_VMRs,"Thyroid_only"=Thy_VMRs,"Thyroid_all"=Thy,"Lung_only"=Lung_VMRs,"Lung_all"=Lung)
list_of_VMRs_w_snps=lapply(list_of_VMRs,function(x) subsetByOverlaps(x,snpspos_gr))
list_of_VMRs_wo_snps=lapply(list_of_VMRs,function(x) subsetByOverlaps(x,snpspos_gr,invert=T))
save(list_of_VMRs,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/List_of_VMRs.rda")

# How many VMRs overlaps DMRs
olap=lapply(list_of_DMRs,function(x){subsetByOverlaps(list_of_VMRs[[1]],x)})


#### Enrichment Functions

FT <- function(x, notx, features, db, offset = 1L) {
  x <- lapply(features, function(feature) {
    if (length(feature)) {
      n11 <- sum(overlapsAny(x, feature))
      n12 <- length(x) - n11
      n21 <- sum(overlapsAny(notx, feature))
      n22 <- length(notx) - n21
      m <- matrix(c(n11, n21, n12, n22), ncol = 2) + offset
      ft <- fisher.test(m)
    } else {
      list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_)
    }
  })

  data.frame(db = db,
             feature = names(features),
             estimate = sapply(x, "[[", "estimate"),
             lower = sapply(lapply(x, "[[", "conf.int"), "[[", 1),
             upper = sapply(lapply(x, "[[", "conf.int"), "[[", 2),
             p.value = sapply(x, "[[", "p.value"),
             stringsAsFactors = FALSE,
             row.names = NULL)
}

##### Load Features #####
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_non-brain_samples.mCG.small_smooth")
BS=keepStandardChromosomes(BS,pruning.mode="coarse")
BS=dropSeqlevels(BS,c("chrX","chrY","chrM"),pruning.mode="coarse")
BS=sort(sortSeqlevels(BS))
cpgs=granges(BS)

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/K27ac_enhancers_hg38.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/CGI-related-features-hg38.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/flattened-GENCODE-v26-features.rda")
CG_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/general_CG-DMRs.rds")
BasGanglia_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/basal_ganglia_CG-DMRs.rds")
HC_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/hippocampus_CG-DMRs.rds")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/Reduced_CH-DMRs.rda")

CG_DMRs=GRanges(CG_DMRs)
mcols(CG_DMRs)=NULL
BasGanglia_DMRs=GRanges(BasGanglia_DMRs)
mcols(BasGanglia_DMRs)=NULL
HC_DMRs=GRanges(HC_DMRs)
mcols(HC_DMRs)=NULL

list_of_DMRs=GRangesList("CG-DMRs"=CG_DMRs,"ganglia_DMRs"=BasGanglia_DMRs,"HC_DMRs"=HC_DMRs,
  "CH-DMRs"=general_CH_DMRs)
cgi_features=GRangesList("CGI"=cgis,"Shores"=shores,"Shelves"=shelves,"Open_Sea"=open_sea)

K27ac_enhancers=GRangesList(K27ac_enhancers)
names(K27ac_enhancers)="K27ac_enhancers"

PsychENCODE_enh=import("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PsychEncode_enh_PFC/DER-04a_hg38lft_PEC_enhancers.bed")
PsychENCODE_highConf_enh=import("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PsychEncode_enh_PFC/DER-04b_hg38lft_high_confidence_PEC_enhancers.bed")
PsychENCODE_list=GRangesList(PsychENCODE_enh,PsychENCODE_highConf_enh)
names(PsychENCODE_list)=c("PsychENCODE_enh","PsychENCODE_highConf_enh")


features <- c(cgi_features,
        flattened_features_pc_transcripts[c("promoter", "five_utr", "exonic",
                                               "intronic", "three_utr",
                                               "intergenic")],
        list_of_DMRs,K27ac_enhancers,PsychENCODE_list)


VMRs_OR_df <- bind_rows(
        lapply(names(list_of_VMRs), function(name) {
                VMR <- list_of_VMRs[[name]]
                bind_rows(FT(x = subsetByOverlaps(cpgs, VMR),
                             notx = subsetByOverlaps(cpgs, VMR, invert = TRUE),
                             features = features,
                             db = name))
        }))

mCG_OR_df <- bind_rows(VMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all <- mCG_OR_df %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]


hmcol <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("VMR_OR_new.pdf")
heatmap.2(t(mCG_OR_matrix_all),
          trace = "none",
          col = hmcol,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()

### Chrom HMM features 
load("chromHMM_brain_hg38.RData") #object created in DMR_enrichments.R script

features_chromHMM <- c(AH46920,AH46921,AH46922,AH46924,AH46925,AH46926,AH46927)

VMRs_OR_df_chromHMM <- bind_rows(
        lapply(names(list_of_VMRs), function(name) {
                VMR <- list_of_VMRs[[name]]
                bind_rows(FT(x = subsetByOverlaps(cpgs, VMR),
                             notx = subsetByOverlaps(cpgs, VMR, invert = TRUE),
                             features = features_chromHMM,
                             db = name))
        }))

mCG_OR_df_chromHMM <- bind_rows(VMRs_OR_df_chromHMM) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all_chromHMM <- mCG_OR_df_chromHMM %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all_chromHMM) <- mCG_OR_df_chromHMM[, 1]

#### Summarize across ChromHMM tissues
x=rownames(t(VMR_mCG_OR_matrix_all_chromHMM))
y=strsplit(x,"[_]")
cat=sapply(y,"[[",2)

summary=as.data.frame(t(VMR_mCG_OR_matrix_all_chromHMM),stringsAsFactors=F)
summary$group=cat
t=as_tibble(summary)
summary_table=t %>% group_by(group) %>% summarize_all(mean)
s=as.matrix(summary_table[,c(2:25)])
rownames(s)=summary_table$group
breaks=seq(-5, 4, length.out=101)
pdf("NEW_VMR_OR_chromHMM_summary.pdf")
heatmap.2(s,
          trace = "none",
          col = hmcol,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()
save(VMR_mCG_OR_matrix_all_chromHMM,features,features_chromHMM,VMR_mCG_OR_matrix_all,file="VMR_OR_Enrich_results.RData")


#########################################
### Check VMRs against those from PMID: 30273333 Supp Table 1 just looking at neuronal

pub=read.xlsx("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/TS1_PMID30273333_VMRs.xlsx",sheet=5,startRow=2)
pub=pub[,c(2:4)]
colnames(pub)=c("chr","start","end")
pub=GRanges(pub) #7076
# convert to hg38
ch = import.chain("/amber3/feinbergLab/personal/lrizzard/genomes/hg19ToHg38.over.chain")
pub=unlist(liftOver(pub,ch))
pub=pub+500
pub=GenomicRanges::reduce(pub)

y=lapply(list_of_VMRs,function(x) subsetByOverlaps(pub,x))
t(as.data.frame(lapply(y,length)))

### check CoRSIV from PMID:31155008 and subtelomeric enrichment......1Mb from ends of chromosomes

#### look at overlap with waterland CoRSIVs
paper=read.xlsx("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/13059_2019_1708_MOESM1_ESM.xlsx",3)
CoRSIV=paper$USCS_Coordinates_CoRSIV
chr=sapply(strsplit(CoRSIV,"[:]"),"[[",1)
start=sapply(strsplit(CoRSIV,"[:]"),"[[",2)
end=sapply(strsplit(start,"[-]"),"[[",2)
start=sapply(strsplit(start,"[-]"),"[[",1)
x=as.data.frame(cbind(chr,as.numeric(start),as.numeric(end)))
colnames(x)=c("chr","start","end")
CoRSIV=unique(GRanges(x))
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/List_of_VMRs.rda") 

features=list("CoRSIV"=CoRSIV,"nVMRs"=pub)
names=names(list_of_VMRs)
VMRs_OR_df <- bind_rows(
        lapply(names(list_of_VMRs), function(name) {
                VMRs <- list_of_VMRs[[name]]
                bind_rows(FT(x = subsetByOverlaps(cpgs, VMRs),
                             notx = subsetByOverlaps(cpgs, VMRs, invert = TRUE),
                             features = features,
                             db = name))
        }))

mCG_OR_df <- bind_rows(VMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)

mCG_OR_matrix_all <- mCG_OR_df %>%
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]

breaks=seq(-2,7, length.out=101)  # set bc min is much greater than max
mid=which.min(abs(breaks-0))
hmcol <- c(colorRampPalette(c("blue", "white"))(mid),colorRampPalette(c("white", "red"))(100-mid))

pdf("VMR_enrich_CoRSIV_nVMRs.pdf")
heatmap.2(t(mCG_OR_matrix_all),
          trace = "none",
          col = hmcol,breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()





