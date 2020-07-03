##### Calculate Odds Ratios for CG- and CH-DMRs
library(gplots)
library(bsseq)
library(dplyr)
library(tidyr)
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/Reduced_CH-DMRs.rda")
x1=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/overall_CG-DMRs.rds"))
x2=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/overall_CG-blocks.rds"))
x3=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/general_CG-DMRs.rds"))
x4=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/hippocampus_CG-DMRs.rds"))
x5=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/basal_ganglia_CG-DMRs.rds"))
x6=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/general_CG-blocks.rds"))
x7=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/hippocampus_CG-blocks.rds"))
x8=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CG-DMRs/basal_ganglia_CG-blocks.rds"))
mcols(x1)=NULL
mcols(x2)=NULL
mcols(x3)=NULL
mcols(x4)=NULL
mcols(x5)=NULL
mcols(x6)=NULL
mcols(x7)=NULL
mcols(x8)=NULL

names=c("general_CH_DMRs", "basal_ganglia_CH_DMRs", "hippocampus_CH_DMRs","overall_CG-DMRs","overall_CG-blocks","general_CG-DMRs",
  "hippocampus_CG-DMRs","basal_ganglia_CG-DMRs","general_CG-blocks","hippocampus_CG-blocks","basal_ganglia_CG-blocks")
list_of_all_DMRs=GRangesList(general_CH_DMRs,basal_ganglia_CH_DMRs,hippocampus_CH_DMRs,x1,x2,x3,x4,x5,x6,x7,x8)
names(list_of_all_DMRs)=names

### Genomic Features

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/K27ac_enhancers_hg38.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/CGI-related-features-hg38.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/flattened-GENCODE-v26-features.rda")
K27ac_enhancers=GRangesList(K27ac_enhancers)
names(K27ac_enhancers)="K27ac_enhancers"

cgi_features=GRangesList("CGI"=cgis,"Shores"=shores,"Shelves"=shelves,"Open_Sea"=open_sea)

PsychENCODE_enh=import("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PsychEncode_enh_PFC/DER-04a_hg38lft_PEC_enhancers.bed")
PsychENCODE_highConf_enh=import("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PsychEncode_enh_PFC/DER-04b_hg38lft_high_confidence_PEC_enhancers.bed")
PsychENCODE_list=GRangesList(PsychENCODE_enh,PsychENCODE_highConf_enh)
names(PsychENCODE_list)=c("PsychENCODE_enh","PsychENCODE_highConf_enh")


features <- c(cgi_features,
        flattened_features_pc_transcripts[c("promoter", "five_utr", "exonic",
                                               "intronic", "three_utr",
                                               "intergenic")],K27ac_enhancers,PsychENCODE_list)

BS <- HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_non-brain_samples.mCG.small_smooth")
BS=keepStandardChromosomes(BS,pruning.mode="coarse")
BS=dropSeqlevels(BS,c("chrX","chrY","chrM"),pruning.mode="coarse")
BS=sort(sortSeqlevels(BS))
cpgs=granges(BS)

## original fullard data from the supplementary table S4

load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/Fullard_ATAC_hg38.rda")
mcols(Fullard_ATAC_hg38[[1]])=NULL
mcols(Fullard_ATAC_hg38[[2]])=NULL
mcols(Fullard_ATAC_hg38[[3]])=NULL
mcols(Fullard_ATAC_hg38[[4]])=NULL
mcols(Fullard_ATAC_hg38[[5]])=NULL
mcols(Fullard_ATAC_hg38[[6]])=NULL
mcols(Fullard_ATAC_hg38[[7]])=NULL
features3=Fullard_ATAC_hg38

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

DMRs_OR_df <- bind_rows(
        lapply(names(list_of_all_DMRs), function(name) {
                DMR <- list_of_all_DMRs[[name]]
                dplyr::bind_rows(FT(x = subsetByOverlaps(cpgs, DMR),
                             notx = subsetByOverlaps(cpgs, DMR, invert = TRUE),
                             features = features,
                             db = name))
        }))


mCG_OR_df <- bind_rows(DMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all <- mCG_OR_df %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]

new=mCG_OR_matrix_all[c(5,6,13,14,21,22),c(1,11:19,21:23)]



breaks=seq(-3, 3, length.out=101)  # set bc min is much greater than max
hmcol <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("CG-DMR_OR_newest.pdf")
heatmap.2(t(new),
          trace = "none",
          col = hmcol,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()

DMRs_OR_df <- bind_rows(
        lapply(names(list_of_all_DMRs), function(name) {
                DMR <- list_of_all_DMRs[[name]]
                dplyr::bind_rows(FT(x = subsetByOverlaps(cpgs, DMR),
                             notx = subsetByOverlaps(cpgs, DMR, invert = TRUE),
                             features = features3,
                             db = name))
        }))


mCG_OR_df <- bind_rows(DMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all <- mCG_OR_df %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]
breaks=seq(-2.5, 2.5, length.out=101)  # set bc min is much greater than max
hmcol <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("CG-DMR_OR_Fullard_newest.pdf")
heatmap.2(t(mCG_OR_matrix_all),
          trace = "none",
          col = hmcol,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()

#Enrichment of CH-DMRs
### NEED DIFF BACKGROUND WITH CH_DMRs
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/K27ac_enhancers_hg38.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/CGI-related-features-hg38.rda")
load("/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/flattened-GENCODE-v26-features.rda")
CG_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/general_CG-DMRs.rds")
BasGanglia_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/basal_ganglia_CG-DMRs.rds")
HC_DMRs=readRDS("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR/hippocampus_CG-DMRs.rds")
PsychENCODE_enh=import("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PsychEncode_enh_PFC/DER-04a_hg38lft_PEC_enhancers.bed")
PsychENCODE_highConf_enh=import("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PsychEncode_enh_PFC/DER-04b_hg38lft_high_confidence_PEC_enhancers.bed")
PsychENCODE_list=GRangesList(PsychENCODE_enh,PsychENCODE_highConf_enh)
names(PsychENCODE_list)=c("PsychENCODE_enh","PsychENCODE_highConf_enh")


CA_pos=GRanges(readRDS("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/CH-DMRs/general_CA_pos-DMRs.rds"))
CA_pos_max=CA_pos[which((CA_pos$maxDiff) >= 0.1),] #3195

plot_these=subsetByOverlaps(general_CH_DMRs,CA_pos_max) #3060




cgranges=granges(CA_pos_BSseq)
cgranges=keepStandardChromosomes(cgranges,pruning.mode="coarse")
cgranges=dropSeqlevels(cgranges,c("chrX","chrY","chrM"),pruning.mode="coarse")
names=c("general_CH_DMRs", "basal_ganglia_CH_DMRs", "hippocampus_CH_DMRs","TenPercent")
list_of_all_DMRs=GRangesList(general_CH_DMRs,basal_ganglia_CH_DMRs,hippocampus_CH_DMRs,plot_these)
names(list_of_all_DMRs)=names

features=c("K27ac"=K27ac_enhancers,"CG_DMRs"=CG_DMRs,"BG_DMRs"=BasGanglia_DMRs,"HC_DMRs"=HC_DMRs,"genic"=flattened_features_pc_transcripts$genic,
      "promoter" =flattened_features_pc_transcripts$promoter,"five_UTR"=flattened_features_pc_transcripts$five_utr,"three_UTR"=flattened_features_pc_transcripts$three_utr,
        "exonic"=flattened_features_pc_transcripts$exonic,"intronic"=flattened_features_pc_transcripts$intronic,"intergenic"=flattened_features_pc_transcripts$intergenic,
        "CGI"=cgis,"Shores"=shores,"Shelves"=shelves,"Open_Sea"=open_sea,"PsychENCODE_enh"=PsychENCODE_enh,"PsychENCODE_highConf_enh"=PsychENCODE_highConf_enh)
DMRs_OR_df <- bind_rows(
        lapply(names(list_of_all_DMRs), function(name) {
                DMR <- list_of_all_DMRs[[name]]
                dplyr::bind_rows(FT(x = subsetByOverlaps(cgranges, DMR),
                             notx = subsetByOverlaps(cgranges, DMR, invert = TRUE),
                             features = features,
                             db = name))
        }))


mCG_OR_df <- bind_rows(DMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all <- mCG_OR_df %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]

mCG_OR_matrix_all_ALLCHDMRs=mCG_OR_matrix_all
ALL_CHDMRs_OR_df=DMRs_OR_df


DMRs_OR_df <- dplyr::bind_rows(FT(x = subsetByOverlaps(cgranges, plot_these),
                             notx = subsetByOverlaps(cgranges, plot_these, invert = TRUE),
                             features = features,
                             db = "Ten % Diff CA-DMRs"))
        }))


mCG_OR_df <- bind_rows(DMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all <- mCG_OR_df %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]

x=rbind(mCG_OR_matrix_all_ALLCHDMRs,mCG_OR_matrix_all)

breaks=seq(-6, 2, length.out=101)  # set bc min is much greater than max
mid=which.min(abs(breaks-0))
hmcol <- c(colorRampPalette(c("blue", "white"))(mid),colorRampPalette(c("white", "red"))(100-mid))

pdf("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/Correct_CH_OR.pdf")
heatmap.2(t(x),
          trace = "none",
          col = hmcol,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()


### for atacseq fullard

DMRs_OR_df <- bind_rows(
        lapply(names(list_of_all_DMRs), function(name) {
                DMR <- list_of_all_DMRs[[name]]
                dplyr::bind_rows(FT(x = subsetByOverlaps(cgranges, DMR),
                             notx = subsetByOverlaps(cgranges, DMR, invert = TRUE),
                             features = Fullard_ATAC_hg38,
                             db = name))
        }))


mCG_OR_df <- bind_rows(DMRs_OR_df) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all <- mCG_OR_df %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all) <- mCG_OR_df[, 1]
breaks=seq(-1, 2, length.out=101)  # set bc min is much greater than max
#hmcol <- colorRampPalette(c("blue", "white", "red"))(100)
mid=which.min(abs(breaks-0))
hmcol <- c(colorRampPalette(c("blue", "white"))(mid),colorRampPalette(c("white", "red"))(100-mid))


pdf("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/Correct_CH_OR_Fullard.pdf")
heatmap.2(t(mCG_OR_matrix_all),
          trace = "none",
          col = hmcol,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()

#### ChromHMM enrichments

library(dplyr)
library(gplots)
library(circlize)
library(ggplot2)
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/chromHMM_brain_hg38.RData")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/VMR_OR_Enrich_results.RData")

files=list.files("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR",full.names = T)
names=strsplit
x1=GRanges(readRDS(files[1]))
x2=GRanges(readRDS(files[2]))
x3=GRanges(readRDS(files[3]))
x4=GRanges(readRDS(files[4]))
x5=GRanges(readRDS(files[5]))
x6=GRanges(readRDS(files[6]))
x7=GRanges(readRDS(files[7]))
x8=GRanges(readRDS(files[8]))
x9=GRanges(readRDS(files[9]))
x10=GRanges(readRDS(files[10]))
x11=GRanges(readRDS(files[11]))
x12=GRanges(readRDS(files[12]))
x13=GRanges(readRDS(files[13]))
x14=GRanges(readRDS(files[14]))
x15=GRanges(readRDS(files[15]))
x16=GRanges(readRDS(files[16]))
x17=GRanges(readRDS(files[17]))
x18=GRanges(readRDS(files[18]))
x19=GRanges(readRDS(files[19]))
x20=GRanges(readRDS(files[20]))
x21=GRanges(readRDS(files[21]))
x22=GRanges(readRDS(files[22]))
x23=GRanges(readRDS(files[23]))
x24=GRanges(readRDS(files[24]))

mcols(x1)=NULL
mcols(x2)=NULL
mcols(x3)=NULL
mcols(x4)=NULL
mcols(x5)=NULL
mcols(x6)=NULL
mcols(x7)=NULL
mcols(x8)=NULL
mcols(x9)=NULL
mcols(x10)=NULL
mcols(x11)=NULL
mcols(x12)=NULL
mcols(x13)=NULL
mcols(x14)=NULL
mcols(x15)=NULL
mcols(x16)=NULL
mcols(x17)=NULL
mcols(x18)=NULL
mcols(x19)=NULL
mcols(x20)=NULL
mcols(x21)=NULL
mcols(x22)=NULL
mcols(x23)=NULL
mcols(x24)=NULL

list_of_all_DMRs=GRangesList(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24)
files=list.files("/dcl01/FB2/data/personal/gtex/gtex/methylation/DMR")
names=strsplit(files,"[.]")
names=sapply(names,"[[",1)
names(list_of_all_DMRs)=names

#### ChromHMM Features
library(AnnotationHub)
ah <- AnnotationHub()

# Query AnnotationHub for all brain chromHMM tracks
queries <- query(ah, c("chromHMM", "brain"))   # 15-state model

# Keep the ones you do want
queries <- queries[c(1:3,5:8)]

#E071 Brain Hippocampus Middle    
#E074 Brain Substantia Nigra
#E068 Brain Anterior Caudate
#E069 Brain Cingulate Gyrus
#E072 Brain Inferior Temporal Lobe
#E067 Brain Angular Gyrus
#E073 Brain_Dorsolateral_Prefrontal_Cortex

                             
#  AH46920 | E067_15_coreMarks_mnemonics.bed.gz
#  AH46921 | E068_15_coreMarks_mnemonics.bed.gz
#  AH46922 | E069_15_coreMarks_mnemonics.bed.gz
#  AH46924 | E071_15_coreMarks_mnemonics.bed.gz
#  AH46925 | E072_15_coreMarks_mnemonics.bed.gz
#  AH46926 | E073_15_coreMarks_mnemonics.bed.gz
# Download data
chromHMM_tracks <- lapply(queries, function(x) {
  x[[1]]
})
AH46920=as.data.frame(chromHMM_tracks[[1]],stringsAsFactors=F)
AH46921=as.data.frame(chromHMM_tracks[[2]],stringsAsFactors=F)
AH46922=as.data.frame(chromHMM_tracks[[3]],stringsAsFactors=F)
AH46924=as.data.frame(chromHMM_tracks[[4]],stringsAsFactors=F)
AH46925=as.data.frame(chromHMM_tracks[[5]],stringsAsFactors=F)
AH46926=as.data.frame(chromHMM_tracks[[6]],stringsAsFactors=F)
AH46927=as.data.frame(chromHMM_tracks[[7]],stringsAsFactors=F)

ch = import.chain("/amber3/feinbergLab/personal/lrizzard/genomes/hg19ToHg38.over.chain")
x=AH46920 %>% split( .$name ) 
AH46920=lapply(x,function(x) GRanges(x))
names(AH46920)=paste0(names(x),"_E067")
AH46920=lapply(AH46920,function(x) unlist(liftOver(x,ch)))

x=AH46921 %>% split( .$name ) 
AH46921=lapply(x,function(x) GRanges(x))
names(AH46921)=paste0(names(x),"_E068")
AH46921=lapply(AH46921,function(x) unlist(liftOver(x,ch)))

x=AH46922 %>% split( .$name ) 
AH46922=lapply(x,function(x) GRanges(x))
names(AH46922)=paste0(names(x),"_E069")
AH46922=lapply(AH46922,function(x) unlist(liftOver(x,ch)))

x=AH46924 %>% split( .$name ) 
AH46924=lapply(x,function(x) GRanges(x))
names(AH46924)=paste0(names(x),"_E071")
AH46924=lapply(AH46924,function(x) unlist(liftOver(x,ch)))

x=AH46925 %>% split( .$name) 
AH46925=lapply(x,function(x) GRanges(x))
names(AH46925)=paste0(names(x),"_E072")
AH46925=lapply(AH46925,function(x) unlist(liftOver(x,ch)))

x=AH46926 %>% split( .$name ) 
AH46926=lapply(x,function(x) GRanges(x))
names(AH46926)=paste0(names(x),"_E073")
AH46926=lapply(AH46926,function(x) unlist(liftOver(x,ch)))

x=AH46927 %>% split( .$name ) 
AH46927=lapply(x,function(x) GRanges(x))
names(AH46927)=paste0(names(x),"_E074")
AH46927=lapply(AH46927,function(x) unlist(liftOver(x,ch)))

save(AH46920,AH46921,AH46922,AH46924,AH46925,AH46926,AH46927,file="chromHMM_brain_hg38.RData")

features_chromHMM <- c(AH46920,AH46921,AH46922,AH46924,AH46925,AH46926,AH46927)


DMRs_OR_df_chromHMM <- bind_rows(
        lapply(names(list_of_all_DMRs), function(name) {
                DMR <- list_of_all_DMRs[[name]]
                bind_rows(FT(x = subsetByOverlaps(cpgs, DMR),
                             notx = subsetByOverlaps(cpgs, DMR, invert = TRUE),
                             features = features_chromHMM,
                             db = name))
        }))

mCG_OR_df_chromHMM <- bind_rows(DMRs_OR_df_chromHMM) %>%
        mutate(`log2(OR)` = ifelse(is.na(`p.value`), NA, log2(estimate))) %>%
        dplyr::select(db, feature, `log2(OR)`) %>%
        spread(feature, `log2(OR)`)


mCG_OR_matrix_all_chromHMM <- mCG_OR_df_chromHMM %>% 
        dplyr::select(-db) %>%
        as.matrix()
rownames(mCG_OR_matrix_all_chromHMM) <- mCG_OR_df_chromHMM[, 1]


x=rownames(t(mCG_OR_matrix_all_chromHMM))
y=strsplit(x,"[_]")
cat=sapply(y,"[[",2)

summary=as.data.frame(t(mCG_OR_matrix_all_chromHMM),stringsAsFactors=F)
summary$group=cat
t=as_tibble(summary)
summary_table=t %>% group_by(group) %>% summarize_all(mean)
s=as.matrix(summary_table[,c(2:12)])
rownames(s)=summary_table$group
breaks=seq(-5, 4, length.out=101)
pdf("NEW_DMR_OR_new_chromHMM_summary.pdf")
heatmap.2(s,
          trace = "none",
          col = hmcol,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()

