## GREAT analysis of top 2000 HC-DMRs
library(rGREAT) # v4.0
dmrs=readRDS("/Users/lrizzardi/Desktop/Transcend/eGTEx/DMRs/objects/hippocampus_CG-DMRs.rds")
dmrs=dmrs[order(abs(dmrs$areaStat),decreasing=T),]
dmrs=GRanges(dmrs)
job = submitGreatJob(dmrs[1:2000,],species="hg38",adv_span=100,rule="basalPlusExt")
tb = getEnrichmentTables(job)

GO_BP=tb[[2]]
GO_BP=GO_BP[order(GO_BP$Hyper_Adjp_BH),]
GO_BP_sig=GO_BP[which(GO_BP$Hyper_Adjp_BH<0.01),]
GO_BP_sig=GO_BP_sig[which(GO_BP_sig$Binom_Adjp_BH<0.01),]
GO_BP_sig2=GO_BP_sig[order(GO_BP_sig$Hyper_Fold_Enrichment,decreasing=T),]
write.csv(GO_BP_sig2,file="/Users/lrizzardi/Dropbox/eGTEx/tables/HC_DMRs_2000_GREAT.csv")

### GREAT analysis of 2295 unique BG-DMRs
dmrs=readRDS("/Users/lrizzardi/Desktop/Transcend/eGTEx/DMRs/objects/basal_ganglia_CG-DMRs.rds")
dmrs=GRanges(dmrs)
all_dmrs=readRDS("/Users/lrizzardi/Desktop/Transcend/eGTEx/DMRs/objects/basal_ganglia_CG-DMRs.rds")
all_dmrs=GRanges(all_dmrs)
bg=subsetByOverlaps(dmrs,all_dmrs,invert=T)

job = submitGreatJob(bg,species="hg38",adv_span=100,rule="basalPlusExt")
tb = getEnrichmentTables(job)

GO_BP=tb[[2]]
GO_BP=GO_BP[order(GO_BP$Hyper_Adjp_BH),]
GO_BP_sig=GO_BP[which(GO_BP$Hyper_Adjp_BH<0.01),]
GO_BP_sig=GO_BP_sig[which(GO_BP_sig$Binom_Adjp_BH<0.01),]
GO_BP_sig2=GO_BP_sig[order(GO_BP_sig$Hyper_Fold_Enrichment,decreasing=T),]
write.csv(GO_BP_sig2,file="/Users/lrizzardi/Dropbox/eGTEx/tables/unique_BG_DMRs_GREAT.csv")
