### Gene Annotations for Plotting 
library(rtracklayer)
# Get Gencode v26 gene models from GTEx
gv26_gtf=import("/dcl01/FB2/data/personal/gtex/dbgap/dbGaP-7882/59386/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf")
# keep protein coding and small RNA genes
gv26_gtf2=gv26_gtf[gv26_gtf$transcript_type %in% c("protein_coding","miRNA","lincRNA","snoRNA","snRNA","macro_lncRNA","scaRNA"),]
gv26_gtf2 <- keepSeqlevels(
  x = gv26_gtf2,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")
exons=gv26_gtf2[gv26_gtf2$type=="exon",]
num=strsplit(exons$exon_id, "[ ]")
num=sapply(num,"[[",3)
exons$exon_number=num
exon=as.data.frame(exons,stringsAsFactors=F)

colnames(exon)=c("chr","start","end","width","strand","source","type","score","phase","gene_id",
  "transcript_id","gene_type","gene_name","transcript_type","transcript_name","level",
  "havana_gene","exon_id","tag","exon_number")
exon$isoforms=1

save(exon,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/Gencodev26_exons_plotting.rda")


