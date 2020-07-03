#####################
## meQTL Analysis ###
##                ###
#####################

#### Run all the LD stuff on the main GTEx dataset then filtering individuals
#### allowed missingness 0.1 rather than 0 bc so many samples in dataset and don't do maf filter until after subset individuals
### Outside of R:
/amber3/feinbergLab/personal/lrizzard/bin/plink2 --make-bpgen --autosome --geno 0.1 --export vcf --indep-pairwise 100'kb' 1 0.2 --vcf \
/dcl01/FB2/data/personal/gtex/dbgap/dbGaP-7882/59386/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz \
--out /dcl01/FB2/data/personal/gtex/dbgap/dbGaP-7882/59386/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_v8_MAF_phased_LD
### kept 699654 SNPs ih the GTEx_v8_MAF_phased_LD.prune.in file   use this to subset snp_matrix
## now subset the new VCF by samples
/amber3/feinbergLab/personal/lrizzard/bin/plink2 --keep /dcl01/FB2/data/personal/gtex/gtex/Individuals_Sorted_GTEx.txt --maf 0.1 --export vcf --vcf GTEx_v8_MAF_phased_LD.vcf --out ./GTEx_v8_MAF_phased_LD_Sorted_Ind_Only
# 4914120 SNPs letf


# Create snp matrix from the new vcf
library(VariantAnnotation)
library(bsseq)
library(HDF5Array)
library(SummarizedExperiment)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

BS=HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/QC_failures.rda")
bsseq_filtered <- BS[, !colnames(BS) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]
BS <- bsseq_filtered
BS=BS[,!pData(BS)$Participant_ID %in% c("PT-13NYS","PT-12WSC")]

samples=sampleNames(BS)
samples=paste0("GTEX-",sapply(strsplit(samples,"[-]"),"[[",2))
samples=unique(samples)

file="/dcl01/FB2/data/personal/gtex/dbgap/dbGaP-7882/59386/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_v8_MAF_phased_LD_Sorted_Ind_Only.vcf"
param = ScanVcfParam(samples=samples)
vcf=readVcf(file,genome="hg38",param=param)
snp_matrix=genotypeToSnpMatrix(vcf)

## Filter LD based on GTEx_v8_MAF_phased_LD.prune.in file
keep=readLines("/dcl01/FB2/data/personal/gtex/dbgap/dbGaP-7882/59386/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_v8_MAF_phased_LD.prune.in")
snp_matrix2=snp_matrix
snp_matrix2$genotypes=snp_matrix2$genotypes[,which(colnames(snp_matrix2$genotypes) %in% keep)]
snp_matrix2$map=snp_matrix2$map[which(snp_matrix2$map$snp.names %in% keep),]
save(snp_matrix2, file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/snpMatrix_VCF_LDadj_200124.rda") 
### leaves us with 149,185 SNPs to work with


## create data for running MatrixEQTL
library(MatrixEQTL)
write.table(t(snp_matrix2$genotypes),"LD_adjusted_snps_for_MatrixEQTL_2020.txt",sep="\t")
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/list_of_filtered_autosomal_vmrs_wSDS.rda")
vmrs=lapply(list_of_filtered_autosomal_vmrs_wSDS,GRanges)
# read in methylation data (this is just from Caudate)
 M_HC1=as.matrix(getMeth(BS[,pData(BS)$Sample_ID %in% HC1],regions=vmrs[[11]],type="smooth",what="perRegion"))
 M_HC2=as.matrix(getMeth(BS[,pData(BS)$Sample_ID %in% HC2],regions=vmrs[[12]],type="smooth",what="perRegion"))
 M_CAU=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Caudate (basal ganglia)"],regions=vmrs[[4]],type="smooth",what="perRegion"))
 M_AMY=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Amygdala"],regions=vmrs[[6]],type="smooth",what="perRegion"))
 M_HYP=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Hypothalamus"],regions=vmrs[[2]],type="smooth",what="perRegion"))
M_NACC=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Nucleus accumbens (basal ganglia)"],regions=vmrs[[5]],type="smooth",what="perRegion"))
 M_PUT=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Putamen (basal ganglia)"],regions=vmrs[[3]],type="smooth",what="perRegion"))
 M_BA9=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Frontal Cortex (BA9)"],regions=vmrs[[1]],type="smooth",what="perRegion"))
M_BA24=as.matrix(getMeth(BS[,BS$Tissue=="Brain - Anterior cingulate cortex (BA24)"],regions=vmrs[[7]],type="smooth",what="perRegion"))

  
vmrs=vmrs[c(1:7,11:12)]
M_list=list(M_BA9,M_HYP, M_PUT,M_CAU, M_NACC,M_AMY,M_BA24,M_HC1,M_HC2)
names(M_list)=names(vmrs)
M_list2=M_list

for (i in 1:9){
	a=colnames(M_list[[i]])
	names=strsplit(a,"[-]")
names3=sapply(names,"[[",2)
colnames(M_list[[i]])=paste0("GTEX-",names3)
loc=as.data.frame(vmrs[i])
rownames(M_list[[i]])=paste(loc[,1],loc[,2],sep="_")
write.table(M_list[[i]],file=paste0(names(M_list[i]),"_VMR_meanMeth.txt"),col.names=T,row.names=T,quote=F,sep="\t")
}
save(M_list,file="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/M_list_for_meQTL.rda")
## now need to order colnames 
#so they match up

all_names=rownames(snp_matrix2$genotypes)
meth_names=lapply(M_list,function(x){which(all_names %in% colnames(x))})
# reorder meth_names to match order of files (below)
meth_names=meth_names[c(6,7,4,1,2,5,3,8,9)]
files=list.files(pattern="meanMeth.txt",full.names=T)

a=strsplit(snp_matrix2$map$snp.names,"[_]")
snpspos=as.data.frame(cbind(snp_matrix2$map$snp.names,sapply(a,"[[",1),as.numeric(sapply(a,"[[",2))),stringsAsFactors=FALSE)
colnames(snpspos)=c("snp","chr","start")
snpspos$start=as.numeric(snpspos$start)
SNP_file_name="/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/VMR/LD_adjusted_snps_for_MatrixEQTL_2020.txt"


####################
## RUN MatrixEQTL ##
####################
me=list(1,2,3,4,5,6,7,8,9)
useModel=modelLINEAR
cisDist = 5e4 
pvOutputThreshold_cis=2e-05
errorCovariance = numeric()

## Below will be in for loop:
for (i in 1:9){
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$LoadFile( files[i])

a=strsplit(rownames(gene),"[_]")

genepos = as.data.frame(cbind(rownames(gene),sapply(a,"[[",1),as.numeric(sapply(a,"[[",2)),as.numeric(sapply(a,"[[",2))),stringsAsFactors=FALSE)
colnames(genepos)=c("id","chr","s1","s2")
genepos$s1=as.numeric(genepos$s1)
genepos$s2=as.numeric(genepos$s2)
x=meth_names[[i]]
snps = SlicedData$new()
snps$LoadFile( SNP_file_name )
snps$ColumnSubsample(x);

me2 = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name = NULL,
pvOutputThreshold = 0,
pvOutputThreshold.cis = pvOutputThreshold_cis,
output_file_name.cis=NULL,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE)
me[i]=list(me2)

}
## Determine appropriate pvalue cutoffs based on number of tests run:
pvalues=unlist(lapply(me,function(x){0.05/t(x$cis$ntests)})) 

### now rerun this with the appropriate pvalue cutoffs (could have just subset them, but some were bigger than the current reported cutoff.....)
me=list(1,2,3,4,5,6,7,8,9)
for (i in 1:9){
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$LoadFile( files[i])

a=strsplit(rownames(gene),"[_]")

genepos = as.data.frame(cbind(rownames(gene),sapply(a,"[[",1),as.numeric(sapply(a,"[[",2)),as.numeric(sapply(a,"[[",2))),stringsAsFactors=FALSE)
colnames(genepos)=c("id","chr","s1","s2")
genepos$s1=as.numeric(genepos$s1)
genepos$s2=as.numeric(genepos$s2)
x=meth_names[[i]]
snps$LoadFile( SNP_file_name )
snps$ColumnSubsample(x);

me2 = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name = NULL,
pvOutputThreshold = 0,
pvOutputThreshold.cis = pvalues[i],
output_file_name.cis=NULL,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE)
me[i]=list(me2)
}

# how many meQTLs identified:
names(me)=names(meth_names)
t(t(unlist(lapply(me,function(x){x$cis$neqtls})) ))
#                                         cisDist
#                                           50kb  totalVMRs  % w meQTL
# Brain - Amygdala                            6		34876		0.017%
# Brain - Anterior cingulate cortex (BA24)   21		7327		0.287%
# Brain - Caudate (basal ganglia)            96		8137		1.18
# Brain - Frontal Cortex (BA9)               60		8719		0.688
# Brain - Hypothalamus                        9		14602		0.062
# Brain - Nucleus accumbens (basal ganglia)  37		17144		0.216
# Brain - Putamen (basal ganglia)            42		10654		0.4
# HC1                                        17		17731		0.096
# HC2                                         0
results=lapply(me,function(x){x$cis$eqtls})
results_10kb=results
me_10kb=me
results_500kb=results
me_500kb=me
results_50kb=results
me_50kb=me
###SAVE TO TABLES FOR PUBLICATION 
library(openxlsx)
library(tidyr)
#make short names
short=c("AMY","BA24","Caudate","BA9","HYPO","NAcc","Putamen","HC1","HC2")
wb <- createWorkbook()
Map(function(data, name){
 
    addWorksheet(wb, name)
    writeData(wb, name, data)
 
}, results, short)

saveWorkbook(wb, file = "meQTL_results_50kb_200124.xlsx", overwrite = TRUE)

pdf("QQplot_test2.pdf")
lapply(me,function(x){plot(x, pch = 16, cex = 0.7)})
dev.off()

## Make one results GRanges not list
res2=results[1:8]

for (i in 1:8){
		res2[[i]]$Tissue=rep(short[i],times=nrow(results[[i]]))

}

res=rbind(res2[[1]],res2[[2]],res2[[3]],res2[[4]],res2[[5]],res2[[6]],res2[[7]],res2[[8]])
## make into GRanges based on VMR and SNPs
SNP_res=res
SNP_res$chr=sapply(strsplit(as.character(SNP_res$snps),"[_]"),"[[",1)
SNP_res$start=sapply(strsplit(as.character(SNP_res$snps),"[_]"),"[[",2)
SNP_res$end=SNP_res$start
SNP_res=GRanges(SNP_res)## 301 ranges


VMR_res=res
VMR_res$chr=sapply(strsplit(as.character(VMR_res$gene),"[_]"),"[[",1)
VMR_res$start=sapply(strsplit(as.character(VMR_res$gene),"[_]"),"[[",2)
VMR_res$end=VMR_res$start
VMR_res=GRanges(VMR_res)## 301 ranges
VMR_res_dup=VMR_res[duplicated(VMR_res),] # Which VMRs have multiple snps
dups=as.data.frame(VMR_res_dup[order(VMR_res_dup),],stringsAsFactors=F)[,c(6:7,12)] ## some are duplicated bc same VMR-SNP in multiple tissues, some multiple SNPs per tissue (77 total)
dups$gene=as.character(dups$gene)
t(t(table(dups$gene)))

# so there are 2 VMRs with SNPs in 6/8 tissues tested
# chr10_26059678
# chr4_32198541
# 
# # shared SNPs
# chr10_26066170_A_G_b38 #in 7/8 tissues
# chr4_32193194_A_G_b38. # in 6/8 tissues

### plot karyoplot
SNP_res_df=as.data.frame(SNP_res,stringsAsFactors=F)
SNP_res$color="#e6ab02"
SNP_res$color[which(SNP_res$Tissue=="HYPO")]="#66a61e"
SNP_res$color[which(SNP_res$Tissue=="Putamen")]="#7570b3"
SNP_res$color[which(SNP_res$Tissue=="Caudate")]="#666666"
SNP_res$color[which(SNP_res$Tissue=="NAcc")]="#d95f02"
SNP_res$color[which(SNP_res$Tissue=="BA24")]="#a6761d"
SNP_res$color[which(SNP_res$Tissue=="BA9")]="#e6ab02"
SNP_res$color[which(SNP_res$Tissue=="HC1")]="#1b9e77"
SNP_res$color[which(SNP_res$Tissue=="AMY")]="#e7298a"

VMR_res$color="#e6ab02"
VMR_res$color[which(VMR_res$Tissue=="HYPO")]="#66a61e"
VMR_res$color[which(VMR_res$Tissue=="Putamen")]="#7570b3"
VMR_res$color[which(VMR_res$Tissue=="Caudate")]="#666666"
VMR_res$color[which(VMR_res$Tissue=="NAcc")]="#d95f02"
VMR_res$color[which(VMR_res$Tissue=="BA24")]="#a6761d"
VMR_res$color[which(VMR_res$Tissue=="BA9")]="#e6ab02"
VMR_res$color[which(VMR_res$Tissue=="HC1")]="#1b9e77"
VMR_res$color[which(VMR_res$Tissue=="AMY")]="#e7298a"


pdf("mQTL_karyoplot.pdf")
kp <- plotKaryotype(genome="hg38",plot.type=2)
kpPlotRegions(kp, VMR_res,data.panel=2,col=VMR_res$color)
kpAxis(kp,data.panel=1,ymin=0,ymax=1)
kpPoints(kp, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)

kp <- plotKaryotype(genome="hg38",plot.type=2,chromosomes=c("chr1","chr2","chr3","chr4"))
kpPlotRegions(kp, VMR_res,data.panel=2,col=VMR_res$color)
kpAxis(kp,data.panel=1,ymin=0,ymax=1)
kpPoints(kp, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)

kp2 <- plotKaryotype(genome="hg38",plot.type=2,chromosomes=c("chr5","chr6","chr7","chr8"))
kpPlotRegions(kp2, VMR_res,data.panel=2,col=VMR_res$color)
kpAxis(kp2,data.panel=1,ymin=0,ymax=1)
kpPoints(kp2, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)

kp3 <- plotKaryotype(genome="hg38",plot.type=2,chromosomes=c("chr9","chr10","chr11","chr12"))
kpPlotRegions(kp3, VMR_res,data.panel=2,col=VMR_res$color)
kpAxis(kp3,data.panel=1,ymin=0,ymax=1)
kpPoints(kp3, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)

kp4 <- plotKaryotype(genome="hg38",plot.type=2,chromosomes=c("chr13","chr14","chr15","chr16"))
kpPlotRegions(kp4, VMR_res,data.panel=2,col=VMR_res$color)
kpAxis(kp4,data.panel=1,ymin=0,ymax=1)
kpPoints(kp4, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)

kp5 <- plotKaryotype(genome="hg38",plot.type=2,chromosomes=c("chr17","chr18","chr19","chr20","chr21","chr22"))
kpPlotRegions(kp5, VMR_res,data.panel=2,col=VMR_res$color)
kpAxis(kp5,data.panel=1,ymin=0,ymax=1)
kpPoints(kp5, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)

dev.off()


pdf("mQTL_karyoplot_allchr.pdf")
kp <- plotKaryotype(genome="hg38",plot.type=2,chromosomes=c("autosomal"))
kpPlotRegions(kp, VMR_res,data.panel=2,col=VMR_res$color,r0=0,r1=1)
kpPoints(kp, SNP_res, y=SNP_res$log10_FDR/10, col=SNP_res$color,r0=0,r1=1,data.panel=1)
dev.off()

# plot these shared ones
z=VMR_res_dup[which(VMR_res_dup$gene %in% c("chr10_26059678","chr4_32198541")),][1:2,]


BS_HC1=BS[,which(pData(BS)$Sample_ID %in% HC1)]
try=as.data.frame(snp_matrix2$genotypes[,"chr10_26066170_A_G_b38"])
a=strsplit(rownames(try),"[-]")
rownames(try)=paste0("PT-",sapply(a,"[[",2))
a=match(pData(BS_HC1)$Participant_ID,rownames(try))
b=try[a,]
pData(BS_HC1)$GT=b
pData(BS_HC1)$GT_col="red"
pData(BS_HC1)$GT_col=replace(pData(BS_HC1)$GT_col,pData(BS_HC1)$GT==01,"black")
pData(BS_HC1)$GT_col=replace(pData(BS_HC1)$GT_col,pData(BS_HC1)$GT==02,"blue")
#pData(BS_HC1)$GT_col=replace(pData(BS_HC1)$GT_col,pData(BS_HC1)$GT==00,"white") #these individuals dont have GT info

BS_new=BS[,which(pData(BS)$Tissue %in% c("Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)","Brain - Frontal Cortex (BA9)","Brain - Putamen (basal ganglia)"))]
try=as.data.frame(snp_matrix2$genotypes[,"chr10_26066170_A_G_b38"])
a=strsplit(rownames(try),"[-]")
rownames(try)=paste0("PT-",sapply(a,"[[",2))
a=match(pData(BS_new)$Participant_ID,rownames(try))
b=try[a,]
pData(BS_new)$GT=b
pData(BS_new)$GT_col="red"
pData(BS_new)$GT_col=replace(pData(BS_new)$GT_col,pData(BS_new)$GT==01,"black")
pData(BS_new)$GT_col=replace(pData(BS_new)$GT_col,pData(BS_new)$GT==02,"blue")

x=lapply(vmrs,function(x){findOverlaps(x,z[1,])})
y=lapply(vmrs,function(x){findOverlaps(x,z[2,])})
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/Gencodev26_exons_plotting.rda")

### plot hook2 mQTL  in BA9  BA24 Cau PUt
BS_new2=BS[,which(pData(BS)$Tissue %in% c("Brain - Anterior cingulate cortex (BA24)","Brain - Caudate (basal ganglia)","Brain - Frontal Cortex (BA9)","Brain - Putamen (basal ganglia)"))]
try2=as.data.frame(snp_matrix2$genotypes[,"chr19_12766150_G_C_b38"])
a=strsplit(rownames(try2),"[-]")
rownames(try2)=paste0("PT-",sapply(a,"[[",2))
a=match(pData(BS_new2)$Participant_ID,rownames(try2))
b=try2[a,]
pData(BS_new2)$GT=b
pData(BS_new2)$GT_col="red"
pData(BS_new2)$GT_col=replace(pData(BS_new2)$GT_col,pData(BS_new2)$GT==01,"black")
pData(BS_new2)$GT_col=replace(pData(BS_new2)$GT_col,pData(BS_new2)$GT==02,"blue")

z=VMR_res[which(VMR_res$snps %in% c("chr19_12766150_G_C_b38")),][1,]

x=lapply(list_of_VMRs[[1]],function(x){findOverlaps(x,z[1,])})
y=lapply(vmrs,function(x){findOverlaps(x,z[2,])})
pdf("HOOK2_mQTL_Plot.pdf")
plotRegion(BS_new2,
                z,col=pData(BS_new2)$GT_col,
                extend = 5000,annoTrack=list("mQTL"=SNP_res),geneTrack=exon,
                addRegions = list_of_VMRs[[1]])
plotRegion(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Anterior cingulate cortex (BA24)")],
                z,col=pData(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Anterior cingulate cortex (BA24)")])$GT_col,
                extend = 5000,annoTrack=list("mQTL"=SNP_res),geneTrack=exon,
                addRegions = list_of_VMRs[[1]])
plotRegion(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Caudate (basal ganglia)")],
                z,col=pData(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Caudate (basal ganglia)")])$GT_col,
                extend = 5000,annoTrack=list("mQTL"=SNP_res),geneTrack=exon,
                addRegions = list_of_VMRs[[1]])

plotRegion(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Frontal Cortex (BA9)")],
                z,col=pData(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Frontal Cortex (BA9)")])$GT_col,
                extend = 5000,annoTrack=list("mQTL"=SNP_res),geneTrack=exon,
                addRegions = list_of_VMRs[[1]])

plotRegion(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Putamen (basal ganglia)")],
                z,col=pData(BS_new2[,which(pData(BS_new2)$Tissue=="Brain - Putamen (basal ganglia)")])$GT_col,
                extend = 5000,annoTrack=list("mQTL"=SNP_res),geneTrack=exon,
                addRegions = list_of_VMRs[[1]])

dev.off()


## ALL THE TISSUES TOGETHER
a=match(pData(BS)$Participant_ID,rownames(try))
b=try[a,]
pData(BS)$GT=b
pData(BS)$GT_col="red"
pData(BS)$GT_col=replace(pData(BS)$GT_col,pData(BS)$GT==01,"black")
pData(BS)$GT_col=replace(pData(BS)$GT_col,pData(BS)$GT==02,"blue")
pdf("Common_mQTL_PlotAllTissues.pdf")
plotRegion(BS,
                vmrs[[1]][4063,],col=pData(BS)$GT_col,
                extend = 10000,annoTrack=list("mQTL"=SNP_res),geneTrack=exon,
                addRegions = vmrs[[1]])
dev.off()



### two questions: 1) how many of the snps per meQTL region are C_T 2) how many meQTLs overlap bw tissues

# 1 first

x=lapply(results,function(y){
	length(grep("C_T",y$snps))
	})

#
y=lapply(results,function(y){
	length(y$snps)
	})

# put in table per tissue
a1=unlist(t(t(y[[1]])))
a2=unlist(t(t(y[[2]])))
a3=unlist(t(t(y[[3]])))
a4=unlist(t(t(y[[4]])))
a5=unlist(t(t(y[[5]])))
a6=unlist(t(t(y[[6]])))
a7=unlist(t(t(y[[7]])))
a8=unlist(t(t(y[[8]])))


b1=unlist(t(t(x[[1]])))
b2=unlist(t(t(x[[2]])))
b3=unlist(t(t(x[[3]])))
b4=unlist(t(t(x[[4]])))
b5=unlist(t(t(x[[5]])))
b6=unlist(t(t(x[[6]])))
b7=unlist(t(t(x[[7]])))
b8=unlist(t(t(x[[8]])))

c1=b1/a1*100
c2=b2/a2*100
c3=b3/a3*100
c4=b4/a4*100
c5=b5/a5*100
c6=b6/a6*100
c7=b7/a7*100
c8=b8/a8*100

# made tables manually
d=list()
d[[1]]=cbind(a1,b1,c1)
d[[2]]=cbind(a2,b2,c2)
d[[3]]=cbind(a3,b3,c3)
d[[4]]=cbind(a4,b4,c4)
d[[5]]=cbind(a5,b5,c5)
d[[6]]=cbind(a6,b6,c6)
d[[7]]=cbind(a7,b7,c7)
d[[8]]=cbind(a8,b8,c8)
colnames(d[[1]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[2]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[3]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[4]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[5]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[6]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[7]])=c("all_snps","C_T_snps","Pct_CT_snps")
colnames(d[[8]])=c("all_snps","C_T_snps","Pct_CT_snps")

names(d)=names(results[1:8])

short=c("AMY","BA24","Caudate","BA9","HYPO","NAcc","Putamen","HC1")
wb <- createWorkbook()
Map(function(data, name){
 
    addWorksheet(wb, name)
    writeData(wb, name, data)
 
}, d, short)

saveWorkbook(wb, file = "meQTL_C_T_pct_50kb_200124.xlsx", overwrite = TRUE)