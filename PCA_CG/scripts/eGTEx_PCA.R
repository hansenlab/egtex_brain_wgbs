library(bsseq)
library(HDF5Array)
library(DelayedArray) 
library(DelayedMatrixStats)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(plotly)
library(htmlwidgets)
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/PCA/")
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)


PCA <- function(BSseq, min_width = 1000, min_loci = 10, cov_variable= NULL,type=c("raw","smooth")) {
  # TODO: This means the min_cov requirement must hold for all samples; is
  #       this what we want?
  if (!is.null(cov_variable)) {
    message("Getting coverage...")
    M<-as.matrix(getMeth(BSseq,type=type))
    M<-M[cov_variable,]
    BSseq=BSseq[cov_variable,]
  }
  else {
  message("Generating M without coverage...\n")
  M<-as.matrix(getMeth(BSseq,type=type)) }
  # Function to compute methylation level in bins along genome
                          
    if (is.unsorted(BSseq)) {
      message("Sorting BSseq ...\n")
      BSseq <- sort(BSseq)
    }
    # Make regions containing at least min_loci and at least min_width wide
    # TODO: This is godawful slow
    grl <- split(rowRanges(BSseq), seqnames(BSseq))
    message("Generating regions...\n")
    regions <- unlist(endoapply(grl, function(gr) {
      # TODO: Print actual seqnames (not integer representation underlying
      #       the factor).
      # TODO: Add verbose option to report this output
      message(unique(seqnames(gr)))
      ss <- start(gr)
      s <- rep(NA_integer_, length(gr))
      e <- rep(NA_integer_, length(gr))
      i <- 1L
      j <- 1L
      while (i <= length(gr)) {
        k <- min_loci
        while ((i + k) <= length(gr)) {
          # cat("i =", i, "\n")
          # cat("j =", j, "\n")
          # cat("k =", k, "\n")
          # cat("i + k =", i + k, "\n")
          if ((ss[i + k - 1L] - ss[i]) >= min_width) {
            s[j] <- ss[i]
            e[j] <- ss[i + k - 1L]
            j <- j + 1L
            i <- i + k
            # Sys.sleep(5)
            break
          } else {
            k <- k + 1L
          }
        }
        i <- i + 1L
      }
      GRanges(unique(seqnames(gr)),
              IRanges(na.omit(s), na.omit(e)))
    }))

    # Compute avereage methylation level in each bin
    message("Computing methylation\n")

    ov <- findOverlaps(granges(BSseq), regions)
    meth=M[queryHits(ov),]
    meth_list=split(meth, subjectHits(ov))

    out <- lapply(meth_list, matrix, ncol = ncol(M))

    out <- do.call(rbind, lapply(out, colMeans2, na.rm = TRUE))
    outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
    colnames(outMatrix) <- sampleNames(BSseq)
    outMatrix[as.integer(rownames(out)), ] <- out
    .DelayedMatrix <- function(x) {
    x_name <- deparse(substitute(x))
    X <- try(DelayedArray(x), silent = TRUE)
    if (is(X, "try-error")) {
        stop("Could not construct DelayedMatrix from '", x_name, "'",
             call. = FALSE)
    }
    if (!is(X, "DelayedMatrix")) {
        stop("'", x_name, "' must be matrix-like", call. = FALSE)
    }
    X
}

    meth<-.DelayedMatrix(outMatrix)

  message("Computing PCA")
  stats::prcomp(t(meth))
}




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

BS.cov=as.matrix(getCoverage(bsseq))
colnames(BS.cov)=colnames(bsseq)

### Subset by CpGs with at least 1 read in 90% of samples from each tissue 
pData=pData(bsseq)
all=NULL
tissue=as.character(levels(factor(pData$Tissue)))
for (t in 1:length(tissue)){

tiss=pData[which(pData$Tissue==tissue[t]),]
BS.cov_tiss=BS.cov[,which(colnames(BS.cov) %in% rownames(tiss))]
tiss.all <- rowSums(BS.cov_tiss >= 1) >= (round(0.9*ncol(BS.cov_tiss))) 
all=as.data.frame(cbind(all,tiss.all))
}
x=matrixStats::rowAlls(as.matrix(all))
pca=PCA(bsseq,type="smooth",cov_variable=x)

percent_variance <- setNames(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2),
                               colnames(pca$x))
library(plotly)
df <- data.frame(x = pca$x[, "PC1"],y = pca$x[, "PC2"])

pdf("eGTEx_Sorted_PCA_QC_Pass_1kb.pdf",width=12,height=10)

theme_set(theme_bw(base_size=16))
g <- ggplot(data = df, aes(x = x, y = y, colour = pData(bsseq)$Tissue, label=pData(bsseq)$Sample_ID)) + coord_fixed(1)+
      geom_point(size=3.5)+scale_colour_manual(values = c("Brain - Hippocampus"="#1b9e77","Brain - Amygdala" ="#e7298a","Brain - Putamen (basal ganglia)" = "#7570b3","Brain - Hypothalamus" = "#66a61e", "Brain - Nucleus accumbens (basal ganglia)" = "#d95f02", "Brain - Anterior cingulate cortex (BA24)" = "#a6761d","Brain - Frontal Cortex (BA9)" = "#e6ab02", "Brain - Caudate (basal ganglia)" = "#666666"))+
      labs(colour="Tissue")
g + xlab(paste0("PC1", " (", round(percent_variance["PC1"], 1), " %)")) +
    ylab(paste0("PC2", " (", round(percent_variance["PC2"], 1), " %)")) 

htmlwidgets::saveWidget(ggplotly(g),"eGTEx_Sorted_PCA_QC_Pass_1kb.html")
dev.off()
pca_GTEx_sorted=pca
################
#### now do with NN HC and GTEx HC
BS_HC=bsseq[,bsseq$group=="Hippocampus"]


## Fix pData to look the same
pData(BS_HC)$Sample_Type="NeuN_pos"
pData(BS_HC)$Project="GTEx"
pData_BS_original=pData(BS_HC)
pData_GTEx=pData(BS_HC)[,c(34,39,41)]
pData(BS_HC)=pData_GTEx
colnames(pData(BS_HC))=c("Sample_ID","Tissue","Project")

BS2=HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/BrainEpigenome_hg38/BSseq/BrainEpigenome_hg38.mCG.small_smooth")
BS2=keepStandardChromosomes(BS2,pruning.mode="coarse")
BS2=dropSeqlevels(BS2,c("chrX","chrY","chrM"),pruning.mode="coarse")
BS2=sort(sortSeqlevels(BS2))

pData(BS2)$Project="NatNeuro"
pData(BS2)$Tissue=replace(pData(BS2)$Tissue,pData(BS2)$Tissue=="BA9","Brain - Frontal Cortex (BA9)")
pData(BS2)$Tissue=replace(pData(BS2)$Tissue,pData(BS2)$Tissue=="HC","Brain - Hippocampus")
pData(BS2)$Tissue=replace(pData(BS2)$Tissue,pData(BS2)$Tissue=="BA24","Brain - Anterior cingulate cortex (BA24)")
pData(BS2)$Tissue=replace(pData(BS2)$Tissue,pData(BS2)$Tissue=="NA","Brain - Nucleus accumbens (basal ganglia)")

pData(BS2)=pData(BS2)[,c("Sample_ID","Tissue","Project")]
pData(BS2)$Type=sapply(strsplit(sampleNames(BS2),"[_]"),"[[",3)
## Just get sorted pos
BS2=BS2[,pData(BS2)$Type =="pos"]
BS2=BS2[,pData(BS2)$Tissue =="Brain - Hippocampus"]
pData(BS2)$Type=NULL

#Combine objects:
BS3=combineList(BS_HC,BS2)
BS.cov=as.matrix(getCoverage(BS3))
### Subset by CpGs with at least 1 read in 90% of samples from each tissue 

pData=pData(BS3)
all=NULL
tissue=as.character(levels(factor(pData$Tissue)))
for (t in 1:length(tissue)){

tiss=pData[which(pData$Tissue==tissue[t]),]
BS.cov_tiss=BS.cov[,which(colnames(BS.cov) %in% rownames(tiss))]
tiss.all <- rowSums(BS.cov_tiss >= 1) >= (round(0.9*ncol(BS.cov_tiss)))
all=as.data.frame(cbind(all,tiss.all))
}

x=matrixStats::rowAlls(as.matrix(all))

pca=PCA(BS3,cov_variable=x,type="smooth")
pca_eGTEx_NN_HC=pca

percent_variance <- setNames(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2),
                               colnames(pca$x))

df <- data.frame(x = pca$x[, "PC1"],y = pca$x[, "PC2"])
pdf("HC_eGTEx_NN_Sorted_PCA_QC_1kb.pdf",width=12,height=10)
theme_set(theme_bw(base_size=16))
g <- ggplot(data = df, aes(x = x, y = y, label = pData(BS3)$Sample_ID, shape=pData(BS3)$Project,colour = pData(BS3)$Tissue)) +coord_fixed(1)+
      geom_point(size=3.5)+scale_colour_manual(values = c("Brain - Hippocampus"="#1b9e77"))+
      labs(colour="Tissue",shape="Project")
g + xlab(paste0("PC1", " (", round(percent_variance["PC1"], 1), " %)")) +
    ylab(paste0("PC2", " (", round(percent_variance["PC2"], 1), " %)")) 
dev.off()

saveRDS(pca_eGTEx_NN_HC,pca_GTEx_sorted,file="New_CG_PCAs.rds")





