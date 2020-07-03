# WGBS summary statistics
# Peter Hickey
# 2017-01-26
# Edited for eGTEx by Lindsay Rizzardi
# 2018-07-10
# need conda_R/3.5.x
library(dplyr)
library(purrr)
library(bsseq)
library(readr)
library(HDF5Array)
library(matrixStats)
library(BSgenome.Hsapiens.UCSC.hg38)

#-------------------------------------------------------------------------------
# Parse files
#

# Parse Bismark output and import into R
sequenced <- system("grep 'Sequence pairs analysed in total' /dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/HiSeq*/*/*PE_report.txt", intern = TRUE)
aligned <- system("grep 'Number of paired-end alignments with a unique best hit' /dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/HiSeq*/*/*PE_report.txt", intern = TRUE)

# More parsing of Bismark output (now within R)
sequenced_df <- map_df(strsplit(sequenced, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 10L),
             `Number sequenced PE reads` = as.numeric(x[[2]]))
})
aligned_df <- map_df(strsplit(aligned, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 10L),
             `Number aligned PE reads` = as.numeric(x[[2]]))
})

duplicated <- system("grep 'Total number duplicated alignments removed' /dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/HiSeq*/*/*deduplication_report.txt", intern = TRUE)
duplicated_df <- map_df(strsplit(duplicated, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 10L),
             `Number duplicated alignments` = map_chr(strsplit(x[[2]]," "),`[[`,1L))
})
duplicated_df$`Number duplicated alignments`=as.numeric(duplicated_df$`Number duplicated alignments`)

#-------------------------------------------------------------------------------
# Aggregate at sample level
#

sequenced_df <- sequenced_df %>%
  group_by(ID) %>%
  summarise(`Number sequenced PE reads` = sum(`Number sequenced PE reads`))
aligned_df <- aligned_df %>%
  group_by(ID) %>%
  summarise(`Number aligned PE reads` = sum(`Number aligned PE reads`))

#-------------------------------------------------------------------------------
# Alignment rate
#

df <- inner_join(sequenced_df, aligned_df, duplicated_df, by = c("ID" = "ID"))
df <- df %>%
  mutate(`Alignment rate (%)` =
           round(100 * `Number aligned PE reads` / `Number sequenced PE reads`,
                 0))
df <- inner_join(df, duplicated_df, by = c("ID" = "ID")) 
#-------------------------------------------------------------------------------
# Bisulfite-conversion rate
#

# Load objects
BS2=HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_brain_samples.mCG.small_smooth")
BS3=HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Sorted_GTEx_Brain_hg38/BSseq/eGTEx.Phase2_non-brain_samples.mCG.small_smooth")


# Bisulfite-conversion rates
lambda <- GRanges("lambda",
                  IRanges(1, 48502))
bs_conversion_rate_df <-
  data_frame(ID = c(colnames(BS2),
                    colnames(BS3)),
             `Bisulfite conversion rate (%)` =
               round(100 *
                       c(1 - colSums(getCoverage(BSseq = BS2,
                                                 regions = lambda,
                                                 type = "M")[[1L]]) /
                           colSums(getCoverage(BSseq = BS2,
                                               regions = lambda,
                                               type = "Cov")[[1L]]),
                         1 - colSums(getCoverage(BSseq = BS3,
                                                 regions = lambda,
                                                 type = "M")[[1L]]) /
                           colSums(getCoverage(BSseq = BS3,
                                               regions = lambda,
                                               type = "Cov")[[1L]])),
                     2))
df <- inner_join(df, bs_conversion_rate_df)


#-------------------------------------------------------------------------------
# Coverage-related metrics
#

# Keep std chr only
BS2=keepStandardChromosomes(BS2,pruning.mode="coarse")
BS2=dropSeqlevels(BS2,"chrM",pruning.mode="coarse")

BS3=keepStandardChromosomes(BS3,pruning.mode="coarse")
BS3=dropSeqlevels(BS3,"chrM",pruning.mode="coarse")

# Number of CpGs on autosomes in hg38
n_CpGs <- Reduce(sum, bsapply(BSParams = new("BSParams",
                                             X = BSgenome.Hsapiens.UCSC.hg38,
                                             FUN = countPattern,
                                             exclude = c("M", "_")),
                              pattern = "CG"))

# Compute coverage-related metrics
coverage_df <-
  data_frame(ID = c(colnames(BS2),
                    colnames(BS3)),
             `Number covered CpGs` =
               c(colSums(getCoverage(BS2) > 0),
                 colSums(getCoverage(BS3) > 0)),
             `Covered CpGs (%)` =
               round(100 * `Number covered CpGs` / n_CpGs, 0),
             `Mean depth` = c(colMeans(getCoverage(BS2)),
                              colMeans(getCoverage(BS3))))

df <- inner_join(df, coverage_df, c("ID" = "ID"))

#-------------------------------------------------------------------------------
# Arrange columns and split into sorted and unsorted data
#
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis")
df <- select(df, ID, `Number sequenced PE reads`, `Number aligned PE reads`,
             `Alignment rate (%)`, `Number duplicated alignments`,`Number covered CpGs`, `Covered CpGs (%)`,
             `Mean depth`, `Bisulfite conversion rate (%)`)

df <- df %>%
  mutate(`Duplication rate (%)` =
           round(100 * `Number duplicated alignments` / `Number aligned PE reads`,
                 0))


pData=rbind(pData(BS2),pData(BS3))

df_list <- split(df, grepl("Lung|Thyroid", pData$Tissue))
#write_csv(x = df,
#          path = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/tables/Summary_of_Sorted_WGBS.csv",
#          col_names = TRUE)

write_csv(x = df_list[[1]],
          path = "tables/Summary_of_Sorted_WGBS.csv",
          col_names = TRUE)

df_list1=df_list[[2]]


#####################################################################
#########  Repeat for Bulk_GTEx_Brain and GTEx_nonbrain  ############
#####################################################################


# Parse Bismark output and import into R
sequenced <- system("grep 'Sequence pairs analysed in total' /dcl01/FB2/data/personal/gtex/gtex/Bulk_GTEx_Brain_hg38/HiSeq*/*/*PE_report.txt", intern = TRUE)
aligned <- system("grep 'Number of paired-end alignments with a unique best hit' /dcl01/FB2/data/personal/gtex/gtex/Bulk_GTEx_Brain_hg38/HiSeq*/*/*PE_report.txt", intern = TRUE)

sequenced2 <- system("grep 'Sequence pairs analysed in total' /dcl01/FB2/data/personal/gtex/gtex/GTEx_nonbrain_hg38/HiSeq*/*/*PE_report.txt", intern = TRUE)
aligned2 <- system("grep 'Number of paired-end alignments with a unique best hit' /dcl01/FB2/data/personal/gtex/gtex/GTEx_nonbrain_hg38/HiSeq*/*/*PE_report.txt", intern = TRUE)

sequenced<-c(sequenced,sequenced2)
aligned<-c(aligned,aligned2)

# More parsing of Bismark output (now within R)
sequenced_df <- map_df(strsplit(sequenced, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 10L),
             `Number sequenced PE reads` = as.numeric(x[[2]]))
})
aligned_df <- map_df(strsplit(aligned, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 10L),
             `Number aligned PE reads` = as.numeric(x[[2]]))
})

duplicated <- system("grep 'Total number duplicated alignments removed' /dcl01/FB2/data/personal/gtex/gtex/Bulk_GTEx_Brain_hg38/HiSeq*/*/*deduplication_report.txt", intern = TRUE)
duplicated2 <- system("grep 'Total number duplicated alignments removed' /dcl01/FB2/data/personal/gtex/gtex/GTEx_nonbrain_hg38/HiSeq*/*/*deduplication_report.txt", intern = TRUE)
duplicated<-c(duplicated,duplicated2)

duplicated_df <- map_df(strsplit(duplicated, "\t"), function(x) {
  data_frame(ID = map_chr(strsplit(x[[1]], "/"), `[[`, 10L),
             `Number duplicated alignments` = map_chr(strsplit(x[[2]]," "),`[[`,1L))
})
duplicated_df$`Number duplicated alignments`=as.numeric(duplicated_df$`Number duplicated alignments`)

#-------------------------------------------------------------------------------
# Aggregate at sample level
#

sequenced_df <- sequenced_df %>%
  group_by(ID) %>%
  summarise(`Number sequenced PE reads` = sum(`Number sequenced PE reads`))
aligned_df <- aligned_df %>%
  group_by(ID) %>%
  summarise(`Number aligned PE reads` = sum(`Number aligned PE reads`))
duplicated_df <- duplicated_df %>%
  group_by(ID) %>%
  summarise(`Number duplicated alignments` = sum(`Number duplicated alignments`))

#-------------------------------------------------------------------------------
# Alignment rate
#

df <- inner_join(sequenced_df, aligned_df, duplicated_df, by = c("ID" = "ID"))
df <- df %>%
  mutate(`Alignment rate (%)` =
           round(100 * `Number aligned PE reads` / `Number sequenced PE reads`,
                 0))
df <- inner_join(df, duplicated_df, by = c("ID" = "ID")) 
#-------------------------------------------------------------------------------
# Bisulfite-conversion rate
#

# Load objects
BS4=HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/GTEx_nonbrain_hg38/BSseq/eGTEx.Phase1_non-brain_samples.mCG.small_smooth")
BS5=HDF5Array::loadHDF5SummarizedExperiment("/dcl01/FB2/data/personal/gtex/gtex/Bulk_GTEx_Brain_hg38/BSseq/eGTEx.Phase1_brain_samples.mCG.small_smooth")


# Bisulfite-conversion rates
lambda <- GRanges("lambda",
                  IRanges(1, 48502))
bs_conversion_rate_df <-
  data_frame(ID = c(colnames(BS4),
                    colnames(BS5)),
             `Bisulfite conversion rate (%)` =
               round(100 *
                       c(1 - colSums(getCoverage(BSseq = BS4,
                                                 regions = lambda,
                                                 type = "M")[[1L]]) /
                           colSums(getCoverage(BSseq = BS4,
                                               regions = lambda,
                                               type = "Cov")[[1L]]),
                         1 - colSums(getCoverage(BSseq = BS5,
                                                 regions = lambda,
                                                 type = "M")[[1L]]) /
                           colSums(getCoverage(BSseq = BS5,
                                               regions = lambda,
                                               type = "Cov")[[1L]])),
                     2))
df <- inner_join(df, bs_conversion_rate_df)


#-------------------------------------------------------------------------------
# Coverage-related metrics
#

# Keep std chr only
BS4=keepStandardChromosomes(BS4,pruning.mode="coarse")
BS4=dropSeqlevels(BS4,"chrM",pruning.mode="coarse")

BS5=keepStandardChromosomes(BS5,pruning.mode="coarse")
BS5=dropSeqlevels(BS5,"chrM",pruning.mode="coarse")

# Compute coverage-related metrics
coverage_df <-
  data_frame(ID = c(colnames(BS4),
                    colnames(BS5)),
             `Number covered CpGs` =
               c(colSums(getCoverage(BS4) > 0),
                 colSums(getCoverage(BS5) > 0)),
             `Covered CpGs (%)` =
               round(100 * `Number covered CpGs` / n_CpGs, 0),
             `Mean depth` = c(colMeans(getCoverage(BS4)),
                              colMeans(getCoverage(BS5))))

df <- inner_join(df, coverage_df, c("ID" = "ID"))

#-------------------------------------------------------------------------------
# Arrange columns and split into sorted and unsorted data
#
setwd("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis")
df <- select(df, ID, `Number sequenced PE reads`, `Number aligned PE reads`,
             `Alignment rate (%)`, `Number duplicated alignments`,`Number covered CpGs`, `Covered CpGs (%)`,
             `Mean depth`, `Bisulfite conversion rate (%)`)

df <- df %>%
  mutate(`Duplication rate (%)` =
           round(100 * `Number duplicated alignments` / `Number aligned PE reads`,
                 0))

pData2=rbind(pData(BS4),pData(BS5))

df_list <- split(df, grepl("Brain", pData2$Tissue))

df_list3<-rbind(df_list[[1]],df_list1)

write_csv(x = df_list3,
          path = "tables/Summary_of_eGTEx_nonbrain_WGBS.csv",
          col_names = TRUE)
write_csv(x = df_list[[2]],
          path = "tables/Summary_of_eGTEx_Bulk_Brain_WGBS.csv",
          col_names = TRUE)


