# Make gene models and 'flattened' genomic features from GENCODE 26
# annotations
# Peter Hickey
# 2016-10-26
# edited by Lindsay Rizzardi
# 2018-08-30

library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)

###################
### GENCODE v26 ###
###################


extdir <- "/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/"

# NOTE: makeTxDbFromGFF() only imports the "type", "gene_id", "transcript_id",
#       and "exon_id" columns, which is insufficient for our needs. We have to
#       therefore use import() + makeTxDbFromGRanges()
#       (see https://support.bioconductor.org/p/73158/#73288)

###-----------------------------------------------------------------------------
### Import FASTA files (protein-coding transcripts and lncRNAs tranxcripts)
###

fasta_pc <- import(file.path(extdir, "gencode.v26.pc_transcripts.fa.gz"),
                   format = "FASTA")
# > length(fasta_pc)
# [1] 95146

fasta_pc_tx_name <- sapply(strsplit(names(fasta_pc), "\\|"), "[[", 1)
fasta_pc_gene_name <- sapply(strsplit(names(fasta_pc), "\\|"), "[[", 2)

fasta_lnc <- import(file.path(extdir, 
                              "gencode.v26.lncRNA_transcripts.fa.gz"),
                    format = "FASTA")
# > length(fasta_lnc)
# [1] 27720

fasta_lnc_tx_name <- sapply(strsplit(names(fasta_lnc), "\\|"), "[[", 1)
fasta_lnc_gene_name <- sapply(strsplit(names(fasta_lnc), "\\|"), "[[", 2)

###-----------------------------------------------------------------------------
### Import GTF file as GRanges objects, only retaining those transcripts in the
### FASTA files
###

gtf <- import(file.path(extdir,  "gencode.v26.annotation.gtf.gz"))
# > length(gtf)
# [1] 2604486
       
fasta_id <- ifelse(gtf$transcript_id %in% sub("\\|.*", "", names(fasta_pc)),
                   "pc_transcripts",
                   ifelse(gtf$transcript_id %in% sub("\\|.*", "",
                                                     names(fasta_lnc)),
                          "lnc_transcripts",
                          NA_character_))
names(fasta_id) <- gtf$transcript_id
gtf$fasta_id <- fasta_id
gtf <- gtf[!is.na(gtf$fasta_id)]
# > length(gtf)
# [1]  2198897

###-----------------------------------------------------------------------------
### Add SeqInfo and drop non-autosomal seqlevels and
###

seqinfo(gtf) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
gtf <- keepSeqlevels(gtf, paste0("chr", 1:22),pruning.mode="coarse")
# > length(gtf)
# [1] 2123734
#
# > table(gtf$type)
#
#          gene     transcript           exon            CDS    start_codon 
#             0         118986         893194         685355          80181 
#    stop_codon            UTR Selenocysteine 
#         72221         273678            119 
#
#
# > length(unique(gtf$transcript_id))
# [1] 118986
#
# > length(unique(gtf$gene_id))
# [1] 34761

# NOTE: fasta_pc and fasta_lnc may both contain non-autosomal transcripts. We
#       don't want to remove these because it would affect the quasi-mapping.
#       However, we need to bear this in mind when looking for transcripts
#       that are 'missing' GTF annotations (since we only retain autosomal
#       GTF annotations).

###-----------------------------------------------------------------------------
### Convert to TxDb object
###

txdb <- makeTxDbFromGRanges(gr = gtf,
                            drop.stop.codons = FALSE,
                            metadata = NULL,
                            taxonomyId = 9606)
# TxDb object:
# Db type: TxDb
# Supporting package: GenomicFeatures
# Genome: hg38
# Organism: Homo sapiens
# Taxonomy ID: 9606
# transcript_nrow: 118986
# exon_nrow: 467642
# cds_nrow: 261781
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2018-10-08 10:16:27 -0400 (Mon, 08 Oct 2018)
# GenomicFeatures version at creation time: 1.33.2
# RSQLite version at creation time: 2.1.1
# DBSCHEMAVERSION: 1.2


###-----------------------------------------------------------------------------
### Extract features of interest
###

transcripts <- transcripts(txdb, columns = columns(txdb))
# > length(transcripts)
# [1] 118986
genes <- genes(txdb, columns = columns(txdb))
# > length(genes)
# [1] 34761
promoters <- promoters(txdb,
                       columns = columns(txdb),
                       upstream = 2000,
                       downstream = 2000)
# > length(promoters)
# [1] 118986
transcripts_by_gene <- transcriptsBy(txdb, by = "gene")
# > length(transcripts_by_gene)
# [1] 34761
five_utrs_by_transcript <- fiveUTRsByTranscript(txdb, use.names = TRUE)
# > length(five_utrs_by_transcript)
# [1] 76577
three_utrs_by_transcript <- threeUTRsByTranscript(txdb, use.names = TRUE)
# > length(three_utrs_by_transcript)
# [1] 69316
exons_by_transcript <- exonsBy(txdb, by = "tx", use.names = TRUE)
# > length(exons_by_transcript)
# [1]  118986
introns_by_transcript <- intronsByTranscript(txdb, use.names = TRUE)
# > length(introns_by_transcript)
# [1] 118986

###-----------------------------------------------------------------------------
### 'Unflattened' features
###

unflattened_features <- list(transcripts = transcripts,
                             genes = genes,
                             promoters = promoters,
                             transcripts_by_gene = transcripts_by_gene,
                             five_utrs_by_transcript = five_utrs_by_transcript,
                             three_utrs_by_transcript = three_utrs_by_transcript,
                             exons_by_transcript = exons_by_transcript,
                             introns_by_transcript = introns_by_transcript)

###-----------------------------------------------------------------------------
### 'Unflattened features' based on protein-coding transcript sequences
###

unflattened_features_pc_transcripts <-
  list(transcripts = transcripts[transcripts$TXNAME %in% fasta_pc_tx_name],
       genes = genes[unlist(genes$GENEID) %in% fasta_pc_gene_name],
       promoters = promoters[promoters$TXNAME %in% fasta_pc_tx_name],
       transcripts_by_gene = transcripts_by_gene[
         na.omit(match(fasta_pc_gene_name,
                       names(transcripts_by_gene)))],
       five_utrs_by_transcript = five_utrs_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(five_utrs_by_transcript)))],
       three_utrs_by_transcript = three_utrs_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(three_utrs_by_transcript)))],
       exons_by_transcript = exons_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(exons_by_transcript)))],
       introns_by_transcript = introns_by_transcript[
         na.omit(match(fasta_pc_tx_name, names(introns_by_transcript)))])

###-----------------------------------------------------------------------------
### 'Unflattened features' based on lncRNA transcript sequences
###

unflattened_features_lnc_transcripts <-
  list(transcripts = transcripts[transcripts$TXNAME %in% fasta_lnc_tx_name],
       genes = genes[unlist(genes$GENEID) %in% fasta_lnc_gene_name],
       promoters = promoters[promoters$TXNAME %in% fasta_lnc_tx_name],
       transcripts_by_gene = transcripts_by_gene[
         na.omit(match(fasta_lnc_gene_name,
                       names(transcripts_by_gene)))],
       five_utrs_by_transcript = five_utrs_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(five_utrs_by_transcript)))],
       three_utrs_by_transcript = three_utrs_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(three_utrs_by_transcript)))],
       exons_by_transcript = exons_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(exons_by_transcript)))],
       introns_by_transcript = introns_by_transcript[
         na.omit(match(fasta_lnc_tx_name, names(introns_by_transcript)))])

###-----------------------------------------------------------------------------
### Sanity checking UTRs

# Check that no lncRNA has 5' or 3' UTRs
# > any(names(five_utrs_by_transcript) %in% sub("\\|.*", "", names(fasta_lnc)))
# [1] FALSE
#
# > any(names(three_utrs_by_transcript) %in% sub("\\|.*", "", names(fasta_lnc)))
# [1] FALSE

###-----------------------------------------------------------------------------
### Find the transcripts that lack 5' UTRs
###

no_five_utr <- transcripts[!transcripts$TXNAME %in%
                             names(five_utrs_by_transcript)]

# > table(gtf$fasta_id, gtf$transcript_id %in% no_five_utr$TXNAME)
#
#                      FALSE    TRUE
#  lnc_transcripts       0  108694
#  pc_transcripts  1743890  271150


# > table(gtf$transcript_type, gtf$transcript_id %in% no_five_utr$TXNAME)
#
#                                FALSE    TRUE
#  3prime_overlapping_ncRNA            0     103
#  antisense                           0   44184
#  bidirectional_promoter_lncRNA       0      87
#  IG_C_gene                           0     175
#  IG_D_gene                           0      67
#  IG_J_gene                          10      48
#  IG_V_gene                         846     140
#  lincRNA                             0   52615
#  macro_lncRNA                        0       2
#  misc_RNA                            0      24
#  non_coding                          0       6
#  non_stop_decay                   1029     325
#  nonsense_mediated_decay        228634   73329
#  polymorphic_pseudogene             70      78
#  processed_pseudogene                0       2
#  processed_transcript                0    3139
#  protein_coding                1512716  196536
#  pseudogene                          0      74
#  retained_intron                     0    2321
#  sense_intronic                      0    2652
#  sense_overlapping                   0    1288
#  TEC                                 0    2197
#  TR_C_gene                           0      64
#  TR_D_gene                           0      12
#  TR_J_gene                           0     234
#  TR_V_gene                         585     142

# > table(is.na(gtf$ccdsid), gtf$transcript_id %in% no_five_utr$TXNAME)
#
#            FALSE    TRUE
#  FALSE 1000499   46437
#  TRUE   743391  333407

###-----------------------------------------------------------------------------
### Find the transcripts that lack 3' UTRs
###

no_three_utr <- transcripts[!transcripts$TXNAME %in%
                              names(three_utrs_by_transcript)]

# > table(gtf$fasta_id, gtf$transcript_id %in% no_three_utr$TXNAME)
#
#                   FALSE    TRUE
# lnc_transcripts       0  108694
# pc_transcripts  1691317  323723

# > table(gtf$transcript_type, gtf$transcript_id %in% no_three_utr$TXNAME)
#
#                               FALSE    TRUE
# 3prime_overlapping_ncRNA            0     103
# antisense                           0   44184
# bidirectional_promoter_lncRNA       0      87
# IG_C_gene                         175       0
# IG_D_gene                          10      57
# IG_J_gene                           0      58
# IG_V_gene                          17     969
# lincRNA                             0   52615
# macro_lncRNA                        0       2
# misc_RNA                            0      24
# non_coding                          0       6
# non_stop_decay                      0    1354
# nonsense_mediated_decay        301892      71
# polymorphic_pseudogene             61      87
# processed_pseudogene                0       2
# processed_transcript                0    3139
# protein_coding                1389098  320154
# pseudogene                          0      74
# retained_intron                     0    2321
# sense_intronic                      0    2652
# sense_overlapping                   0    1288
# TEC                                 0    2197
# TR_C_gene                          64       0
# TR_D_gene                           0      12
# TR_J_gene                           0     234
# TR_V_gene                           0     727

# > table(is.na(gtf$ccdsid), gtf$transcript_id %in% no_three_utr$TXNAME)
#
#          FALSE    TRUE
#  FALSE 1011861   35075
#  TRUE   679456  397342


###-----------------------------------------------------------------------------
### Find the transcripts that have 5' UTRs but lack 3' UTRs
###

have_five_utr_no_three_utr <- transcripts[transcripts$TXNAME %in%
                                            names(five_utrs_by_transcript) &
                                            !transcripts$TXNAME %in%
                                            names(three_utrs_by_transcript)]

# > length(have_five_utr_no_three_utr)
# [1] 18331

# > table(gtf$fasta_id, gtf$transcript_id %in% have_five_utr_no_three_utr$TXNAME)
#
#                     FALSE    TRUE
#  lnc_transcripts  108694       0
#  pc_transcripts  1760689  254351

# > table(gtf$transcript_type, gtf$transcript_id %in% have_five_utr_no_three_utr$TXNAME)
#
#                                FALSE    TRUE
# 3prime_overlapping_ncRNA          103       0
# antisense                       44184       0
# bidirectional_promoter_lncRNA      87       0
# IG_C_gene                         175       0
# IG_D_gene                          67       0
# IG_J_gene                          48      10
# IG_V_gene                         149     837
# lincRNA                         52615       0
# macro_lncRNA                        2       0
# misc_RNA                           24       0
# non_coding                          6       0
# non_stop_decay                    325    1029
# nonsense_mediated_decay        301904      59
# polymorphic_pseudogene            127      21
# processed_pseudogene                2       0
# processed_transcript             3139       0
# protein_coding                1457442  251810
# pseudogene                         74       0
# retained_intron                  2321       0
# sense_intronic                   2652       0
# sense_overlapping                1288       0
# TEC                              2197       0
# TR_C_gene                          64       0
# TR_D_gene                          12       0
# TR_J_gene                         234       0
# TR_V_gene                         142     585
###-----------------------------------------------------------------------------
### 'Flatten' features by unlisting, unstranding, and reducing
### NOTE: This procedure means that a base may belong to multiple categories
###

# NOTE: Intergenic regions are defined as those without a gene on either strand
intergenic <- gaps(reduce(unstrand(genes)))
intergenic <- intergenic[strand(intergenic) == "*"]

flattened_features <- list(genic = reduce(unstrand(genes)),
                           promoter = reduce(unstrand(promoters)),
                           five_utr =
                             reduce(unstrand(unlist(five_utrs_by_transcript))),
                           three_utr =
                             reduce(unstrand(unlist(three_utrs_by_transcript))),
                           exonic =
                             reduce(unstrand(unlist(exons_by_transcript))),
                           intronic =
                             reduce(unstrand(unlist(introns_by_transcript))),
                           intergenic = intergenic)

#-------------------------------------------------------------------------------
# 'Flattened features' based on protein-coding transcript sequences
#

intergenic_pc_transcripts <-
  gaps(reduce(unstrand(genes[names(genes) %in% fasta_pc_gene_name])))
intergenic_pc_transcripts <- intergenic_pc_transcripts[
  strand(intergenic_pc_transcripts) == "*"]

flattened_features_pc_transcripts <-
  list(genic = reduce(unstrand(genes[names(genes) %in% fasta_pc_gene_name])),
       promoter = reduce(unstrand(promoters[promoters$TXNAME %in%
                                              fasta_pc_tx_name])),
       five_utr =
         reduce(unstrand(unlist(five_utrs_by_transcript[
           names(five_utrs_by_transcript) %in% fasta_pc_tx_name]))),
       three_utr =
         reduce(unstrand(unlist(three_utrs_by_transcript[
           names(three_utrs_by_transcript) %in% fasta_pc_tx_name]))),
       exonic =
         reduce(unstrand(unlist(exons_by_transcript[
           names(exons_by_transcript) %in% fasta_pc_tx_name]))),
       intronic =
         reduce(unstrand(unlist(introns_by_transcript[
           names(introns_by_transcript) %in% fasta_pc_tx_name]))),
       intergenic = intergenic_pc_transcripts)

#-------------------------------------------------------------------------------
# 'Flattened features' based on long non-coding RNA transcript sequences
#

intergenic_lnc_transcripts <-
  gaps(reduce(unstrand(genes[names(genes) %in% fasta_lnc_gene_name])))
intergenic_lnc_transcripts <- intergenic_lnc_transcripts[
  strand(intergenic_lnc_transcripts) == "*"]

flattened_features_lnc_transcripts <-
  list(genic = reduce(unstrand(genes[names(genes) %in% fasta_lnc_gene_name])),
       promoter = reduce(unstrand(promoters[promoters$TXNAME %in%
                                              fasta_lnc_tx_name])),
       five_utr =
         reduce(unstrand(unlist(five_utrs_by_transcript[
           names(five_utrs_by_transcript) %in% fasta_lnc_tx_name]))),
       three_utr =
         reduce(unstrand(unlist(three_utrs_by_transcript[
           names(three_utrs_by_transcript) %in% fasta_lnc_tx_name]))),
       exonic =
         reduce(unstrand(unlist(exons_by_transcript[
           names(exons_by_transcript) %in% fasta_lnc_tx_name]))),
       intronic =
         reduce(unstrand(unlist(introns_by_transcript[
           names(introns_by_transcript) %in% fasta_lnc_tx_name]))),
       intergenic = intergenic_lnc_transcripts)

#-------------------------------------------------------------------------------
# Save unflattened features objects
#

save(unflattened_features,
     unflattened_features_pc_transcripts,
     unflattened_features_lnc_transcripts,
     file = "/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/unflattened-GENCODE-v26-features.rda",
     compress = "xz")

#-------------------------------------------------------------------------------
# Save flattened features objects
#

save(flattened_features,
     flattened_features_lnc_transcripts,
     flattened_features_pc_transcripts,
     file = "/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/flattened-GENCODE-v26-features.rda")

#-------------------------------------------------------------------------------
# Make object of genes/transcripts for use in plotting routines
#

saveRDS(transcripts_by_gene,
        file = "/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/GENCODE-v26-transcripts-by-gene.rds")
saveRDS(sort(granges(genes)),
        file = "/dcl01/FB2/data/personal/gtex/gtex/gencode_v26/GENCODE-v26-genes.rds")

