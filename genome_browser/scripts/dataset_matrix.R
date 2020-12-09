# Create a dataset matrix
# Peter Hickey
# Edited by Lindsay Rizzardi 
# 12-09-20

# ------------------------------------------------------------------------------
# Setup
#

library(SummarizedExperiment)
library(HDF5Array)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

extdir <- "../extdata"

# ------------------------------------------------------------------------------
# Load data
#

CG_BSseq_small_smooth <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(extdir, "eGTEx.Phase2_brain_samples.mCG.small_smooth"))

# ------------------------------------------------------------------------------
# Create data frames
#

datasets <- bind_rows(
  
  as.data.frame(
    colData(CG_BSseq_small_smooth)) %>%
    select(Participant_ID, Tissue) %>%
    mutate(Assay = "WGBS"))
  

tissue_colours <- as.data.frame(
  colData(CG_BSseq_small_smooth)) %>%
  select(Tissue, Brain_color) %>%
  distinct() %>%
  rename(Colour = Brain_color)

# ------------------------------------------------------------------------------
# Plot
#

individual_wide <- datasets %>%
  mutate(Assayed = TRUE) %>%
  complete(Participant_ID, Assay, Tissue, fill = list(Assayed = FALSE)) %>%
  group_by(Tissue) %>%
  filter(any(Assayed)) %>%
  ungroup() %>%
  mutate(Assayed = as.integer(Assayed)) %>%
  spread(-Assayed, Assayed)

individual_wide_matrix <- select(individual_wide, -Assay, -Tissue) %>%
  as.matrix()
rownames(individual_wide_matrix) <- select(
  individual_wide, Assay, Tissue) %>%
  mutate(Rowname = Tissue) %>%
  pull(Rowname)



hm <- Heatmap(
  matrix = t(individual_wide_matrix),
  col = c("grey", "darkblue"),
  name = "WGBS",
  rect_gp = gpar(col = "white"),
  row_title = "Donor",
  column_title = "Brain Region (NeuN+)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  heatmap_legend_param = list(labels = c("False", "True")))

pdf("../figures/individual-level_dataset_matrix.pdf")
hm 
dev.off()


