# Make data frame of metadata for eGTEx samples
# Peter Hickey
# 2018-07-26

# NOTE: Creates one row per combination of
#       `Sample_ID_for_Data_Sharing_and_Public_Release` and `HiSeq_Run`; this
#       corresponds to having one row per Bismark output file.

# Setup ------------------------------------------------------------------------

library(S4Vectors)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(jsonlite)
library(tidyr)
library(here)

# Read CSVs --------------------------------------------------------------------

# Phase 1 data
cd1 <- read_csv(
  file = here("metadata", "data", "GTEx_pData.csv"),
  col_types = "ccccccccccciiiii")
# Phase 2 data
cd2 <- read_csv(here("metadata", "data", "GTEx_Flow_Sorted_pData.csv"))

# Normalise colnames -----------------------------------------------------------

# NOTE: Some columns renamed for consistency with Phase 2 data
# NOTE: Some columns dropped due to redundancy
cd1 <- cd1 %>%
  rename(
    Sample_ID_for_Data_Sharing_and_Public_Release = Sample_Name,
    Average_Library_Length = `Library_Length(avg)`,
    Self_Reported_Sex = Gender) %>%
  select(-Tissue_Name)

# NOTE Some columns renamed for consistency with Phase 1 data
# NOTE: Some columns added for consistency with Phase 2 data
colnames(cd2) <- colnames(cd2) %>%
  gsub(pattern = "\\.+", " ", x = .) %>%
  str_trim() %>%
  gsub(pattern = "\\s", replacement = "_", x = .)
cd2 <- cd2 %>%
  rename(
    Sample_ID_for_Data_Sharing_and_Public_Release =
      `Sample_ID_for_data_sharing_and_public_release`,
    Tissue = Tissue_Site_Detail,
    Root_Sample = Root_Sample_s,
    Participant_ID = Participant_ID_s,
    MasterPure_Lot = MasterPure_lot,
    Average_Library_Length = `Library_Length(avg)`,
    Self_Reported_Sex = Gender) %>%
  mutate(`Library_Prep_Batch` = NA_integer_)

# Fix column classes and formatting --------------------------------------------

cd1 <- cd1 %>%
  mutate(
    Library_Prep_Date = as.Date(
      as.character(Library_Prep_Date),
      format = "%Y%m%d"),
    Cluster_Date = as.Date(as.character(Cluster_Date), format = "%Y%m%d")) %>%
  separate_rows(HiSeq_Run, sep = ",") %>%
  mutate(HiSeq_Run = paste0("HiSeq", HiSeq_Run))

cd2 <- cd2 %>%
  mutate(
    DNA_Extraction_Date = as.Date(
      as.character(DNA_Extraction_Date),
      format = "%y%m%d"),
    Sort_Date = as.Date(as.character(Sort_Date), format = "%y%m%d"),
    Tissue_Remaining = case_when(
      .$Tissue_Remaining == "y" ~ TRUE,
      .$Tissue_Remaining == "n" ~ FALSE,
      TRUE ~ NA),
    Cluster_Date = as.Date(Cluster_Date, format = "%d-%b-%y"),
    Library_Prep_Date = as.Date(Library_Prep_Date, format = "%d-%b-%y"))

# Move GTEX-U3ZN-0626-SM-9VXYL from Phase2 -------------------------------------

# NOTE: GTEX-U3ZN-0626-SM-9VXYL is a Phase 1 sample that was resequenced in
#       Phase 2 to obtain deeper coverage (see
#       https://github.com/feinberglabepigenetics/eGTEx/issues/4)

by <- setdiff(
  x = intersect(colnames(cd1), colnames(cd2)),
  y = c("HiSeq_Run", "Library_Prep_Batch"))
cd1 <- left_join(
  x = cd1,
  y = cd2,
  by = by,
  suffix = c(".Phase1", ".Phase2")) %>%
  group_by(Sample_ID_for_Data_Sharing_and_Public_Release) %>%
  mutate(
    HiSeq_Run = map2(
      .x = HiSeq_Run.Phase1,
      .y = HiSeq_Run.Phase2,
      ~ union(.x, .y)),
    Library_Prep_Batch = map2(
      .x = Library_Prep_Batch.Phase1,
      .y = Library_Prep_Batch.Phase2,
      ~ union(.x, .y))) %>%
  mutate(
    HiSeq_Run = map(HiSeq_Run, function(x) x[!is.na(x)]),
    Library_Prep_Batch = map(Library_Prep_Batch, function(x) x[!is.na(x)])) %>%
  select(
    -HiSeq_Run.Phase1, -HiSeq_Run.Phase2, -Library_Prep_Batch.Phase1,
    -Library_Prep_Batch.Phase2) %>%
  unnest(HiSeq_Run)

cd2 <- anti_join(cd2, cd1, by = by)

# Add colors -------------------------------------------------------------------

# Official GTEx colors
gtex_colours <- fromJSON(
  txt = here("metadata", "data", "GTEx_Official_Colors.json"),
  simplifyVector = FALSE) %>%
  map_df(~ data_frame(
    Tissue = .x[["tissue_name"]],
    GTEx_Color = paste0("#", .x[["tissue_color_hex"]])))
cd1 <- inner_join(cd1, gtex_colours)
cd2 <- inner_join(cd2, gtex_colours)

# Our brain region colors
brain_colors <- read_csv(here("metadata", "data", "brain_colors.csv"))
cd1 <- left_join(cd1, brain_colors) %>%
  rename(Brain_color = Color)
cd2 <- left_join(cd2, brain_colors) %>%
  rename(Brain_color = Color)

# Reorder columns --------------------------------------------------------------

col_order <- sort(colnames(cd1))
cd1 <- select(cd1, col_order)
cd2 <- select(cd2, col_order)

# Sanity check -----------------------------------------------------------------

stopifnot(identical(colnames(cd1), colnames(cd2)))

# Convert to DataFrame ---------------------------------------------------------

CD1 <- as(cd1, "DataFrame")
rownames(CD1) <- CD1[["Sample_ID_for_Data_Sharing_and_Public_Release"]]
CD2 <- as(cd2, "DataFrame")
rownames(CD2) <- CD2[["Sample_ID_for_Data_Sharing_and_Public_Release"]]

# Save colData -----------------------------------------------------------------

saveRDS(CD1, here("metadata", "objects", "colData.eGTEx.Phase1.rds"))
saveRDS(CD2, here("metadata", "objects", "colData.eGTEx.Phase2.rds"))
