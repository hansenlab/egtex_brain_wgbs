eGTEx sample metadata
================
Edited: 15 August, 2018; Compiled: 15 August 2018

Metadata for the eGTEx project, including the [Phase 1 sample
spreadsheet (CSV)](extdata/GTEx_pData.csv), [Phase 2 sample spreadsheet
(CSV)](extdata/GTEx_Flow_Sorted_pData.csv), and the [GTEx official
tissue colours(JSON)](extdata/GTEx_Official_Colors.json).

The `objects/` directory contains *S4Vectors::DataFrame* objects
constructed from the sample spreadsheets.

# Duplicate samples

`GTEX-U3ZN-0626-SM-9VXYL` is a Phase1 sample
(<https://github.com/feinberglabepigenetics/eGTEx/blob/master/metadata/extdata/GTEx_pData.csv#L61>)
that was resequenced (from the same library) in Phase2
(<https://github.com/feinberglabepigenetics/eGTEx/blob/master/metadata/extdata/GTEx_Flow_Sorted_pData.csv#L2>).
Since it is considered a Phase 1 sample, it is excluded from the Phase 2
object.

# Adding dbGaP-protected metadata to an eGTEx object

Some of the metadata for these samples is protected and only available
via dbGaP. **These data must not be added to the GitHub repository**. If
you require dbGaP-protected metadata in your analysis, you will need to
first add it to the relevant object (e.g., the *BSseq* object). Below is
an example of how to do this.

``` r
library(bsseq)
library(HDF5Array)
library(dplyr)

# Load the object you need for analysis ----------------------------------------

bsseq <- loadHDF5SummarizedExperiment("/path/to/object")

# Load the GTEx dbGaP-protected phenotype data ---------------------------------

# NOTE: You will only be able to access this file if you have the appropriate 
#       permissions.
load("/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/GTEx_pheno_data.rda")

# Add the dbGaP-protected metadata to the BSseq object -------------------------

# NOTE: Need to take care to sensible combine with existing colData.
colData(bsseq) <- left_join(
  x = as.data.frame(colData(bsseq)),
  y = GTEx_pheno_data,
  by = c(
    "Collaborator_Participant_ID" = "Individual_ID",
    "Tissue" = "SMTSD")) %>%
  # NOTE: Remove columns that are specific to the RNA-seq assay in 
  #       'GTEx_pheno_data'.
  select(-starts_with("SM", ignore.case = FALSE), -ANALYTE_TYPE, -SAMPID) %>%
  DataFrame(row.names = .$Sample_ID_for_Data_Sharing_and_Public_Release)
```
