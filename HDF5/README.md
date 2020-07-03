HDF5Array-backed BSseq objects
================
Edited: 20 September, 2018; Compiled: 20 September 2018

# Important tips for running on JHPCE

1.  If requesting an interactive session via `qrsh`, be sure to use
    `qrsh bash`, e.g., `qrsh -l mem_free=80G,h_vmem=80.1G -pe local 4
    bash`. This ensures a large-enough `TMPDIR` environment variable is
    set. This is not required if submitting a job via `qsub`.
2.  Use `module load conda_R/3.5.x`
3.  Copy the directory containing the HDF5-backed *BSseq* object to
    scratch disk and use this copy in your analysis (loaded with
    `HDF5Array::loadHDF5SummarizedExperiment()`). This may improve
    performance if you are doing lots of reading from (or writing to)
    the HDF5 file. This step shouldn’t by necessary if you’re doing
    something small (e.g., plotting a few DMRs).
4.  Include this at the top of your R script or at the beginning of your
    R session.

<!-- end list -->

``` r
library(bsseq)
library(HDF5Array) # Needed for loadHDF5SummarizedExperiment()
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE) # Needed for HDF5 file stored on some JHPCE file systems
```

The final tip is a workaround to ensure you can load HDF5 file stored on
`/dcl01/FB2/data/personal/gtex/` (see
<https://github.com/grimbough/Rhdf5lib/issues/11> for technical
discussion). If you first copy the HDF5-backed *BSseq* object to
`/scratch` then you don’t need to do this step, but it shouldn’t hurt.
