---
output: html_document
title: Whole Genome Bisulfite Sequencing (WGBS)
---

## Introduction

This document describes the whole genome bisulfite sequencing (WGBS) genome browser tracks for the project 'Human brain region-specific variably methylated regions (VMRs) are enriched for heritability of distinct neuropsychiatric traits'. 
For an overview of the project, please see the [README](https://s3.us-east-2.amazonaws.com/brainepigenome/hg19/docs/README.html).

## Read pre-processing, alignment, and post-processing

We trimmed reads of their adapter sequences using [**Trim Galore!** (v0.4.0)](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and quality-trimmed using the following parameters: `trim_galore -q 30 --paired ${READ1} ${READ2}`. 
We then aligned these trimmed reads to the hg38 build of the human genome (including autosomes, sex chromosomes, mitochondrial sequence, and lambda phage (accession NC_001416.1) but excluding non-chromosomal sequences) using [**Bismark** (v0.19.0)](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) with the following alignment parameters: `bismark --bam --bowtie2 -1 ${READ1} -2 ${READ2}`. 
Using the reads aligned to the lambda phage genome, we estimated that all libraries had a bisulfite conversion rate > 99%.

We then used `bismark_methylation_extractor` to summarize the number of reads supporting a methylated cytosine and the number of reads supported a unmethylated cytosine for every cytosine in the reference genome. 
Specifically, we first computed and visually inspected the M-bias of our libraries. 
Based on these results, we decided to ignore the first 2 bp and last 1 bp of read1 and the first and last 3 bp of read2 in the subsequent call to bismark_methylation_extractor with parameters: `--ignore 2 --ignore_r2 3 --ignore_3prime 1 --ignore_3prime_r2 3`. Additional parameters were: `-p --bedGraph --counts --cytosine_report --comprehensive --merge_non_CpG --multicore 4 --gzip`.
The final cytosine report file summarizes the methylation evidence at each cytosine in the reference genome.

## Smoothing WGBS

We used [**BSmooth**](https://bioconductor.org/packages/bsseq/) to estimate CpG methylation levels as previously described. 
Specifically, we ran a 'small' smooth to identify small DMRs (smoothing over windows of at least 1 kb containing at least 70 CpGs) and a 'large' smooth to identify large-scale blocks (smoothing over windows of at least 20 kb containing at least 500 CpGs). 
Following smoothing, we analyzed all CpGs that had a sequencing coverage of at least 1 in all samples.

Unlike CpGs, CpAs and CpTs are not palindromic, so were analyzed separately for each strand, for a total of 4 strand/dinucleotide combinations: 

- mCA (forward strand)
- mCA (reverse strand)
- mCT (forward strand)
- mCT (reverse strand)

For display purposes, mCA and mCT on the reverse strand are shown as negative values.

For each dinucleotide/strand combination we ran a single 'small-ish' smooth to identify DMRs (smoothing over windows of at least 3 kb containing at least 200 CpAs or CpTs). 
Following smoothing, we analyzed all CpAs and CpTs regardless of sequencing coverage.

## Identification of differentially methylated regions (DMRs)

We ran separate analyses to identify 10 types of differentially methylated regions (DMRs):

1. CG-DMRs: Using data from the 'small' smooth of CpG methylation levels
2. CG-blocks: Using data from the 'large' smooth of CpG methylation levels
3. CA-DMRs (forward strand): Using data from the 'small-ish' smooth of CpA methylation levels on the forward strand
4. CA-DMRs (reverse strand): Using data from the 'small-ish' smooth of CpA methylation levels on the reverse strand
5. CT-DMRs (forward strand): Using data from the 'small-ish' smooth of CpT methylation levels on the forward strand
6. CT-DMRs (reverse strand): Using data from the 'small-ish' smooth of CpT methylation levels on the reverse strand

Previously, we have used [**BSmooth**](https://bioconductor.org/packages/bsseq/) to perform pairwise (two-group) comparisons. 
In the present study, we had 8 brain regions to compare: Frontal cortex, Anterior cingulate cortex, Putamen, Caudate, Nucleus accumbens, Hypothalamus, Hippocampus, and Amygdala. 
We used our previously described Rather F-statistic method to perform multi-group comparisons among all 8 brain regions. Additionally, we collapsed several similar brain regions together for a 5-group comparison: Cortical (Frontal cortex, Anterior cingulate cortex), Basal ganglia (Putamen, Caudate, Nucleus accumbens), Hypothalamus, Hippocampus, and Amygdala.

For the F-statistic method, we constructed a design matrix with a term for each group. 
For each CpX (CpG, CpT, or CpA), we then fitted a linear model of the smoothed methylation levels against the design matrix. 
To improve standard error estimates, we thresholded the residual standard deviations at the 75% percentile and smoothed these using a running mean over windows containing 101 CpXs.

We then combined the estimated coefficients from the linear model, their estimated correlations, and the smoothed residual standard deviations to form F-statistics to summarize the evidence that methylation differs between the groups at each of the CpXs.

Next, we identified runs of CpXs where the F-statistic exceeded a cutoff and where each CpX was within a maximum distance of the next. 
Specifically, we used cutoffs of F = 4.6^2 for CG-DMRs, F = 2^2 for CG-blocks (following10), and F = 4^2 for CA-DMRs and CT-DMRs, and required that the CpXs were within 300 bp of one another 
for DMRs and 1000 bp of one another for blocks. 
For blocks, we also required that the average methylation in the block varied by at least 0.1 across the groups. 
These runs of CpXs formed our candidate DMRs and blocks.
Each candidate DMR and block was summarized by the area under the curve formed when treating the F-statistic as a function along the genome (`areaStat`).

We used permutation testing to assign a measure of statistical significance to each candidate DMR/block. 
We randomly permuted the design matrix, effectively permuting the sample labels, and repeated the F-statistic analysis with the same cutoffs using the permuted design matrix, resulting in a set of null DMRs/blocks for each permutation. 
We performed 1000 such permutations. 
We then asked, for each candidate DMR/block, in how many permutations did we see a null DMR/block anywhere in the genome with the same or better `areaStat` as the candidate DMR/block; dividing this number by the total number of permutations gives a permutation P-value for each DMR/block. 
Since we are comparing each candidate block/DMR against anything found anywhere in the genome in the permutation set, we are also correcting for multiple testing by controlling the family-wise error rate. 
Those candidates DMRs/blocks with a permutation P-value ≤ 0.05 form our set of DMRs/blocks.

## Annotation of differentially methylated regions

The F-statistic approach allows us to jointly use all samples for the identification of differentially methylated regions. 
However, it does not tell us which group(s) are hypomethylated or hypermethylated for the region. 
To assign such labels to our F-statistic CG-DMRs and CG-blocks, we used a post-hoc analysis for all pairwise comparisons.

We identified small DMRs and blocks using the original t-statistic method of BSmooth; an F-statistic DMR was assigned a label (e.g., hypermethylated in BA9 and hypomethylated in AMY) if the corresponding t-statistic DMR or block overlapped at least 50% of the F-statistic DMR or block. 
This procedure does not change the coordinates of the DMR/block and means an F-statistic DMR/block may be assigned multiple labels. 

No such analysis was performed for CA-DMRs or CT-DMRs due to computational costs.

## Credits

Questions should be directed to [Peter Hickey](mailto:peter.hickey@gmail.com).

