#==========================================================
# README for DNA Methylation Data for the eGTEx release
#==========================================================

This document describes the output format and column headers
for bigwig, bigBED, and fastq files from the eGTEx data release.

#-----------------------------------------------
# bigWig Files
#-----------------------------------------------
  Naming convention for individual samples: GTEX-Sample-ID.context_strand.smoothing_type
  Example: GTEX-11WQC-0011-R10a-SM-AFUIX.mCT_pos.small_smooth.bw
  Naming convention for samples grouped by tissue:  TissueGroup.context_strand.smoothing_type
  Example: Cortical.mCG.small_smooth.bw
  File extension: *.bw
  Column headers:

    chr:					chromosome
    start:					start position
    score:                  mean methylation value for each sample(s)

  There are 8 bigwig files for each brain tissue sample (2 mCG files for Lung and Thyroid tissues) – 
	•	2 for mCG (large and small smoothed – corresponds to finding Blocks of differential methylation vs DMRs)
	•	6 for nonCpG methylation – positive and negative strands for mCA, mCC, and mCT all small smoothed for DMR analysis.
  There are also 8 bigwig files per tissue profiled (except for Lung/Thyroid). These are merged values from each sample for each tissue. Same description applies.
Metadata column headers:
	filename:										Name of bigwig file.
	Description:									Description of data in file.	
	Sample_ID_for_Data_Sharing_and_Public_Release:	Sample ID used across all of GTEx
	Tissue:											Tissue type.							
    Hippocampus Group:								Which group hippocampus samples were assigned - HC1 (hippocampus) and HC2 (dentate gyrus). 



#-----------------------------------------------
# bigBED Files
#-----------------------------------------------

  File extension: *.bb
  Column headers:

    chr:					chromosome
    start:					start position
    strand:                 no strand information (*)
    name:					random number
    score:					all are 0 by default; no meaning
  
  We performed 4 differential DNA methylation analyses:
	•	Overall analysis = differential analysis among NeuN+ samples from all 8 brain regions
		o	2 files = CG-DMRs and large blocks of differential mCG
	•	General analysis = differential analysis among NeuN+ samples from 5 brain tissue groups (Cortical, Basal ganglia,   
		Hippocampus, Amygdala, Hypothalamus)
		o	6 files = CH-DMRs identified for each strand and context (CA,CC,CT; pos and neg strands)
		o	2 files = CG-DMRs and large blocks of differential mCG
	•	Basal ganglia analysis = differential analysis among NeuN+ samples from Nucleus accumbens, Putamen, and Caudate.
		o	6 files = CH-DMRs identified for each strand and context (CA,CC,CT; pos and neg strands)
		o	2 files = CG-DMRs and large blocks of differential mCG
	•	Hippocampus analysis = differential analysis between NeuN+ samples from the two hippocampal groups.
		o	6 files = CH-DMRs identified for each strand and context (CA,CC,CT; pos and neg strands)
		o	2 files = CG-DMRs and large blocks of differential mCG


#-----------------------------------------------
# fastq Files (Bulk and Sorted)
#-----------------------------------------------
	File extension: *.fastq.gz

  Sorted = DNA extracted from NeuN+ nuclei after fluorsecence activated nuclei sorting.
 
  Metadata column headers:
	Read1:											Name of fastq file for read 1.
	Read2:											Name of fastq file for read 2.
	Sample_ID_for_Data_Sharing_and_Public_Release:	Sample ID used across all of GTEx
	Average_Library_Size (bp):						Average size (bp) of completed WGBS library.
	Container:										ID of container tissue or DNA was shipped in from LDACC.
	Position:										Position of tissue/DNA sample in the container.	
	Position_ColWise:								Column position of tissue/DNA sample in container.	
	HiSeq_Run:										HiSeq run ID for each library.
	Input_DNA (ng):									Amount of DNA (ng) used as input for WGBS.
	Library_Prep_Date:								Date of WGBS library preparation.
	Mass (mg):										Total mass (mg) of tissue received from LDACC as reported by LDACC.
	NeuNpos_DNA_Concentration (ng/uL):				Concentration (ng/uL) of DNA isolated from FANS NeuN positive nuclei.
	NeuNpos_Nuclei_Count:							Total number of NeuN positive nuclei sorted by FANS and used to isolate DNA.
	NeuNpos_Total_DNA (ng):							Total DNA (ng) isolated from FANS NeuN positive nuclei.		
	Participant_ID:									ID for the individual donor.	
	Library Conc.(ng/uL):							Concentration of WGBS library measured by Qubit Broad Range kit.			
	Sample_ID:										GTEx Sample ID	
	Sample_Type:									All samples were "Normal" (non-disease).	
	Self_Reported_Sex:								Sex as reported by donor.	
	Sort_Date:										Date of nuclei isolation and FANS.
	Stock_Sample:									ID of DNA stock from which our aliquot was sampled from.	
	Tissue:											Tissue type.							
    Hippocampus Group:								Which group hippocampus samples were assigned to - HC1 (hippocampus) and HC2 (dentate gyrus). 

Bulk = DNA extracted directly from tissue or received as DNA from LDACC.

Metadata column headers:
	Read1:											Name of fastq file for read 1.
	Read2:											Name of fastq file for read 2.
	Sample_ID_for_Data_Sharing_and_Public_Release:	Sample ID used across all of GTEx
	Stock_Sample:									ID of DNA stock from which our aliquot was sampled from.
	Root_Sample:									Root sample ID.
	Participant_ID:									ID for the individual donor.	
	Self_Reported_Sex:								Sex as reported by donor.	
	Tissue:											Tissue type.							
    HiSeq_Run:										HiSeq run ID for each library. Some samples were sequenced across two flowcells and are indicated.
	Average_Library_Size (bp):						Average size (bp) of completed WGBS library.
	Library_Prep_Date:								Date of WGBS library preparation.
	







	

