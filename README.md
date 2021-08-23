# Investigating therapy-driven changes in isoform expression in glioblastoma - supporting code :dna: :brain:

This repository contains the code for my MSc project titled "Investigating therapy-driven changes in isoform expression in glioblastoma".

### Structure

* **bashscripts** - .sh scripts submitted to the HPC at the University of Leeds
* **dataexploration** - .R scripts for running PCA on in-house gene and in-house/GLASS isoform log2FC data
* **dataprocessing** - .R scripts for filtering patients and gene/transcripts in the raw in-house and GLASS datasets. Three sub-folders:
  * Sample_Filter_Logic - .R scripts for selecting patients to remove based on their Metadata and Sequencing metric values
  * Transcript_Filter_Logic - .R scripts for selecting genes/transcripts to remove based on low expression values
  * Filtered_Data - .R scripts to produce final filtered data ready for DEA
* **dea** - .R scripts for running DESeq2 and processing results. Sub-folders:
  * in-house - DEA based on in-house data (includes scripts for paired and responder-type DEA, also includes script for DEA based on sub-samples of up-responders)
  * glass - DEA based on glass data (includes script for DEA based on sub-sample of up-responders)
  * tpm_fpkm_check - per-patient bar plots of TPM (GLASS) and FPKM (in-house) expression values for specific isoforms
  * master_deseq2_isoform_results.R - script used to produce Volcano Plots
* **resultsanalysis** - .R scripts for further exploration of results including:
  * GO analysis - .R scripts for extracting lists of interesting genes and visualising overlap in GO terms between different GO enrichment analyses
  * summary.dataframe.isoforms.R - summary of all DEA results
  * isoform_switch.R - script used to identify cases of isoform switching
  * TSS.R - scripts used to map JBSgene transcripts to JARID2 TSSs, and to calculate % of JBS genes with transcripts without a JARID2 TSS


N.B. All folders titled "archive" contain scripts not directly used to produce reported results.

