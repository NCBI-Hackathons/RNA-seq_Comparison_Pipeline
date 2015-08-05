# RNA-seq Comparison Tool

__RNA-seq Comparison Tool__ is a small collection of scripts that a user with limited experience with a terminal/command line environment could use to easily compare RNA-seq data from several tissues or conditions. These scripts are meant to be downloaded, edited and forked, and users are welcome to pull a subset or all of them. The workflow is shown below, but a user would be best served by following the short tutorial.

## Workflow
__Step 0__ 
Collect and parse metadata and SRP/ERP number. Use SRAdb or fetch9606+ to go through metadata and find which studies are of interest.

__Step 1__
Enter tissue and SRP/ERP number. fetchAndRun+ (bash) performs this.

__Step 2__
Assemble data. This is returned by fetchAndRun+ if transcript data are desired. If gene-wise data are desired GetTargetID+ (bash) or collapseGene+ (Python) collapse transcript-specific counts onto their respective genes. collapseGene+ offers a way to examine the data interactively in Python.

__Step 3__
Analyses can be run in R or with the RFCA+ package. In R, bionconductoR packages edgeR and DESeq2 provide powerful statistical tools in cases where adequate sample sizes have been attained. 

## Dependencies
__Typical__ linux (i.e.: ubuntu, KDE) should have most of the software necessary to run these scripts. Your system should have bash, Perl, Python (>= 2.7), and R installed (analysis only).

__R packages__: DESeq2, edgeR, SRAdb (and all of their respective dependencies).


## Tutorial


