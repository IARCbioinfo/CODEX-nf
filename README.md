# CODEX-nf
Pipeline for copy number variant calling from Whole Exome sequencing data using CODEX

## Prerequisites
- *R* with package *CODEX*
- *Rscript*

and the following files:
- a bed file with the positions to consider
- a folder with Normal BAM files
- a folder with Tumor BAM files

## Usage
```bash
nextflow run iarcbioinfo/CODEX-nf --dirNormal Normal --dirTumor Tumor --bedFile positions.bed --outdir output
```

## Output
A series of **files *projectname*\_*chromosome*\_*optK*\_*mode*.txt**, where *chromosome* takes its values in [1,2,...,22,X,Y], while *projectname*, *optK*, and *mode* are the same for all files. *projectname* and *mode* are specified by the user with options *--project* and *--mode*, respectively, and *optK* is computed by package CODEX.

## All parameters
| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*--help*         |null | Print usage and parameters
*--dirNormal*    |Normal | Path to Normal BAM files
*--dirTumor*     |Tumor | Path to Tumor BAM files
*--outdir*       |. | Path to output directiry
*--bedFile*      |positions.bed | Path to bed file with positions to consider
*--rem\_from\_bed* |\_random\|chrUn\|GL000209R\|GL000191R\|GL000194R | Strings to exclude from the bed file chromosome list
*--project*      |  | Project name
*--mem*          |5 | Memory requested
*--cpus*         | 1 | CPUs requested
*--seg_mode*     | fraction | Mode for the segmentation algorithm (*fraction* or *integer*)
