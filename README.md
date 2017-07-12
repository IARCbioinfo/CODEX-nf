# CODEX-nf

## Copy number variation calling from whole exome sequencing data

## Description
Pipeline to detect copy number variations from Whole Exome sequencing data using CODEX, annotate the results, and produce an IGV-readable output

## Dependencies
- *R* with package *CODEX*
- *Rscript*

and the following files:
- a bed file with the positions to consider
- a folder with Normal BAM files
- a folder with Tumor BAM files

## Input 
| Type      | Description     |
|-----------|---------------|
|--dirNormal    | Path to Normal BAM files
|--dirTumor     | Path to Tumor BAM files
  

## Parameters

* #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
|--dirNormal    |Normal | Path to Normal BAM files
|--dirTumor     |Tumor | Path to Tumor BAM files
|--bedFile      |positions.bed | Path to bed file with positions to consider
  
* #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------| 
|--outdir       |. | Path to output directiry
|--rem\_from\_bed |\_random\|chrUn\|GL000209R\|GL000191R\|GL000194R | Strings to exclude from the bed file chromosome list
|--project      |  | Project name
|--mem          |5 | Memory requested
|--cpus         | 1 | CPUs requested
|--seg_mode     | fraction | Mode for the segmentation algorithm (*fraction* or *integer*)

  * #### Flags
  
Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------| 
| --help    | Display help |

	
## Usage
```bash
nextflow run iarcbioinfo/CODEX-nf --dirNormal Normal --dirTumor Tumor --bedFile positions.bed --outdir output
```
  
## Output 
| Type      | Description     |
|-----------|---------------|
| allChr\_annotated.txt | output files with segments, associated copy number, modified BIC, and annotation |
| codex.seg.txt  | IGV-readable segmentation file wih CNVs and locations|
| MarkerFile.txt | List of markers and location |

Note that files codex.seg.txt and MarkerFile.txt can be used as input of software GISTIC

## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/CODEX-nf/blob/dev/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Noemie Leblay*    |    LeblayN@students.iarc.fr | Developer to contact for support|
  | Nicolas Alcala*   |   AlcalaN@fellows.iarc.fr | DDeveloper to contact for support |
  
## References

	
## FAQ

### How much memory should I have?
The amount of memory depends mainly on the number of exons and the windows size during the segmentation step. Considering a typical bed file with most known exons and a window size *l*=200, the default 5G memory "--mem 5" seems to be enough. For *l*=5000, 50G (--mem 50) is preferable to handle the largest chromosomes.

