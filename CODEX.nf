#! /usr/bin/env nextflow
// usage : ./CODEX.nf --input_folder input/ --cpu 8 --mem 32 
/*
vim: syntax=groovy
-*- mode: groovy;-*- */

// requirement:
// - Rscript

//default values
params.dirNormal    = "Normal"
params.dirTumor     = "Tumor"
params.outdir       = "."
params.bedFile      = "positions.bed"
params.rem_from_bed = "_random|chrUn|GL000209R|GL000191R|GL000194R"
params.project      = ""

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW COPY NUMBER VARIATION CALLING FROM WHOLE EXOME PIPELINE'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run CODEX-nf --input_folder input/ [--cpu 8] [--mem 32] [--RG "PL:ILLUMINA"] [--fastq_ext fastq.gz] [--suffix1 _1] [--suffix2 _2] [--out_folder output/]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '    --fasta_ref          FILE                    Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --indel_realignment                    Performs local indel realignment (default: no).'
    log.info ''
    exit 1
}

//read files
fasta_ref = file(params.fasta_ref)
fasta_ref_fai = file( params.fasta_ref+'.fai' )
chrs = Channel.from( 'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY' )

process CODEX_normalize_perchr {
        cpus params.cpu
        memory params.mem+'G'
        tag { file_tag }
        
        input:
	env chr from chrs
     
        output:
	file("*.txt") into chr_files
	
        shell:
        file_tag = chr
        '''
	Rscript !{baseDir}/bin/codex_run.R !{params.dirNormal} !{params.dirTumor} !{params.outdir} !{params.bedFile} !{params.rem_from_bed} !{params.project} $chr
        '''
}

process CODEX_segmentation_allchr {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }
        
    input:
    chr from chr_files.toList()
    
    output:
    file("*.txt") into outfiles
    publishDir params.out_folder, mode: 'move'

    shell:
    '''
    
    '''
}