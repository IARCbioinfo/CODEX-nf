#!/usr/bin/env Rscript

################################################################
####          Annotation CNV
####                Noémie Leblay
################################################################

args <- commandArgs(trailingOnly=TRUE)
GTF     = args[1]
segfile = args[2]
projectname = args[3]

gene2              = read.table(GTF,quote="\"",stringsAsFactors=F,sep="\t",header=F)
colnames(gene2)[1] = "Chr"
colnames(gene2)[3] = "info"
colnames(gene2)[4] = "Start"
colnames(gene2)[5] = "End"

library (reshape2)
info2      = colsplit(gene2$V9, ";", c(1:10))
gene_info2 = cbind(gene2[,-9], info2[,])
gene_info2 = gene_info2[gene_info2$info=="gene",]
colnames(gene_info2)[9]  = "gene_id"
colnames(gene_info2)[10] = "gene_name"
colnames(gene_info2)[11] = "gene_source"
colnames(gene_info2)[2]  = "gene_biotype"
gene_info2               = gene_info2[,1:11]

gene_info2[,"gene_name"]=gsub("^\\s+|\\s+$", "", gene_info2[,10])
gene_info2$gene_name=unlist(lapply(gene_info2$gene_name,function(x){strsplit(x, " ")[[1]][2]}))

#download the CODEX results files and take a single results file with all the chr

CODEX=read.table(segfile,quote="\"",stringsAsFactors=F,sep="\t",header=T)

#Annotation start and end gene involved in the CNV
CODEX$Start_gene=unlist(lapply(1:length(CODEX$Start), function(x){res=gene_info2[gene_info2$Start<= CODEX[x, "st_bp"] & gene_info2$End>=CODEX[x, "st_bp"] & gene_info2$Chr==CODEX[x, "chr"],"gene_name"] if( length(res)==0){NA} else{paste(res, collapse=",")}}))


CODEX$End_gene=unlist(lapply(1:length(CODEX$End), function(x){res=gene_info2[gene_info2$Start<= CODEX[x, "ed_bp"] &
                                                                                 gene_info2$End>=CODEX[x, "ed_bp"] &
                                                                                 gene_info2$Chr==CODEX[x, "chr"],"gene_name"]
                                                                if( length(res)==0){NA} else{paste(res, collapse=",")}}))

#Annotation genes involved in the CNV

CODEX$Involved_genes=unlist(lapply(1:length(CODEX$Start), function(x){res=gene_info2[gene_info2$Start>= CODEX[x, "st_bp"] &
                                                                                       gene_info2$Start<=CODEX[x, "ed_bp"] &
                                                                                       gene_info2$End<=CODEX[x, "ed_bp"] &
                                                                                       gene_info2$End>=CODEX[x, "st_bp"] &
                                                                                       gene_info2$Chr==CODEX[x, "chr"],"gene_name"]
                                                                      if( length(res)==0){NA} else{paste(res, collapse=",")}}))

#Annotation exon involved in the CNV for the starting and ending genes

CODEX$Percentage_Start_gene=unlist(lapply(1:length(CODEX$Start), function(x){ if(is.na(CODEX[x,"Start_gene"])) { return(NA)}
                                                                              print(x)
                                                                              res_start=gene_info2[gene_info2$Start<= CODEX[x, "st_bp"] &
                                                                                                     gene_info2$End>=CODEX[x, "st_bp"] &
                                                                                                     gene_info2$Chr==CODEX[x, "chr"],"Start"]
                                                                              res_end=gene_info2[gene_info2$Start<= CODEX[x, "st_bp"] &
                                                                                                   gene_info2$End>=CODEX[x, "st_bp"] &
                                                                                                   gene_info2$Chr==CODEX[x, "chr"],"End"]
                                                                              res=unlist(lapply(1: length(res_start), function(i){res_int=intersect(seq(res_start[i], res_end[i], 1), seq (CODEX[x,"st_bp"], CODEX[x, "ed_bp"],1))
                                                                                                                                  length(res_int)/length(seq(res_start[i], res_end[i], 1))*100
                                                                              }))
                                                                              if( length(res)==0){NA} else{paste(res, collapse=",")}}))

CODEX$Percentage_End_gene=unlist(lapply(1:length(CODEX$End), function(x){ if(is.na(CODEX[x,"End_gene"])) { return(NA)}
                                                                            print(x)
                                                                            res_start=gene_info2[gene_info2$Start<= CODEX[x, "ed_bp"] &
                                                                                                   gene_info2$End>=CODEX[x, "ed_bp"] &
                                                                                                   gene_info2$Chr==CODEX[x, "chr"],"Start"]
                                                                            res_end=gene_info2[gene_info2$Start<= CODEX[x, "ed_bp"] &
                                                                                                 gene_info2$End>=CODEX[x, "ed_bp"] &
                                                                                                 gene_info2$Chr==CODEX[x, "chr"],"End"]
                                                                            res=unlist(lapply(1: length(res_start), function(i){res_int=intersect(seq(res_start[i], res_end[i], 1), seq (CODEX[x,"st_bp"], CODEX[x, "ed_bp"],1))
                                                                                                                                length(res_int)/length(seq(res_start[i], res_end[i], 1))*100
                                                                            }))
                                                                            if( length(res)==0){NA} else{paste(res, collapse=",")}}))

write.table(CODEX, file = paste(projectname,"allChr_annotated.txt",sep="_"), sep='\t', quote=FALSE, row.names=FALSE)

