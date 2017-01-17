################################################################
####          CNV AC using CODEX all chromosome 
####                  04/01/2017
####                Noémie Leblay
################################################################

args <- commandArgs(trailingOnly=TRUE)
dirNormal    = args[1]
dirTumor     = args[2]
bedFile      = args[3]
rem_from_bed = args[4]
project      = args[5]
cur_chr      = args[6]

#get bam files paths
bamFile_Normal <- list.files(dirNormal, pattern = '*.bam$')
bamFile_Normal<- bamFile_Normal

bamdir_Normal <- file.path(dirNormal, bamFile_Normal)
bamdir_Normal<- bamdir_Normal

bamFile_Tumor <- list.files(dirTumor, pattern = '*.bam$')
bamdir_Tumor <- file.path(dirTumor, bamFile_Tumor)

bamFile <- c(bamFile_Normal, bamFile_Tumor )
bamdir <- c(bamdir_Normal, bamdir_Tumor )

#get sample names
sampname <- as.matrix(rep(0,length(bamFile)),length(bamFile),1)
for (i in 1:length(bamFile)){
  sampname[i,]=strsplit(bamFile,".bam")[[i]][1]
}

sampnameN <- as.matrix(rep(0,length(bamFile_Normal)),length(bamFile_Normal),1)
for (i in 1:length(bamFile_Normal)){
  sampnameN[i,]=strsplit(bamFile_Normal,".bam")[[i]][1]
}

#read bed file
bed<- read.table(bedFile, sep = '\t')
chr_list = unique(bed[,1])
chr_list = chr_list[grep(rem_from_bed,chr_list, invert=TRUE)]

library(CODEX)
 
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, sampname = sampname, projectname = project, cur_chr)
  
bamdir <- bambedObj$bamdir
sampname <- bambedObj$sampname
ref<- bambedObj$ref
projectname <- bambedObj$projectname
chr <- bambedObj$chr

coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
readlength<- coverageObj$readlength

gc <- getgc(chr, ref)
mapp <- getmapp(chr, ref)

qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 4000),length_thresh = c(20, 2000), mapp_thresh = 0.9, gc_thresh = c(20, 80))
sampname_qc <- qcObj$sampname_qc
qcmat <- qcObj$qcmat
normal_index = which(sapply(1:length(sampname_qc), function(i) (sampname_qc[i] %in% sampnameN) ) )
    
norm.no.reads = which(apply((qcObj$Y_qc)[, normal_index], 1, median)==0)
if( length(norm.no.reads)>0 ){#in case some exons have a low read count
    Y_qc<- qcObj$Y_qc[-norm.no.reads,]
    gc_qc <- qcObj$gc_qc[-norm.no.reads]
    mapp_qc <- qcObj$mapp_qc[-norm.no.reads]
    ref_qc <- qcObj$ref_qc[-norm.no.reads,]
}else{#otherwise
    Y_qc<- qcObj$Y_qc
    gc_qc <- qcObj$gc_qc
    mapp_qc <- qcObj$mapp_qc
    ref_qc <- qcObj$ref_qc
}
write.table(qcmat, file = paste(projectname,'_',cur_chr,'_qcmat','.txt', sep=''),sep='\t', quote=FALSE, row.names=FALSE)

normObj <- normalize2(Y_qc, gc_qc, K = 1:10, normal_index)
Yhat  <- normObj$Yhat
AIC  <- normObj$AIC
BIC  <- normObj$BIC
RSS  <- normObj$RSS
K  <- normObj$K

choiceofK(AIC, BIC, RSS, K, filename = paste(projectname, "_",cur_chr,"_choiceofK", ".pdf", sep = ""))
optK = K[which.max(BIC)]

# write files
write.table(Y_qc, file = paste(projectname,'_',cur_chr,'_Y_qc','.txt', sep=''),sep='\t', quote=FALSE, row.names=FALSE)
write.table(Yhat, file = paste(projectname,'_',cur_chr,'_Yhat','.txt', sep=''),sep='\t', quote=FALSE, row.names=FALSE)
write(optK,file        = paste(projectname,'_',cur_chr,'_optK','.txt', sep='') )
write.table(ref_qc,file = paste(projectname,'_',cur_chr,'_ref_qc','.txt', sep=''),sep='\t', quote=FALSE, row.names=FALSE)
write.table(sampname_qc,file = paste(projectname,'_',cur_chr,'_sampname_qc','.txt', sep=''),sep='\t', quote=FALSE, row.names=FALSE)


