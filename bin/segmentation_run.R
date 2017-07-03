#!/usr/bin/env Rscript

#get arguments
library(CODEX)
args       <- commandArgs(trailingOnly=TRUE)
mode        = args[1]
Y_qc        = read.table(args[2],h=T)
Yhat        = read.table(args[3],h=T)
optKallchr  = as.numeric(args[4])
sampname_qc = scan(args[5],what="character",skip=1)
ref_qc      = read.table(args[6],h=T)
projectname = args[7]
path        = args[8]
lmax        = as.numeric(args[9])
Kmax        = ncol(Yhat)/length(sampname_qc)
K           = 1:Kmax
Yhat        = lapply(K-1, function(i) Yhat[,(1+length(sampname_qc)*i):(length(sampname_qc)*(i+1))] )
ref_qc      = IRanges(ref_qc$start,ref_qc$end,ref_qc$width)
chr         = grep("chr", strsplit(args[2],"_")[[1]],value=T)

# useful function
source(paste(path,"bin/segmentbis.R",sep=""))

# compute segmentation
finalcall<- segmentbis(Y_qc, Yhat, optK = optKallchr, K = K, sampname_qc,ref_qc, chr, lmax = lmax, mode = mode)
# write results
write.table(finalcall, file = paste(projectname,'_results_',chr,'_K', optKallchr,'_',mode,'.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)
