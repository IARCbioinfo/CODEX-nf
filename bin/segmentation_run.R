args <- commandArgs(trailingOnly=TRUE)
mode = args[1]

# load results from all chr

finalcall<- segmentbis(Y_qc, Yhat, optK = globaloptK, K = K, sampname_qc,ref_qc, chr, lmax = 200, mode = mode)
write.table(finalcall, file = paste(outdir,projectname,'_',cur_chr,'_', optK,'_',mode,'.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)
