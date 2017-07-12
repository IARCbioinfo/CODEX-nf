#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
projectname = args[1]
optKallchr  = args[2]
mode        = args[3]
T_suffix    = args[4]

finalcall=read.table(paste(projectname,'_results_chr1_K', optKallchr,'_',mode,'.txt', sep=''),head=T,sep='\t')

for(chr in c(2:22,'X','Y')){
  chr=paste('chr',chr, sep='')
  cat('chr =',chr,'\n')
  finalcall.chr=read.table(paste(projectname,'_results_',chr,'_K', optKallchr,'_',mode,'.txt', sep=''),head=T,sep='\t')
  finalcall=rbind(finalcall,finalcall.chr)
}

write.table(finalcall,file=paste(projectname,'_results_allchr_K', optKallchr,'_',mode,'.txt', sep=''),sep='\t',quote=F,col.names=T,row.names=F)

tumorsampname=as.matrix(unique(finalcall[,1]))

SEG = c()
for(i in 1:nrow(tumorsampname)){
  cat(i,'\t')
  finalcall.temp=finalcall[which(finalcall$sample_name==tumorsampname[i,1]),]
  sampname.temp=paste(rep(tumorsampname[i,1],nrow(finalcall.temp)),'_codex',sep='')
  output=cbind(sampname.temp,finalcall.temp$chr, finalcall.temp$st_bp,finalcall.temp$ed_bp, finalcall.temp$ed_exon-finalcall.temp$st_exon+1, signif(pmax(log(finalcall.temp$copy_no/2,2),-4),4))
  SEG = rbind(SEG,output)
}

SEG = as.data.frame(SEG)
colnames(SEG)=c('Sample', 'Chromosome','Start','End','Num_Probes', 'Segment_Mean')

write.table(SEG,file=paste('all_samples_igv.codex.seg.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)

# select T and N
SEG = SEG[grep(T_suffix,SEG$Sample),]

write.table(SEG,file=paste('T_igv.codex.seg.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)

markerfile_1             = SEG[,c(2:3)]
markerfile_2             = SEG[,c(2,4)]
colnames(markerfile_1)[2]= "Position"
colnames(markerfile_2)[2]= "Position"
markerfile               = rbind(markerfile_1,markerfile_2)
markerfile               = markerfile[!(duplicated(markerfile)),]
markerfile$V3            = c(1:nrow(markerfile))
markerfile$ID            = paste("ID_",markerfile$V3,sep="")
markerfile               = markerfile[,c(4,1:2)]
 
write.table(markerfile,file = paste('MarkerFile.txt',sep=''),sep='\t',quote=F,col.names=T,row.names=F)

 
