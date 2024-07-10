#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#  Rscript 03_Metadata/aggreg.R 04_Results/lexofwd_a 01_Raw_files/lexo_slist &



#
#  setwd("/mnt/citprojects/CIT_METAPROJECTS/PANCREAS/PDACpatients/PilotData")

# rezdir="04_Results/RNA_access_STAR_toler"
# initSlist="01_Raw_files/RNA_access_slist"
# outdir="04_Results/allCounts2"
# prefix= "RNA_access_STAR_toler"
# initSlist=args[2]

# rezdir=args[1]
# outdir=args[3]
# prefix=args[4]
THREADS=12

# datdir="/Volumes/genomic/RNA/BjnBiliaryCancer_Beaufrere/PEsort_251022/05_Process/RNAv1.5_Ensemblv107"
# refgtf="/Volumes/genomic/REFDATA/ensembl_v107_GRCh38/geneInfo.tab"


datdir=args[1]
refgtf=args[2]
outdir=datdir

rezdir=file.path(datdir,"FeatureCounts_output")

allfiles=list.files(rezdir,recursive=F)
allexoncounts=grep("exonscount.txt$",allfiles,value=T)
allgenecounts=grep("genecount.txt$",allfiles,value=T)
allstarlog=grep("StarOutLog.final.out$",allfiles,value=T)



a=sub("exonscount.txt$","",allexoncounts)
b=sub("genecount.txt$","",allgenecounts)
d=sub("StarOutLog.final.out$","",allstarlog)
if(!setequal(a,b) |!setequal(b,d) ){stop("OMG NOT SAME SAMPLES")}

names(allexoncounts)=a
names(allgenecounts)=b
names(allstarlog)=d
allsampl=b

# pattern="_R1.fastq.gz|_R1_001.fastq.gz"

library(parallel)



allFCGcountsL=mclapply(allgenecounts[allsampl],function(f){
  x=NULL
  try({
    x=read.delim(file.path(rezdir,f),sep="\t",header=T,as.is=T,skip=1)
    rownames(x)=x[,1]
  })
  x[ncol(x)]
},mc.cores=THREADS)

allFCEcountsL=mclapply(allexoncounts[allsampl],function(f){
  x=NULL
  try({
    x=unique(read.delim(file.path(rezdir,f),sep="\t",header=T,as.is=T,skip=1))
    rownames(x)=x[,1]
  })
  x[ncol(x)]
},mc.cores=THREADS)

allStarLog=mclapply(allstarlog[allsampl],function(f){
  x=readLines(file.path(rezdir,f))
  x=x[which(x!="")]
  x=gsub("%|\t","",gsub(" ","",grep(":$",gsub("^[ \t]+","",x),value=T,invert=T)))
  x=gsub(":|\\/",".",x)
  do.call(cbind,lapply(strsplit(x,"\\|"),function(y){
    setNames(data.frame(y[2]),y[1])
  }))
})


#
# missingFQ=unique(c(
#   # which(sapply(CntLlo,nrow)==0),
#   # which(sapply(CntLO,nrow)==0),
#   unlist(lapply(allFCcounts,function(a){which(sapply(a,nrow)==0)})),
#   which(sapply(STARcountL,nrow)==0)
# ))
# print(paste("Number of missing sampels : ",length(missingFQ)))
# if(length(missingFQ)==0){print("allgood")}


print("Saving FC")
g=rownames(allFCGcountsL[[1]])
Y=data.frame(do.call(cbind,mclapply(allFCGcountsL,function(x)x[g,1])))
dimnames(Y)=list(g,allsampl)
# saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
write.table( Y, file=file.path(outdir,"StarFCGeneCount.tsv"),quote=F,row.names=T,col.names=T, sep="\t")
system(paste0("gzip ",file.path(outdir,"StarFCGeneCount.tsv")))


g=rownames(allFCEcountsL[[1]])
Y=data.frame(do.call(cbind,mclapply(allFCEcountsL,function(x)x[g,1])))
dimnames(Y)=list(g,allsampl)
# saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
write.table( Y, file=file.path(outdir,"StarFCExonCount.tsv"),quote=F,row.names=T,col.names=T, sep="\t")
system(paste0("gzip ",file.path(outdir,"StarFCExonCount.tsv")))






# saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
write.table( do.call(rbind,allStarLog), file=file.path(outdir,"StarMetrics.tsv"),quote=F,row.names=T,col.names=T, sep="\t")


# FCountLO=data.frame(do.call(cbind,mclapply(CntLlo,function(x)x[g,2])))
# dimnames(FCountLO)=list(g,allsampl)
#
#
# FCountAO=data.frame(do.call(cbind,mclapply(CntLO,function(x)x[g,2])))
# dimnames(FCountAO)=list(g,allsampl)
#


# sed "s/ | */;/" StarOutLog.final.out | sed "s/^  *//" |sed "s/ /_/g" | sed "s/^%/Prop/" \
# | sed "s/:\(\)//g" |sed "s/%$//" |grep ";" |grep -v "^Start\|^Fini\|^Mapping_\|chimer" |\
# tr -d ";" >$OUTdir"/"$ID"summary.txt"
#
# print("Exatrcting STAR summary")
# STARcountsum=data.frame(t(do.call(cbind,mclapply(STARcountL,function(x){
#   data.frame(N=c(x[1:3,2],x[1:3,3],x[1:3,4]),row.names=paste(x[1:3,1],rep(1:3,each=3),sep="_"))
# }))))
# dimnames(STARcountsum)=list(allsampl,paste0("STAR",colnames(STARcountsum)))
#
# print("Merging STAR summary")
# SumSum=data.frame(do.call(rbind,mclapply(sums,function(x){
#   x[,2]
# })),stringsAsFactors=F)
# for(j in 1:(ncol(SumSum)-1)){SumSum[,j]=as.numeric(SumSum[,j])}
# dimnames(SumSum)=list(allsampl,rownames(sums[[1]]))
#
# SumSum=data.frame(SumSum,STARcountsum[rownames(SumSum),])
