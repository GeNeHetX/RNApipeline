#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#  Rscript 03_Metadata/aggreg.R 04_Results/lexofwd_a 01_Raw_files/lexo_slist &



#
#  setwd("/mnt/citprojects/CIT_METAPROJECTS/PANCREAS/PDACpatients/PilotData")

# rezdir="04_Results/RNA_access_STAR_toler"
# initSlist="01_Raw_files/RNA_access_slist"
# outdir="04_Results/allCounts2"
# prefix= "RNA_access_STAR_toler"

rezdir=args[1]

initSlist=args[2]

outdir=args[3]

prefix=args[4]




#  rezdir="04_Results/lexofwd_strict"; initSlist="01_Raw_files/lexo_slist"; outdir="04_Results/allCounts2"; prefix="lexofwd_strict"


# rezdir="04_Results/lexofwd_toler"
# initSlist="01_Raw_files/lexo_slist"
# outdir="04_Results/allCounts"
# prefix="lexofwd_toler"

# if(is.na(args[3])){
#   pattern="_R1.fastq.gz"
# }else{
#   pattern=args[3]
# }

pattern="_R1.fastq.gz|_R1_001.fastq.gz"

library(parallel)


allfq=scan(initSlist,what="character",sep="\n")
allsampl=gsub(pattern,"",allfq)

I=1:length(allfq)





fcsuffixes=c(
#"_genecounts_FCcounts_s0_Largo.tsv.gz","_genecounts_FCcounts_s2_Largo.tsv.gz","_genecounts_FCcounts_s1_Largo.tsv.gz",
  "_genecounts_FCcounts_s0_O.tsv.gz",
"_genecounts_FCcounts_s1_O.tsv.gz",
"_genecounts_FCcounts_s2_O.tsv.gz")

allFCcounts=lapply(fcsuffixes,function(suff){
  mclapply(I,function(i){
    x=NULL
    try({
      x=read.delim(file.path(rezdir,paste0(allsampl[i],suff)),
      sep="\t",header=F,skip=1)
      rownames(x)=x[,1]
    })
    x
  })
})

# CntLlo=mclapply(I,function(i){
#   x=NULL
#   try({
#     x=read.delim(file.path(rezdir,paste0(allsampl[i],"_genecounts_FCLo.tsv.gz")),
#     sep="\t",header=F,skip=1)
#     rownames(x)=x[,1]
#   })
#   x
# })
#
# CntLO=mclapply(I,function(i){
#   x=NULL
#   try({
#     x=read.delim(file.path(rezdir,paste0(allsampl[i],"_genecounts_FCO.tsv.gz")),
#     sep="\t",header=F,skip=1)
#     rownames(x)=x[,1]
#   })
#   x
# })


STARcountL=mclapply(I,function(i){
  x=NULL
  try({
    x=read.delim(file.path(rezdir,paste0(allsampl[i],"_genecounts_STAR.tsv.gz")),
    sep="\t",header=F)
    rownames(x)=x[,1]
  })
  x
})



sums=mclapply(I,function(i){
  x=NULL
  try({
    x=read.delim(file.path(rezdir,paste0(allsampl[i],"summary.txt")),
    sep="\t",header=F,as.is=T)
    rownames(x)=x[,1]
  })
  x
})


missingFQ=unique(c(
  # which(sapply(CntLlo,nrow)==0),
  # which(sapply(CntLO,nrow)==0),
  unlist(lapply(allFCcounts,function(a){which(sapply(a,nrow)==0)})),
  which(sapply(STARcountL,nrow)==0)
))
print(paste("Number of missing sampels : ",length(missingFQ)))
if(length(missingFQ)==0){print("allgood")}


print("Saving FC")
g=rownames(allFCcounts[[1]][[1]])
for(i in 1:length(allFCcounts) ){

  Y=data.frame(do.call(cbind,mclapply(allFCcounts[[i]],function(x)x[g,2])))
  dimnames(Y)=list(g,allsampl)
  saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
}
#
#
#
# FCountLO=data.frame(do.call(cbind,mclapply(CntLlo,function(x)x[g,2])))
# dimnames(FCountLO)=list(g,allsampl)
#
#
# FCountAO=data.frame(do.call(cbind,mclapply(CntLO,function(x)x[g,2])))
# dimnames(FCountAO)=list(g,allsampl)
#

print("Exatrcting STAR count")
STARcount1=data.frame(do.call(cbind,mclapply(STARcountL,function(x)x[g,2])))
dimnames(STARcount1)=list(g,allsampl)

STARcount2=data.frame(do.call(cbind,mclapply(STARcountL,function(x)x[g,3])))
dimnames(STARcount2)=list(g,allsampl)

STARcount3=data.frame(do.call(cbind,mclapply(STARcountL,function(x)x[g,4])))
dimnames(STARcount3)=list(g,allsampl)


print("Exatrcting STAR summary")
STARcountsum=data.frame(t(do.call(cbind,mclapply(STARcountL,function(x){
  data.frame(N=c(x[1:3,2],x[1:3,3],x[1:3,4]),row.names=paste(x[1:3,1],rep(1:3,each=3),sep="_"))
}))))
dimnames(STARcountsum)=list(allsampl,paste0("STAR",colnames(STARcountsum)))

print("Merging STAR summary")
SumSum=data.frame(do.call(rbind,mclapply(sums,function(x){
  x[,2]
})),stringsAsFactors=F)
for(j in 1:(ncol(SumSum)-1)){SumSum[,j]=as.numeric(SumSum[,j])}
dimnames(SumSum)=list(allsampl,rownames(sums[[1]]))

SumSum=data.frame(SumSum,STARcountsum[rownames(SumSum),])



print("Savin STAR counts")
# saveRDS(FCountAO,file=file.path(outdir,paste0(prefix,"_FCountAO.rds")))
# saveRDS(FCountLO,file=file.path(outdir,paste0(prefix,"_FCountLO.rds")))
saveRDS(STARcount1,file=file.path(outdir,paste0(prefix,"_STARcount1.rds")))
saveRDS(STARcount2,file=file.path(outdir,paste0(prefix,"_STARcount2.rds")))
saveRDS(STARcount3,file=file.path(outdir,paste0(prefix,"_STARcount3.rds")))
saveRDS(SumSum,file=file.path(outdir,paste0(prefix,"_CountSummary.rds")))

# cd /mnt/citprojects/CIT_METAPROJECTS/PANCREAS/PDACpatients/PilotData/04_Results/truseq_STARs2
# rename "s/^/Trueseq_STAR_a_/" *.rds
# cp *rds ../allCounts/

# rename "s/Trueseq_STAR_b/Trueseq_STAR_a/" *.rds
