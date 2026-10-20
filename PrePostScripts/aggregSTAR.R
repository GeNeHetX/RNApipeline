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
if (length(args) < 1) stop("Usage: aggregSTAR.R results_directory [threads]")
datdir=args[1]
THREADS=if (length(args) >= 2) as.integer(args[2]) else 1L
if (is.na(THREADS) || THREADS < 1) stop("threads must be a positive integer")
outdir=datdir

rezdirs= file.path(datdir, c("FeatureCounts_output", "StarLog_output"))
rezdirs=rezdirs[dir.exists(rezdirs)]
allfiles=unlist(lapply(rezdirs, function(d) list.files(d, recursive = FALSE, full.names = TRUE)))

allexoncounts=grep("exonscount.txt$", allfiles, value = TRUE)
allgenecounts=grep("genecount.txt$", allfiles, value = TRUE)
allstarlog=grep("StarOutLog.final.out$", allfiles, value = TRUE)

a=basename(sub("exonscount.txt$","",allexoncounts))
b=basename(sub("genecount.txt$","",allgenecounts))
d=basename(sub("StarOutLog.final.out$","",allstarlog))
if (!length(a) || !length(b) || !length(d)) stop("Expected FeatureCounts and STAR outputs were not found")
if (anyDuplicated(a) || anyDuplicated(b) || anyDuplicated(d)) stop("Duplicate sample outputs were found")
if (!setequal(a,b) || !setequal(b,d)) stop("FeatureCounts and STAR outputs contain different samples")

names(allexoncounts)=a
names(allgenecounts)=b
names(allstarlog)=d
allsampl=b

# pattern="_R1.fastq.gz|_R1_001.fastq.gz"

library(parallel)



read_count_vectors <- function(files, sample_ids) {
  counts <- mclapply(files[sample_ids], function(f) {
    x <- read.delim(f, sep="\t", header=TRUE, as.is=TRUE, skip=1, check.names=FALSE)
    if (ncol(x) < 2 || nrow(x) == 0) stop(paste("Invalid FeatureCounts file:", f))
    ids <- x[[1]]
    if (anyDuplicated(ids)) stop(paste("Duplicate feature IDs in:", f))
    setNames(x[[ncol(x)]], ids)
  }, mc.cores=THREADS)
  names(counts) <- sample_ids
  feature_ids <- names(counts[[1]])
  if (any(!vapply(counts, function(x) identical(names(x), feature_ids), logical(1))))
    stop("Feature IDs or order differ between FeatureCounts files")
  counts
}

allFCGcountsL=read_count_vectors(allgenecounts, allsampl)
allFCEcountsL=read_count_vectors(allexoncounts, allsampl)

allStarLog=mclapply(allstarlog[allsampl],function(f){
  x=readLines(f)
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
n=length(allFCGcountsL)

print("Saving FC")
g=names(allFCGcountsL[[1]])
Y=data.frame(do.call(cbind,lapply(allFCGcountsL, identity)), check.names=FALSE)
dimnames(Y)=list(g,allsampl)
# saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
gene_file=file.path(outdir,paste0("StarFCGeneCount_s",n,".tsv"))
write.table(Y, file=gene_file, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
if (system2("gzip", gene_file) != 0) stop("Failed to gzip gene-count matrix")


g=names(allFCEcountsL[[1]])
Y=data.frame(do.call(cbind,lapply(allFCEcountsL, identity)), check.names=FALSE)
dimnames(Y)=list(g,allsampl)
# saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
exon_file=file.path(outdir,paste0("StarFCExonCount_s",n,".tsv"))
write.table(Y, file=exon_file, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
if (system2("gzip", exon_file) != 0) stop("Failed to gzip exon-count matrix")





# saveRDS(Y,file=file.path(outdir,paste0(prefix,sub(".tsv.gz","",fcsuffixes)[i],".rds")))
write.table( do.call(rbind,allStarLog), file=file.path(outdir, paste0("StarMetrics_s", n, '.tsv')),quote=F,row.names=T,col.names=T, sep="\t")


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
