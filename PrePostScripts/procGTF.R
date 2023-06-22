#!/usr/bin/env Rscript


# args[1] "Homo_sapiens.GRCh38.101.GeneLvlOnly.gtf"
# args[2] "Homo_sapiens.GRCh38.91.all.GENES"
args = commandArgs(trailingOnly=TRUE)

library(parallel)
gtf=read.delim(args[1],sep="\t",as.is=T,header=F,comment.char="#")

gtf.getmeta=function(apieceofgtf){strsplit((apieceofgtf)[,9],";| ")}
gtf.getmetavalue=function(apieceofgtf,field){
   metadata=gtf.getmeta(apieceofgtf)
   unlist(mclapply(metadata,function(x){
       if(field%in% x){return(x[which(x==field)+1])}else{
           return(NA)
       }
   }))
}
geneTab =gtf[which(gtf[,3] =="gene"),]
geneTab$GeneID=gtf.getmetavalue(geneTab,"gene_id")
geneTab$GeneName=gtf.getmetavalue(geneTab,"gene_name")
inoname=which(is.na(geneTab$GeneName))
geneTab$GeneName[inoname]=geneTab$GeneID[inoname]
geneTab$biotype=gtf.getmetavalue(geneTab,"gene_biotype")
# geneTab$source=gtf.getmetavalue(geneTab,"gene_source")

geneTab=unique(geneTab)
rownames(geneTab)=geneTab$GeneID

geneTab=unique(geneTab[,- c(9,6,2,3,8)])
colnames(geneTab)[1:4]=c("seqname","start","end", "strand")

saveRDS(geneTab,file=paste0(args[2],".rds"))
write.table(geneTab,file=paste0(args[2],".tsv"),quote=F,sep="\t")
