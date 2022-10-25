#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# R -e 'library(GenomicFeatures);tfdb=makeTxDbFromGFF("Homo_sapiens.GRCh38.101.gtf");saveDb(tfdb,"Homo_sapiens.GRCh38.101.txdb");k <- keys(tfdb, keytype = "TXNAME");tx2gene <- select(tfdb, k, "GENEID", "TXNAME");saveRDS(tx2gene,"Homo_sapiens.GRCh38.101.tx2gene.rds")'


#  cd  /mnt/citprojects/CIT_METAPROJECTS/PANCREAS/PDACpatients/PilotData
#  Rscript 03_Metadata/aggregKallisto.R 04_Results/PacaomicsCellLines_kallisto 01_Raw_files/PacaomicsCellLines_slist.txt &

# $aggkallisto

rezdir=args[1]
refgtf=args[2]


library(GenomicFeatures);tfdb=makeTxDbFromGFF(refgtf);
# saveDb(tfdb,"Homo_sapiens.GRCh38.101.txdb");
k <- keys(tfdb, keytype = "TXNAME");reftx2gene <- select(tfdb, k, "GENEID", "TXNAME");#saveRDS(tx2gene,"Homo_sapiens.GRCh38.101.tx2gene.rds")

# rezdir="/Users/remy.nicolle/Downloads/PilotMillisectSmarter_CROS_090522results/kallistoOut"
# reftxgene="/Users/remy.nicolle/Downloads/PilotMillisectSmarter_CROS_090522results/reftx2gene.rds"

# reftx2gene=readRDS(reftxgene)


pattern="_R1.fastq.gz|_R1_001.fastq.gz"

library(parallel)
library(rjson)
# library(tximportData)
library(tximport)

alldir=list.dirs(rezdir,recursive=F)
# allfq=scan(initSlist,what="character",sep="\n")
allsampl=basename(alldir)

I=1:length(alldir)



runinfo=do.call(rbind,lapply(I,function(i){
  # print(i)
  as.data.frame(fromJSON(file=file.path(alldir[i],"run_info.json")),stringsAsFactors=F,row.names= allsampl[i])
}))

files <- setNames(file.path(alldir, "abundance.tsv"),allsampl)
txiout <- tximport(files, type = "kallisto", txOut=T)
# txgout <- tximport(files, type = "kallisto", tx2gene=reftx2gene,ignoreTxVersion=T)
txgoutraw=summarizeToGene(txiout, reftx2gene,  ignoreTxVersion = T, countsFromAbundance = "no")
# txgoutsctpm=summarizeToGene(txiout, reftx2gene,  ignoreTxVersion = T, countsFromAbundance = "scaledTPM")
# txgoutlsctpm=summarizeToGene(txiout, reftx2gene,  ignoreTxVersion = T, countsFromAbundance = "lengthScaledTPM")
  #c("no", "scaledTPM","lengthScaledTPM"))


write.table( txgoutraw$counts, file="kallistoGeneCount.tsv",quote=F,row.names=T,col.names=T, sep="\t")
system("gzip kallistoGeneCount.tsv")

write.table( txgoutraw$abundance, file="kallistoGeneAbund.tsv",quote=F,row.names=T,col.names=T, sep="\t")
system("gzip kallistoGeneAbund.tsv")


write.table( txiout$abundance , file="kallistoTranscriptAbund.tsv",quote=F,row.names=T,col.names=T, sep="\t")
system("gzip kallistoTranscriptAbund.tsv")

write.table( runinfo[colnames(txgoutraw$abundance),],file="kallistoInfoMetrics.tsv",quote=F,row.names=T,col.names=T, sep="\t")



# kallexp=list(
#   countraw=txgoutraw$counts,
#   countbias=txgbiasoutraw$counts,
#   TPMraw=txgoutraw$abundance,
#   TPMbias=txgbiasoutraw$abundance
# )




# outdir,paste0(prefix,"AllKallistoAllExp.rds")))


# saveRDS(kallexp,file=file.path(outdir,paste0(prefix,"AllKallistoAllExp.rds")))
# saveRDS(kallexp$countraw,file=file.path(outdir,paste0(prefix,"Kallisto_CountRaw.rds")))
# saveRDS(kallexp$TPMraw,file=file.path(outdir,paste0(prefix,"Kallisto_TPM.rds")))
# saveRDS(runinfo[colnames(kallexp$countraw),],file=file.path(outdir,paste0(prefix,"CountSummary.rds")))
