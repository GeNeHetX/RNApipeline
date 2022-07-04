# RNA-seq pipeline 
This pipeline is developed using Nextflow.
Nextflow and Docker installation is required 
the two docker images containning all the tools required by the pipeline are available on Docker Hub : 
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq 
https://hub.docker.com/repository/docker/genehetx/vep_hs

## Data preparation 
* downoald the single_end_pipeline directory 
* modify the config file by changing the following parameters if necessary :
params.outputdir="/path/to your/outputdir"
params.sampleInputDir = "/path/to your/inputdir" #the directory that contains your fastq files 
params.sampleList = "/path/to your/samlist.txt" #file thta contains your samples names 
params.samPsuffix = #specify the suffix name of your fastq file (ex: _R1_001)
params.ref="/PATH/to/ensembl_v105_GRCh38_p13"#specify the path to this directory thta contains all the reference data (indexes) very important !
* The folllowing parameters are for STAR aligner you can specify the values you want or keep the default ones :
params.alignIntronMax =val 
params.alignMatesGapMax= val  
params.limitOutSJcollapsed =val  
params.limitSjdbInsertNsj =val
params.outFilterMultimapNmax =val 
params.winAnchorMultimapNmax =val  
params.alignSJoverhangMin =val 
params.alignSJDBoverhangMin =val  
paramsalignIntronMin =val 
params.outFilterMatchNminOverLread =val 
params.outFilterScoreMinOverLread = val 
params.outFilterMismatchNmax = val  
params.outFilterMismatchNoverLmax = val  


## Running the pipeline on Google Cloud :

To run this pipeline Nextflow installation in required 

# Running the pipe 

cd single_end_pipeline

nextflow run single_end_pipe.nf -c nextflowGC.config -with-trace  -w gs://genehetxtst/<dirName>


