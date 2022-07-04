# RNA-seq pipeline 
This pipeline is developed using Nextflow. \
Nextflow and Docker installation is required \
the two docker images containning all the tools required by the pipeline are available on Docker Hub : \
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq \
https://hub.docker.com/repository/docker/genehetx/vep_hs

## clone this repository 

git clone https://github.com/GeNeHetX/Rna-seq_pipeline \
cd Rna-seq_pipeline 
### Data preparation 
* modify the config file by changing the following parameters if necessary :
params.outputdir="/path/to your/outputdir" -> put the path tp your output directory(you should create a directory) \
params.sampleInputDir = "/path/to your/inputdir" #the directory that contains your raw fastqc files \
params.sampleList = "/path/to your/samlist.txt" #file that contains a list of your fastqc sample names \
params.samPsuffix1= #specify the suffix name of your fastq file (ex: _R1_001) \
params.samPsuffix2=#specify the suffix name of your fastq file (ex: _R2_001) -> for paired end \
params.ref="/PATH/to/ensembl_v105_GRCh38_p13"#specify the path to this directory that all the reference data for the pipeline execution very important!\
* The folllowing parameters are for STAR aligner you can specify the values you want or keep the default ones (available on the config file) 
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

to generate params.ref you can use the bash script "ref_build.sh" \
to generate the sample list : ls fastq_dir|sed -e 's/\_R1.fastq.gz$//' >samlist.txt (for single end) \ 
                              ls fastq_dir|sed -e 's/\_R1.fastq.gz$//'|sed -e 's/\_R2.fastq.gz$//'>samlist.txt (for paired end )

## Pipeline execution 
cd Workflow \
nextflow run single_end.nf -c ../Rna-seq_pipeline/nextflow.config  -w /path/to/your/workdir  -with-report \
nextflow run paired_end_pipe.nf -c ../Rna-seq_pipeline/nextflow.config  -w /path/to/your/workdir  -with-report \
for the -w : you have to specify the name of your work directory otherwise nextflow will name it "work" 


