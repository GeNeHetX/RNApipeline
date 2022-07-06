# RNA-seq pipeline #
This pipeline is developed using Nextflow. \
Nextflow and Docker installation is required \
the two docker images containning all the tools required by the pipeline are available on Docker Hub : \
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq \
https://hub.docker.com/repository/docker/genehetx/vep_hs \
This pipeline has two modules: \ 
* rna_seq_pipe.nf : contains  Fastqc , STAR alignment , FeatureCounts( reads quantification) and MultiQC
* variant_calling.nf : contains the GATK4 workflow for varaiant calling (for RNA -seq data), and VEP (VARIANT EFFECT PREDICTOR /ENSEMBL) for variant annotation \

The aim of this pipeline is : \
check the quality of the fastq files, align to a reference genome (STAR) ,expression quantification with FeatureCounts, Varinat calling (Using GATK4), Variants annotation (VEP) and a summary of the FeatureCounts and STAR quality using MultiQC.


## Clone this repository ##

git clone https://github.com/GeNeHetX/Rna-seq_pipeline \
cd Rna-seq_pipeline 

### Data preparation ###
 1. Generate a sample list containning all of your samples names using the following bash comman line : 
  ls fastq_dir |sed -e 's/\_R1.fastq.gz$//' > samlist.txt (for single end data) \
  ls fastq_dir |sed -e 's/\_R1.fastq.gz$//' |sed -e 's/\_R2.fastq.gz$//'|uniq > samlist.txt (for paired end data) \
  fastq_dir: is directory containning all your fastq files (make sure your only have fastq files you want to analyse) 
  
2. Geneate indexes required for each step of the pipeline 
You will find in the ref_build.sh bash script, all the command lines that will help you to generate them. Please check the ref_build.sh file to understand the aim of each command line. \


3. Modify the config file by changing the following parameters if necessary :
  * Please modify the following parameters :
params.outputdir="/path/to your/outputdir" -> specify the path to your output directory (you should create a directory)\
params.sampleInputDir = "/path/to your/inputdir"  -> the directory that contains your raw fastqc files\
params.sampleList = "/path/to your/samlist.txt"  -> the text file that contains a list of your fastqc sample names  generates in the first step (mentionned above)\
params.samPsuffix1=  -> specify the suffix of your fastq file name (ex: _R1_001)\
params.samPsuffix2=  ->specify the suffix of your fastq file name  (ex: _R2_001) -> (for paired end)\
params.ref="/PATH/to/ensembl_v105_GRCh38_p13" -> specify the path to the directory that  contains all the reference data for the pipeline execution (generated using red_build.sh)
 * optional parameters : The folllowing parameters are for STAR aligner you can specify the values you want or keep the default ones (available on the config file)
  * params.alignIntronMax =val 
  * params.alignMatesGapMax= val  
  * params.limitOutSJcollapsed =val  
  * params.limitSjdbInsertNsj =val
  * params.outFilterMultimapNmax =val 
  * params.winAnchorMultimapNmax =val  
  * params.alignSJoverhangMin =val 
  * params.alignSJDBoverhangMin =val  
  * paramsalignIntronMin =val 
  * params.outFilterMatchNminOverLread =val 
  * params.outFilterScoreMinOverLread = val 
  * params.outFilterMismatchNmax = val  
  * params.outFilterMismatchNoverLmax = val  


## Local Pipeline execution ##

cd workflow \
nextflow run single_end.nf -c ../Rna-seq_pipeline/nextflow.config  -w /path/to/your/workdir  -with-report \
nextflow run paired_end_pipe.nf -c ../Rna-seq_pipeline/nextflow.config  -w /path/to/your/workdir  -with-report 

for the -w : you have to specify the name of your work directory otherwise nextflow will name it "work"\
-c : specify the path to the config file\
-with-report : allows you to generate a report about the pipeline execution 

## Google Cloud pipeline execution ##
 * install Google Cloud SDK 
 * Create a project on Google Cloud Life Science 
 * Create a Bucket for your project 
 * generate a Json Key 
For more information about the previous steps please Check the Google Cloud Documentation (https://cloud.google.com/life-sciences/docs/tutorials/nextflow) \
 1. Modify the paramaters on the nextflowGCP. config file as mentionned above.
   * modify the machines capacities  (Cpus, RAM, Disk) (if you want)
   * specify the name of  your Project
   * specify the zone  where your data will be stored (this parameters is specified during the bucket creation) 
 2. log in to Google Cloud
 3. export GOOGLE_APPLICATION_CREDENTIALS=${PWD}/KEY_FILENAME.json (activate the json  key)
 4. pipeline Execution :
 nextflow run single_end.nf -c ../Rna-seq_pipeline/nextflowGCP.config  -w /path/to/your/workdir  -with-report\
 nextflow run paired_end_pipe.nf -c ../Rna-seq_pipeline/nextflowGCP.config  -w /path/to/your/workdir  -with-report






