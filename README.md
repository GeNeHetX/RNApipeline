# RNA-seq pipeline #

This pipeline is developed using Nextflow and this documentation is intended to run the pipeline either locally or on Google Cloud.

Requires Nextflow and Docker installation is required \
the two docker images containning all the tools required by the pipeline are available on Docker Hub : \
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq \
https://hub.docker.com/repository/docker/genehetx/vep_hs \
This pipeline has two modules:
* rna_seq_pipe.nf : contains  Fastqc , STAR alignment , FeatureCounts( reads quantification) and MultiQC
* variant_calling.nf : contains the GATK4 workflow for varaiant calling (for RNA -seq data), and VEP (VARIANT EFFECT PREDICTOR /ENSEMBL) for variant annotation

The aims of this pipeline are :

* check the quality of the fastq files
* align to a reference genome/transcriptome using STAR
* expression quantification with FeatureCounts
* Variant calling (Using GATK4)
* Variants annotation (VEP)
* summarize of the FeatureCounts and STAR quality using MultiQC


## 0. Pre-requisites ##

In order to run correctly, three main variables :
* _fastq_dir_ : the directory with all fastqs (not in sub directories)
* _output_dir_ : the directory in which all results of the pipe should be stored
* _ref_genome_ : the directory with the reference genome information


:warning: __if running on google cloud__ : do this before ... (or see section here called ....)

:warning: __if _ref_genome_  is not aailable __ : do this and that as described in section something



## 1 Clone this repository ##

git clone https://github.com/GeNeHetX/Rna-seq_pipeline \
cd Rna-seq_pipeline

### 2. Data preparation ###
 1. Generate a sample list containning all of your samples names using the following bash comman line: \
  ls fastq_dir |sed -e 's/\_R1.fastq.gz$//' > samlist.txt (for single end data) \
  ls fastq_dir |sed -e 's/\_R1.fastq.gz$//' |sed -e 's/\_R2.fastq.gz$//'|uniq > samlist.txt (for paired end data) \
  fastq_dir: is directory containning all your fastq files (make sure your only have fastq files you want to analyse)

 2. Geneate indexes required for each step of the pipeline
You will find in the ref_build.sh bash script, all the command lines that will help you to generate them. Please check the ref_build.sh file to understand the aim of each command line.


3. Modify the config file by changing the following parameters if necessary :
  * Please modify the following parameters :
    * params.outputdir="/path/to your/outputdir" -> specify the path to your output directory (you should create a directory)\
    * params.sampleInputDir = "/path/to your/inputdir"  -> the directory that contains your raw fastqc files\
    * params.sampleList = "/path/to your/samlist.txt"  -> the text file that contains a list of your fastqc sample names  generates in the first step (mentionned above)\
    * params.samPsuffix1=  -> specify the suffix of your fastq file name (ex: _R1_001)\
    * params.samPsuffix2=  ->specify the suffix of your fastq file name  (ex: _R2_001) -> (for paired end)\
    * params.ref="/PATH/to/ensembl_v105_GRCh38_p13" -> specify the path to the directory that  contains all the reference data for the pipeline execution (generated using red_build.sh)
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


### 1.1.2 Local Pipeline execution ###

cd workflow \
nextflow run single_end.nf -c ../Rna-seq_pipeline/nextflow.config  -w /path/to/your/workdir  -with-report \
nextflow run paired_end_pipe.nf -c ../Rna-seq_pipeline/nextflow.config  -w /path/to/your/workdir  -with-report \
The single_end.nf workflow accepts only single end data and exectues : the FastQC, STAR, FeatureCounts and MultQC processes \
The paired_end_pipe.nf Wokflow accepts only paired end data and executes : the FastQC, STAR, FeatureCounts, GATK4, Vep and MultQC processes \

for the -w : you have to specify the name of your work directory otherwise nextflow will name it "work"\
-c : specify the path to the config file\
-with-report : allows you to generate a report about the pipeline execution

### 1.1.3 Google Cloud pipeline execution ###
This step requires : \
  * Install Google Cloud SDK
  * Create a project on Google Cloud Life Science
  * Create a Bucket for your project
  * Generate a Json Key
For more information about the previous steps please Check the Google Cloud Documentation (https://cloud.google.com/life-sciences/docs/tutorials/nextflow) \
To execute the pipeline please follow these insctruction :
  1. log in to Google Cloud
  2. Export GOOGLE_APPLICATION_CREDENTIALS=${PWD}/KEY_FILENAME.json (activate the json  key)
  3. Copy all your fastq fileS and the directory generated by the ref_build.sh  using the following command : gsutil cp -r dir1/dir2 gs://my-bucket.
  4. Modify the paramaters on the nextflowGCP. config :
   a) params.sampleInputDir = "gs://bucket_name/fastq_dir" -> google cloud paths start with gs, followed by the bucket name
   b) params.ref = "gs://bucket_name/ensembl_v105_GRCh38_p13"
   c) for the rest of the pipeline parameters you can change them by following the isntruction mentionned above.
   d) Modify the machines capacities  (Cpus, RAM, Disk, container) (if you want) \
   withName: doSTAR{ \
        cpus = 16 \
        container = 'genehetx/genehetx-rnaseq:latest' \
        memory = 40.GB \
        disk = 1.TB \
    }

   e) Specify the name of  your Project, the region where your data will be stored \
   google { \
    project = 'Project_name' \
    zone = 'europe-west4-a' \
    lifeSciences.bootDiskSize=80.GB
    google.lifeSciences.preemptible=true \
}

 5. Pipeline Execution :
 nextflow run single_end.nf -c ../Rna-seq_pipeline/nextflowGCP.config  -w /path/to/your/workdir  -with-report\
 nextflow run paired_end_pipe.nf -c ../Rna-seq_pipeline/nextflowGCP.config  -w /path/to/your/workdir  -with-report
