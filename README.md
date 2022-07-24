# RNA-seq pipeline #

This pipeline is developed using Nextflow and this documentation is intended to run the pipeline either locally or on Google Cloud.

Requires Nextflow and Docker installation is required \
the two docker images containning all the tools required by the pipeline are available on Docker Hub : \
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq \
https://hub.docker.com/repository/docker/genehetx/vep_hs \
This pipeline has two modules:
* rna_seq_pipe.nf : contains Fastqc, STAR alignment, FeatureCounts(reads quantification) and MultiQC
* variant_calling.nf : contains the GATK4 workflow for varaiant calling (for RNA-Seq data), and VEP (VARIANT EFFECT PREDICTOR /ENSEMBL) for variant annotation

The aims of this pipeline are :

* To check the quality of the fastq files (FatsQC)
* To align to a reference genome/transcriptome using (STAR)
* To quantify the expression quantification using (FeatureCounts)
* To quantify abundances of transcripts using (Kallisto)
* To identify variations (Variant calling) using (GATK4)
* To functionaly annotate variants using (VEP)
* To summarize of the FeatureCounts and STAR quality using (MultiQC)


## 0. Pre-requisites ##

In order to run correctly, three main variables :
* _Input_dir_ : the directory with all fastqs (not in sub directories)
* _output_dir_ : the directory in which all results of the pipe should be stored
* _ref_ : the directory with the reference genome information

:warning: __if running on google cloud__ : (Check the section 5 "Google Cloud execution) 

:warning: __if ref is not available__ : (check section 2 "Data preparation" step 2 : generate indexes ...) 



## 1. Clone this repository ##
- The first step is the git cloning of the pipeline directory 

```git clone https://github.com/GeNeHetX/Rna-seq_pipeline``` 

- The second step is considering the Rna-seq_pipeline directory as your work directory 

```cd Rna-seq_pipeline```

## 2. Data preparation ##
 1. Generate a sample list containing all of your samples names without the suffixes using the following bash comman line: \
 * __For single end data:__ \
  ```ls path/to/your/Inputdir |sed -e 's/\_R1.fastq.gz$//' > samlist.txt  ```
  - Inputdir : is the directory containing only your fastq files and nothing else \ 
  - _R1.fastq.gz : is an example of a suffix and it can be diffrenet from a datset to another \
  -samlist.txt: is the output of the command , it is a txt file that contains a list of you samples name \
  ** Example: 
  ```
  Inputdir : \
    * sample1_R1.fastq.gz
    * sample2_R1.fastq.gz 
    
   samlist.txt:
    * sample1
    * sample2
  ```
  
 * __For paired end data__ : \
  ```ls fastq_dir |sed -e 's/\_R1.fastq.gz$//'|sed -e 's/\_R2.fastq.gz$//'|uniq > samlist.txt ``` \
 

 2. Generate indexes required for each step of the pipeline
You will find in the ref_build.sh bash script, all the command lines that will help you to generate them. Please check the ref_build.sh file to understand the aim of each command line. \
For this step you will need: 
     * Reference genome file: (fasta)  GRCh38.p13 you can retrieve it from Ensembl Database.
     * Gene annotation file : (GTF) you can retrieve it from the Ensembl Database.
     * To make sure that the FASTA and the GTF belong to the same genome version !!
     * To get a Known-variants file: You can retrive it from Ensembl Database. 

## 3. setting up  ## 
The pipeline can be executed on a local computer or on Google Cloud Life Science platforme \
For a local execution modify the nextflow.config file and for Google Cloud execution modify the nextflowGCP.config, by changing the following parameters if necessary :
  * To run the piepline correctly, Please modify the following parameters (Mandatory):
    * ```params.kallisto ``` = true , by default, but if yo don't want to execute kallisto put false for this parameter 
    * ```params.variant_calling ``` = true,  by default , but if you don't want to execute Variant calling put "false"
    * ```params.outputdir=``` "/path/to your/outputdir" -> specify the path to your output directory (you should create a directory) where the pipeline outputs will be stored
    * ```params.sampleInputDir =``` "/path/to your/inputdir"  -> the directory that contains your raw fastqc files => for a local eexcution \
    this parameters will change when we executeon Google Cloud platefome to:  ```params.sampleInputDir =``` "gs://bucket_name/sampleInputDir" -> google cloud paths start with gs, followed by the bucket name \
    * ```params.sampleList =``` "/path/to your/samlist.txt"  -> the text file that contains a list of your fastqc sample names  generates in __Data preparation step__
    * ```params.samPsuffix1=```  -> specify the suffix of your fastq file name (ex: _R1_001)(if you have single end data you only need to specify this suffix without the sencond one (params.samPsuffix2) 
    * ```params.samPsuffix2=``` ->specify the suffix of your fastq file name (ex: _R2_001) -> (for paired end you need to specify both suffix1 and suffix2)
    * ```params.ref= ```"/PATH/to/ensembl_v105_GRCh38_p13" -> specify the path to the directory that  contains all the reference data for the pipeline execution) (for a local execution (generated using red_build.sh) in __Data preparation step__, and for Google Cloud execution modify ot this way  ```params.ref = ``` "gs://bucket_name/ensembl_v105_GRCh38_p13" \
 * optional parameters : The following parameters are for STAR aligner and Kallisto you can specify the values you want or keep the default ones (available on the config file)
 * STAR 
   * ```params.alignIntronMax = ```val
   * ```params.alignMatesGapMax=``` val  
   * ```params.limitOutSJcollapsed``` =val  
   * ```params.limitSjdbInsertNsj``` =val
   * ```params.outFilterMultimapNmax``` =val
   * ```params.winAnchorMultimapNmax``` =val  
   * ```params.alignSJoverhangMin``` =val
   * ```params.alignSJDBoverhangMin``` =val  
   * ```paramsalignIntronMin``` =val
   * ```params.outFilterMatchNminOverLread``` =val
   * ```params.outFilterScoreMinOverLread``` = val
   * ```params.outFilterMismatchNmax``` = val  
   * ```params.outFilterMismatchNoverLmax``` = val  
 * Kallisto 
   * ```params.bootstrap``` = 100
   the following two parameters only needed for single end data
   * ```params.read_len``` = 120 
   * ```params.read_sd``` =20 


## 4. Local Pipeline executioon 
For the execution setup your working directory to "workflow" which is the directory containning the workflows 
```cd workflow``` 
* a) For signle end data: \
```nextflow run single_end.nf -c ../nextflow.config  -w /path/to/your/workdir  -with-report``` 

:warning: The single_end.nf workflow accepts only single end data  

* b) For paired end data : \
```nextflow run paired_end_pipe.nf -c ../nextflow.config  -w /path/to/your/workdir  -with-report``` 

:warning:The paired_end_pipe.nf Wokflow accepts only paired end data 

for the -w : you have to specify the name of your work directory otherwise nextflow will name it "work" the work directory is different from the your output directory \
-c : specify the path to the config file\
-with-report : allows you to generate a report about the pipeline execution

## 5. Google Cloud pipeline execution ##
This step requires : 
  * Install Google Cloud SDK
  * Create a project on Google Cloud Life Science
  * Create a Bucket for your project
  * Generate a Json Key
For more information about the previous steps please Check the Google Cloud Documentation (https://cloud.google.com/life-sciences/docs/tutorials/nextflow) \
To execute the pipeline please follow these insctructions:
  1. Log in to Google Cloud
  2. ```Export GOOGLE_APPLICATION_CREDENTIALS=${PWD}/KEY_FILENAME.json (activate the json  key)``` (activation of the json key on youre work_directory)
  3. Copy all your fastq files and the directory generated by the ref_build.sh  using the following command : gsutil cp -r dir1/dir2 gs://my-bucket.
  
   d) Modify the machines capacities on the netflowGCP.config (Cpus, RAM, Disk, container) (if you want) 
  ``` withName: doSTAR{ \
        cpus = 16 \
        container = 'genehetx/genehetx-rnaseq:latest' \
        memory = 40.GB \
        disk = 1.TB \
    }
 ```
   e) Specify the name of  your Project, and the region where your data will be stored 
   ``` 
  google { \
    project = 'Project_name' \
    zone = 'europe-west4-a' \
    lifeSciences.bootDiskSize=80.GB
    google.lifeSciences.preemptible=true \
}
```
 5. Pipeline Execution :
 * a) Single end data : \
 ```nextflow run single_end.nf -c ../nextflowGCP.config  -w /path/to/your/workdir  -with-report [file name]```
 You can specify the name of your pipeline report [file name]
 * b) Paired end data: \
 ```nextflow run paired_end_pipe.nf -c ../nextflowGCP.config  -w /path/to/your/workdir  -with-report [file name]```
