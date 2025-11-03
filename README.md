# RNA-seq pipeline #

This pipeline is developed using Nextflow and this documentation is intended to run the pipeline either locally or on Google Cloud.

Requires Nextflow and Docker installation is required \
the two docker images containning all the tools required by the pipeline are available on Docker Hub : \
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq \
https://hub.docker.com/repository/docker/genehetx/vep_hs \
This pipeline has two main modules:
* rna_seq_pipe.nf : contains Fastqc, STAR alignment, FeatureCounts(reads quantification) and MultiQC
* variant_calling.nf : contains the GATK4 workflow for varaiant calling (for RNA-Seq data), and VEP (VARIANT EFFECT PREDICTOR /ENSEMBL) for variant annotation

The aims of this pipeline are :

* To check the quality of the fastq files (FatsQC)
* To align to a reference genome/transcriptome using (STAR)
* To quantify the expression quantification using (FeatureCounts)
* To quantify abundances of transcripts using (Kallisto)
* To identify variations (Variant calling) using (GATK4 and/or Mpileup)
* To functionaly annotate variants using (VEP)
* To summarize of the FeatureCounts and STAR quality using (MultiQC)


## 0. Pre-requisites ##

In order to run correctly, four main variables :
* _Input_dir_ : the directory with all fastqs (not in sub directories)
* _output_dir_ : the directory in which all results of the pipe should be stored
* _ref_ : the directory with the reference genome information
* _config_file_ : the file to complete according your sequencing library (template available in **configExamples** directory)

:warning: __if running on google cloud__ : (Check the section 5 "Google Cloud execution) 

:warning: __if ref is not available__ : (check section 2 "Data preparation" step 2 : generate indexes ...) 



## 1. Clone this repository ##
- The first step is the git cloning of the pipeline directory 

```git clone https://github.com/GeNeHetX/RNApipeline``` 


## 2. Data preparation ##
 1. **Generate a sample list** containing all of your samples names without the suffixes and their associated project name using the following bash comman line: 
  * __For single end data:__ \
  ```ls path/to/your/Inputdir | sed -e 's/_R1\.fastq\.gz$//' | awk -v proj="projectName" '{print $0 "," proj}' > samlist.csv ```
  - Inputdir : is the directory containing only your fastq files and nothing else \ 
  - _R1.fastq.gz : is an example of a suffix and it can be diffrenet from a datset to another 
  - projectName : the sample's project name
  - samlist.csv: is the output of the command , it is a txt file that contains a list of you samples name  
  **Example**: 
  ```
  Inputdir : 
    * sample1_R1.fastq.gz
    * sample2_R1.fastq.gz 
    
   samlist.csv:
    * sample1, projectName
    * sample2, projectName
  ```   

   * __For paired end data__ :  
  ```ls path/to/your/Inputdir | sed -E 's/_R[12]\.fastq\.gz$//' | sort -u | awk -v proj="projectName" '{print $0 "," proj}' > samlist.csv ``` \
 

 2. **Generate indexes** required for each step of the pipeline
* You will find in the ref_build.sh bash script, all the command lines that will help you to generate them. Please check the ref_build.sh file to understand the aim of each command line. \
For this step you will need: 
     * Reference genome file: (fasta)  GRCh38.p13 you can retrieve it from Ensembl Database.
     * Gene annotation file : (GTF) you can retrieve it from the Ensembl Database.
     * Transcriptome file: (Cdna) you can retrieve it from the Ensembl Database. 
     * Make sure that the FASTA and the GTF belong to the same genome version !!
     * Get a Known-variants file: You can retrive it from Ensembl Database. 
This step requires 32Go RAM, so it is advised to generate once for the same genome, however you can generate it using the pipeline also by changing parameter in the config files (nextflow.config and nextflowGCP.config)=> explained in the __setting up__ setp.

## 3. setting up  ## 
The pipeline can be executed on a local computer, on a Slurm cluster (like the IFB core) or on Google Cloud Life Science platforme \
For a local execution modify the local.config file and for Google Cloud execution modify the GoogleCloud.config or GCP_PE_minimal.config, by changing the following parameters if necessary :
* To run the piepline correctly, Please modify the following parameters (Mandatory):
  * Library parameters:
    * ```params.single_end``` = true if single-end data, = false if paired-end data
    * ```params.FeatureCountStrand``` = 2 for Smarter Library Kit (you can get this information using salmon)
    * ```params.kallistoStrand``` = "--rf-stranded" for Smarter Library Kit (you can get this information using salmon)
  * Tools selection : 
    * ```params.star``` = true , by default, but if yo don't want to execute STAR put false for this parameter 
    * ```params.fastq``` = true , by default, but if yo don't want to execute Fastqc put false for this parameter 
    * ```params.multiqc``` = true , by default, but if yo don't want to execute Multiqc put false for this parameter 
    * ```params.kallisto``` = true , by default, but if yo don't want to execute Mallisto put false for this parameter 
    * ```params.fcounts``` = true , by default, but if yo don't want to execute FeatueCpounts put false for this parameter 
    * ```params.gatk4``` = true , by default, but if yo don't want to execute GATK4 put false for this parameter 
    * ```params.samtools_depth``` = true , by default, but if yo don't want to execute Samtools depth put false for this parameter 
    * ```params.no_multimapped``` = true , by default, filter multimapped read in Samtools depth command
    * ```params.mpileup``` = true , by default, but if yo don't want to execute Mpileup put false for this parameter 
    * ```params.vep``` = true , by default, but if yo don't want to execute VEP put false for this parameter 
  * Reference paths :
    * ```params.ref= ```"/PATH/to/ensembl_v105_GRCh38_p13" -> specify the path to the directory that  contains all the reference data for the pipeline execution) (for a local execution (generated using red_build.sh) in __Data preparation step__
    * But if you want to include the indexes  generation in the pipeline you have to specify the parameter like this ```params.ref = ``` "no_ref"
    * ```params.vep_cache``` = "/PATH/to/ensembl_v105_GRCh38_p13/VEP"
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
For the execution setup, go to your working directory :

* a) In local or in Google Cloud : \
```bash
  cd workflow_dir

  ### To be customized by user:
  PROJECTID="PROJECT_NAME"
  INPUT_FQ_DIR="/path/to/fastqDir/"
  REZ_DIR="/path/to/resultsDir/"$PROJECTID"/"
  SAMPLE_CSV="/path/to/sample_csv/sample_list.csv"
  CONFIG="/path/to/configFile/IFB_PE_minimal.config"

  BED="/path/to/bedFile/my_hotspots.bed"
  REF="/path/to/pdacrna/ensembl_v107_GRCh38b_kallisto_v0.51"
  RNAPIPE_DIR="/path/to/RNApipeline/"

  ### Dont touch variables below this line
  nextflow -c $CONFIG run $RNAPIPE_DIR/RNApipeline.nf -entry Main\
      -with-report report_$PROJECTID.html -resume \ 
      --csvSample $SAMPLE_CSV \
      --ref $REF \
      --sampleInputDir $INPUT_FQ_DIR \
      --scriptDir ${RNAPIPE_DIR}"/PrePostScripts/" \
      --logBackupDir "." \
      --runNumber $PROJECTID \
      --outputdir $REZ_DIR \
      --bed $BED
``` 

* b) In a slurm cluster (like IFB) : \
copy the runIFBjob.sh exemple in the repository and modify the user parameters (see above)
```bash
  cd workflow_dir

  # modify runIFBjob.sh
  sbatch runIFBjob.sh
``` 

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