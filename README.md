# RNA-seq pipeline #
This pipeline allows performing RNA-seq analysis on single- and paired-end fastq files.
It is developed using Nextflow and can be run either locally, on Google Cloud, or on a slurm server.

This pipeline has two main modules:
* rna_seq_pipe.nf : contains Fastqc, STAR alignment, FeatureCounts(reads quantification) and MultiQC
* variant_calling.nf : contains the GATK4 workflow for variant calling (for RNA-Seq data), and VEP (VARIANT EFFECT PREDICTOR /ENSEMBL) for variant annotation

It is composed of the following steps :

* Check the quality of the fastq files (FatsQC)
* Align to a reference genome/transcriptome using (STAR)
* Quantify the RNA expression using (FeatureCounts)
* Quantify abundances of transcripts using (Kallisto)
* Identify nucleotide variations (Variant calling) using (GATK4 and/or Mpileup)
* Functionally annotate variants using (VEP)
* Summary of the FeatureCounts and STAR quality using (MultiQC)


## 0. Pre-requisites ##
Nextflow and Docker are required \
The two Docker images containing all the tools required by the pipeline are available on Docker Hub : \
https://hub.docker.com/repository/docker/genehetx/genehetx-rnaseq \
https://hub.docker.com/repository/docker/genehetx/vep_hs \

In order to run correctly, four main variables :
* _Input_dir_ : the directory with all fastqs (not in subdirectories)
* _output_dir_ : the directory in which all results of the pipe should be stored
* _ref_ : the directory with the reference genome information
* _config_file_ : the file to complete according to your sequencing library (template available in **configExamples** directory)

:warning: __if running on google cloud__ : (Check the section 5 "Google Cloud execution) 

:warning: __if ref is not available__ : (check section 2 "Data preparation" step 2 : generate indexes ...) 



## 1. Clone this repository ##
- The first step is the git cloning of the pipeline directory 

```git clone https://github.com/GeNeHetX/RNApipeline``` 


## 2. Data preparation ##
 1. **Generate a sample list** containing all of your sample names without the suffixes and their associated project name using the following bash comman line: 
  * __For single end data:__ 
  ```bash
  INPUT_DIR="path/to/your/Inputdir"
  PROJECT_NAME="projectName"
  (echo "ID_Sample,Project"; ls "${INPUT_DIR}" | sed -E 's/_R1\.fastq\.gz$//' | sort -u | awk -v proj="${PROJECT_NAME}" '{print $0 "," proj}') > samlist.csv
  ```
  - Inputdir : is the directory containing only your fastq files and nothing else \ 
  - _R1.fastq.gz : is an example of a suffix, and it can be different from one dataset to another 
  - projectName : the sample's project name
  - samlist.csv: is the output of the command, it is a txt file that contains a list of your sample names  
  **Example**: 
  ```
  Inputdir : 
    * sample1_R1.fastq.gz
    * sample2_R1.fastq.gz 
    
   samlist.csv:
    ID_Sample,Project
    sample1, projectName
    sample2, projectName
  ```   

   * __For paired end data__ :  
  ```bash
INPUT_DIR="path/to/your/Inputdir"
PROJECT_NAME="projectName"

(
echo "ID_Sample,Project"
ls "${INPUT_DIR}" \
| grep -E '_R[12](_[0-9]+)?\.fastq(\.gz)?$' \
| sed -E 's/_R[12](_[0-9]+)?\.fastq(\.gz)?$//' \
| sort -u \
| awk -v proj="${PROJECT_NAME}" '{print $0 "," proj}'
) > ./samlist.csv
  ```

 2. **Generate indexes** required for each step of the pipeline
* Use `PrePostScripts/ref_build.nf` to generate them. The workflow documents each independent reference-building step and validates the complete result before publication. \
For this step, you will need: 
     * Reference genome file: (fasta)  GRCh38.p13 you can retrieve it from Ensembl Database.
     * Gene annotation file : (GTF) you can retrieve it from the Ensembl Database.
     * Transcriptome file: (Cdna) you can retrieve it from the Ensembl Database. 
     * Make sure that the FASTA and the GTF belong to the same genome version !!
     * Get a Known-variants file: You can retrieve it from Ensembl Database. 
This step requires 32Go RAM, so it is advised to generate it once for a given
genome. The standard project parameters are in `nextflow.config`; other
execution environments can use the matching file under `configExamples/`.

### TOD/PAM infrastructure

On the TOD/PAM Slurm infrastructure, `/ref` is a direct CephFS mount. IAC
`/etc/nextflow/site.config` owns Slurm, scratch, Apptainer, caches, and the
explicit queue profiles. The repository `nextflow.config` supplies the standard
reference root, containers, and pipeline parameters. `nf-run NAME` supplies the per-run
`/biojobs/nextflow/NAME/work` directory. Repeating the same name resumes the
existing run; do not create a new random run name for every retry.

Build the Kallisto SIFs once, before the reference. The ARM64 SIF must be built
on PAM; the x86_64 SIF must be built on TOD (Intel and AMD are both x86_64).
The scripts are in `containers/kallisto/` and each takes one argument: its
destination directory. From a controller, submit them to the matching worker
partition; do not build the ARM64 image on `ctra`.

After both SIFs exist under `/ref/tools/rnapipeline/v1.7.0/kallisto/0.51.1/`,
install the architecture-selecting Kallisto launcher once:

```bash
install -D -m 0755 containers/kallisto/kallisto-wrapper.sh \
  /ref/tools/rnapipeline/v1.7.0/kallisto/0.51.1/bin/kallisto
```

Then run the reference workflow. `nf-run` runs the coordinator on `ctra`; the
workflow downloads and validates one homogeneous Ensembl release, then sends
the independent sequence, variant, Kallisto, GTF, and STAR tasks to Slurm in
parallel. Each task publishes its own completed files to the staging directory.
The final task only validates the files and renames the completed staging
directory into place:

```bash
source "$HOME/IAC/infra/scripts/shell-setup"
nf-run rnapipeline-ref-v107-test \
  /biojobs/pipelines/RNApipeline/PrePostScripts/ref_build.nf \
  -profile pam_cpu \
  --reference_id ensembl_v107_GRCh38_test
```

The reference workflow never creates SIFs and never invokes Apptainer directly
from a shell helper. The Kallisto launcher selects the SIF after Slurm assigns
the worker; the site Nextflow config supplies the `/ref` binding.
The same reference can therefore be built on PAM or TOD and used on both
architectures.

Use `-profile pam_cpu` for PAM CPU work, `-profile pam_gpu` for PAM GPU work,
or `-profile burst` after enabling infrastructure mode `burst`. Resume the same run name
after a failed task; do not create a new name for each retry.

For a normal pipeline run, choose a meaningful stable name and add the project
config after the site defaults:

```bash
nf-run RUN-NAME /biojobs/pipelines/RNApipeline/fullPairedEnd.nf \
  -profile pam_cpu
```

This publishes to `s3://process/rnapipeline/RUN-NAME`. Add
`--outputdir results` when durable output should stay in the run's local
`results/` directory instead.

### TOD/PAM S3 test run

Use one stable `nf-run` name per analysis. The run directory is
`/biojobs/nextflow/RUN-NAME`; its `metadata/samples.csv` is the input sheet,
`work/` is resumable Nextflow state, and `logs/` contains the engine log and
reports. The FASTQs may remain in S3; do not copy them to `/biojobs`.

For paired-end files named `SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`, use:

```csv
ID_Sample,Project,suffix1,suffix2
BCPEVP0625_R001,rna-prelimtest,_R1,_R2
BCPRNA0625_R038,rna-prelimtest,_R1,_R2
BCPRNA0625_R056,rna-prelimtest,_R1,_R2
```

Then launch from a controller:

```bash
RUN_NAME=rna-prelimtest
RUN_DIR=/biojobs/nextflow/$RUN_NAME
mkdir -p "$RUN_DIR/metadata"
# Create $RUN_DIR/metadata/samples.csv with the sheet above.

nf-run "$RUN_NAME" /biojobs/pipelines/RNApipeline/fullPairedEnd.nf \
  -profile pam_cpu \
  --sampleInputDir s3://sandbox/Fastq_tests \
  --csvSample "$RUN_DIR/metadata/samples.csv"
```

The profile is explicit. Add `-profile burst` only when infrastructure mode
`burst` is active. The pipeline publishes durable results to
`s3://process/rnapipeline/rna-prelimtest`; the same run's logs, metadata, and
resumable work remain under `/biojobs/nextflow/rna-prelimtest`.

Use `-profile burst` only while the infrastructure is in `burst` mode.
If optional AMD64-only tools such as VEP, DeepVariant, or the BioContainers
bcftools step are enabled, use burst mode. Their `amd64` process labels route
those tasks to TOD while compatible tasks may still use PAM. Architecture is
selected by process labels and containers, not by a user-facing profile.

The build uses the configured Kallisto runtime for the worker architecture and
requires all Ensembl inputs to use the same release. `--force` is only for an
intentional replacement; a normal retry uses the same run name without it.
The staging directory is owned by the stable `nf-run` name, so a different
run cannot accidentally overwrite an active build.

## 3. Setting up  ##
The pipeline can be executed on a local computer, on a Slurm cluster (like the IFB core) or on Google Cloud Life Science platform \
For a local execution, modify the local.config file and for Google Cloud execution, modify the GoogleCloud.config or GCP_PE_minimal.config, by changing the following parameters if necessary :
* To run the pipeline correctly, please modify the following parameters (Mandatory):
  * Library parameters:
    * ```params.single_end``` = true if single-end data, = false if paired-end data
    * ```params.FeatureCountStrand``` = 2 for Smarter Library Kit (you can get this information using salmon)
    * ```params.kallistoStrand``` = "--rf-stranded" for Smarter Library Kit (you can get this information using salmon)
  * Tools selection : 
    * ```params.star``` = true , by default, but if you don't want to execute STAR put false for this parameter 
    * ```params.fastq``` = true , by default, but if you don't want to execute Fastqc put false for this parameter 
    * ```params.multiqc``` = true , by default, but if you don't want to execute Multiqc put false for this parameter 
    * ```params.kallisto``` = true , by default, but if you don't want to execute Mallisto put false for this parameter 
    * ```params.fcounts``` = true , by default, but if you don't want to execute FeatueCpounts put false for this parameter 
    * ```params.gatk4``` = true , by default, but if you don't want to execute GATK4 put false for this parameter 
    * ```params.samtools_depth``` = true , by default, but if you don't want to execute Samtools depth put false for this parameter 
    * ```params.no_multimapped``` = true , by default, filter multimapped read in Samtools depth command
    * ```params.mpileup``` = true , by default, but if you don't want to execute Mpileup put false for this parameter 
    * ```params.vep``` = true , by default, but if you don't want to execute VEP put false for this parameter 
  * Reference paths :
    * ```params.ref= ```"/PATH/to/ensembl_v107_GRCh38" -> specify the path to the directory that contains all the reference data for pipeline execution (generated by the reference Nextflow workflow)
    * But if you want to include the indexes  generation in the pipeline you have to specify the parameter like this ```params.ref = ``` "no_ref"
    * ```params.vep_cache``` = "/PATH/to/ensembl_v107_GRCh38_p13/VEP"
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
   The following two parameters are only needed for single-end data
   * ```params.read_len``` = 120 
   * ```params.read_sd``` =20 


## 4. Local Pipeline execution 
For the execution setup, go to your working directory :

* a) In local or in Google Cloud: \
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

  ### Don't touch variables below this line
  nextflow -c $CONFIG run $RNAPIPE_DIR/fullPairedEnd.nf \
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

* b) In a slurm cluster (like IFB): \
Copy the runIFBjob.sh example in the repository and modify the user parameters (see above)
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
  * Generate a JSON Key
For more information about the previous steps, please check the Google Cloud Documentation (https://cloud.google.com/life-sciences/docs/tutorials/nextflow) \
To execute the pipeline, please follow these instructions:
  1. Log in to Google Cloud
  2. ```Export GOOGLE_APPLICATION_CREDENTIALS=${PWD}/KEY_FILENAME.json (activate the json  key)``` (activation of the json key on youre work_directory)
  3. Copy all your fastq files and the directory generated by the reference workflow using the following command: gsutil cp -r dir1/dir2 gs://my-bucket.
  
   d) Modify the machines' capacities in the selected Google Cloud config (CPUs, RAM, disk, container) (if you want)
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
 * a) Single-end data : \
 ```nextflow run simpleSingleEnd.nf -c ../configExamples/GoogleCloud.config  -w /path/to/your/workdir  -with-report [file name]```
 You can specify the name of your pipeline report [file name]
 * b) Paired end data: \
 ```nextflow run fullPairedEnd.nf -c ../configExamples/GoogleCloud.config  -w /path/to/your/workdir  -with-report [file name]```
