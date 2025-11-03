#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5GB
#SBATCH -c 1
#SBATCH -t 23:59:59

module load nextflow/24.04.1


### To be customized by user:
PROJECTID="PROJECT_NAME"
INPUT_FQ_DIR="/path/to/fastqDir/"
REZ_DIR="/path/to/resultsDir/"$PROJECTID"/"
SAMPLE_CSV="/path/to/sample_csv/sample_list.csv"
CONFIG="/path/to/configFile/IFB_PE_minimal.config"

REF="/path/to/pdacrna/ensembl_v107_GRCh38b_kallisto_v0.51"
RNAPIPE_DIR="/path/to/RNApipeline/"

### Don't touch variables below this line
nextflow -c $CONFIG run $RNAPIPE_DIR/RNApipeline.nf -entry Main\
    -with-report report_$PROJECTID.html -resume \ 
    --csvSample $SAMPLE_CSV \
    --ref $REF \
    --sampleInputDir $INPUT_FQ_DIR \
    --scriptDir ${RNAPIPE_DIR}"/PrePostScripts/" \
    --logBackupDir "." \
    --runNumber $PROJECTID \
    --outputdir $REZ_DIR \