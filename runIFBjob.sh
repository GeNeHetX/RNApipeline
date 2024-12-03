#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5GB
#SBATCH -c 1


module load nextflow/24.04.1

ROOTDIR="/shared/projects/your_project"
PROJECTID="exemple"
INPUT_FQ_DIR=$ROOTDIR"/DUMP_TIPMP/fastq_merged_lanes"
REZ_DIR=$ROOTDIR"/Results/"$PROJECTID"/"
REF="/shared/projects/pdacrna/ensembl_v113_GRCh38"

## Adjust this command by keeping or removing the -resume option
nextflow -c $2 run $1 -resume -with-report --sampleInputDir $INPUT_FQ_DIR --ref $REF --outputdir $REZ_DIR --sampleList $3
