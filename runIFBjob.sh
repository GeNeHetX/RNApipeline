#!/bin/bash
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5GB
#SBATCH -c 1
#SBATCH -t 23:59:59



export JAVA_HOME=/shared/software/conda/envs/nextflow-23.10.1

module load nextflow/24.04.1

## Adjust this command by keeping or removing the -resume option
nextflow -c $2 run $1 -resume -with-report --sampleInputDir $3 --ref $4 --outputdir $5 --sampleList $6
