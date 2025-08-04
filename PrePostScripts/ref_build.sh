#!/bin/bash
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --partition fast
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64GB
#SBATCH -t 23:59:59

### BUILD REF FOR RNAPIPE 1.6.0 ###
######  Parameters to change ######
Ensemblv="113"
newDir=../ensembl_v${Ensemblv}_GRCh38
###################################

######  TOOLS TO PULL ######
#### Do it once
if [ ! -f kallisto_0.51.1--heb0cbe2_0.sif ]; then
    singularity pull docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0
fi
if [ ! -f genehetx-rnaseq_v1.6.0.sif ]; then
    singularity pull docker://genehetx/genehetx-rnaseq:v1.6.0
fi
############################

#create new ref directory and change access
mkdir -p $newDir &&\
chmod +rwx $newDir


##downoalding the gtf,fasta and cdna files
if [ ! -f Homo_sapiens.GRCh38.${Ensemblv}.chr.gtf.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${Ensemblv}.chr.gtf.gz
fi
if [ ! -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi
if [ ! -f Homo_sapiens.GRCh38.cdna.all.fa.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
fi

##downoalding the vcf file
if [ ! -f 1000GENOMES-phase_3.vcf.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz
fi

## downloading the vep file
if [ ! -f homo_sapiens_vep_${Ensemblv}_GRCh38.tar.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/variation/indexed_vep_cache/homo_sapiens_vep_${Ensemblv}_GRCh38.tar.gz
fi

## unzip ref files
gunzip -c Homo_sapiens.GRCh38.${Ensemblv}.chr.gtf.gz >ref.gtf
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ref.fa
gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz > cdna.fa

## unzip VEP caches directory
tar -xf homo_sapiens_vep_${Ensemblv}_GRCh38.tar.gz

##unzip the vcf file
gunzip -c 1000GENOMES-phase_3.vcf.gz > knowns_variants.vcf

##kallisto index
## run docker kallisto 
singularity exec kallisto_0.51.1--heb0cbe2_0.sif kallisto index -i kallisto_index cdna.fa

##create an  samtools index (needed for GATK4)
## run docker samtools
singularity exec genehetx-rnaseq_v1.6.0.sif samtools faidx ref.fa


##create a dict (needed for GATK4)
##you should specify the path to your picard.jar file or run docker
singularity exec genehetx-rnaseq_v1.6.0.sif java -jar /usr/local/bin/picard.jar  CreateSequenceDictionary R= ref.fa O= ref.dict

##indexinf the vcf file (needed for  haplotypeCaller)
##you should specify the path to your gatk jar file of the specific region you use
singularity exec genehetx-rnaseq_v1.6.0.sif java -jar /usr/local/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar IndexFeatureFile -I knowns_variants.vcf

##parsing a gtf file in order to get exon informations (needed for featureCount output analysis)
awk -F "\t" '$3 == "exon" { print $4"\t"$5"\t"$7"\t"$9 }' ref.gtf |awk '{for(i=5;i<=NF;i++){if($i~/^"ENSE/){a=$i}} print a, $1,$2,$3,$5,$15}'| sed 's/\"//g'|sed 's/\;//g'| sort -d | awk 'BEGIN {print "exon_id\tstart\tend\tstrand\tgene_id\tgene_name"} { print }' > Exon_gtf_info.tab

###
#Rajouter ref.ExonOnly.gtf pour BedInterval
###

grep -P "\tgene\t" ref.gtf > ref.GeneLvlOnly.gtf
path_Rscript=$(which Rscript)
${path_Rscript} /shared/projects/pdac_vcpipe/RNApipeline/PrePostScripts/procGTF.R ref.GeneLvlOnly.gtf refGeneID_ensembl_v${Ensemblv}

## Move files in the new directory
mv ref.gtf $newDir/ref.gtf
mv ref.fa $newDir/ref.fa
mv ref.fa.fai $newDir/ref.fa.fai
mv ref.dict $newDir/ref.dict
mv kallisto_index $newDir/kallisto_index
mv knowns_variants.vcf  $newDir/knowns_variants.vcf
mv knowns_variants.vcf.idx $newDir/knowns_variants.vcf.idx
mv Exon_gtf_info.tab $newDir/Exon_gtf_info.tab
mv ref.GeneLvlOnly.gtf $newDir/ref.GeneLvlOnly.gtf
mv refGeneID_ensembl_${Ensemblv}.* $newDir/.

#Transcriptome
mkdir -p $newDir/transcript/
cp cdna.fa $newDir/cdna.fa
mv cdna.fa $newDir/transcript/cdna.fa

#VEP files
mv homo_sapiens_vep_${Ensemblv}_GRCh38.tar.gz $newDir/.
mkdir -p $newDir/VEP/
mv homo_sapiens $newDir/VEP/.

##generate a STAR index
## run docker star
singularity exec genehetx-rnaseq_v1.6.0.sif STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $newDir --genomeFastaFiles $newDir/ref.fa  --sjdbOverhang 100 --sjdbGTFfile $newDir/ref.gtf  --genomeSAindexNbases 11

#copy slurm job
cp ${0} $newDir/build_ref${Ensemblv}.slurm