#!/bin/bash

#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --partition fast
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64GB
#SBATCH -t 23:59:59

### NEED R & apptainer installed on your system
# load module if needed
#module load r

###### Arguments ######

usage() {
    echo "Usage:"
    echo "  bash $0 -v ENSEMBL_VERSION -s SPECIES -p RNAPIPE_PATH"
    echo
    echo "Required arguments:"
    echo "  -v    Ensembl release version, for example 113"
    echo "  -s    Species: homo_sapiens or mus_musculus"
    echo "  -p    Path to the RNApipeline directory"
    echo
    echo "Optional arguments:"
    echo "  -h    Display this help message"
    echo
    echo "Example:"
    echo "  bash $0 -v 113 -s mus_musculus \\"
    echo "      -p /path/to/Github/RNApipeline"
}

Ensemblv=""
species=""
RNApipe=""

# getopts : argument entries and errors
while getopts ":v:s:p:h" option; do
    case "$option" in
        v)
            Ensemblv="$OPTARG"
            ;;
        s)
            species="$OPTARG"
            ;;
        p)
            RNApipe="$OPTARG"
            ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "Error: option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
        \?)
            echo "Error: unknown option -$OPTARG." >&2
            usage
            exit 1
            ;;
    esac
done

if [ -z "$Ensemblv" ] ||
   [ -z "$species" ] ||
   [ -z "$RNApipe" ]; then
    echo "Error: options -v, -s and -p are required." >&2
    echo >&2
    usage
    exit 1
fi

if ! [[ "$Ensemblv" =~ ^[0-9]+$ ]]; then
    echo "Error: Ensembl version must be an integer." >&2
    exit 1
fi

## to modify if you want to handle other species
if [ "$species" != "homo_sapiens" ] &&
   [ "$species" != "mus_musculus" ]; then
    echo "Error: unsupported species '$species'." >&2
    echo "Accepted values: homo_sapiens or mus_musculus." >&2
    exit 1
fi

if [ ! -d "$RNApipe" ]; then
    echo "Error: RNApipeline directory does not exist:" >&2
    echo "  $RNApipe" >&2
    exit 1
fi

#### To change if new ref
if [ "$species" == "homo_sapiens" ]; then
    prefix="Homo_sapiens"
    genome_version="GRCh38"
elif [ "$species" == "mus_musculus" ]; then
    prefix="Mus_musculus"
    genome_version="GRCm39"
fi

newDir="./ensembl_v${Ensemblv}_${species}_${genome_version}"

echo "Ensembl version : $Ensemblv"
echo "Species         : $species"
echo "Genome version  : $genome_version"
echo "RNApipeline     : $RNApipe"
echo "Output directory: $newDir"
echo "\n"

############################
############################

######  TOOLS TO PULL ######
#### Do it once
if [ ! -f kallisto_0.51.1--heb0cbe2_0.sif ]; then
    apptainer pull docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0
fi
if [ ! -f genehetx-rnaseq_v1.6.0.sif ]; then
    apptainer pull docker://genehetx/genehetx-rnaseq:v1.6.0
fi

############################

#create new ref directory and change access
mkdir -p $newDir &&\
chmod +rwx $newDir

##downoalding the gtf,fasta and cdna files
echo "Downloading reference files..."
if [ ! -f ${prefix}.${genome_version}.${Ensemblv}.chr.gtf.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/gtf/${species}/${prefix}.${genome_version}.${Ensemblv}.chr.gtf.gz
fi
if [ ! -f ${prefix}.${genome_version}.dna.primary_assembly.fa.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/fasta/${species}/dna/${prefix}.${genome_version}.dna.primary_assembly.fa.gz
fi
if [ ! -f ${prefix}.${genome_version}.cdna.all.fa.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/fasta/${species}/cdna/${prefix}.${genome_version}.cdna.all.fa.gz
fi

##downoalding the vcf file
if [ "$species" == "homo_sapiens" ] && [ ! -f 1000GENOMES-phase_3.vcf.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/variation/vcf/${species}/1000GENOMES-phase_3.vcf.gz
fi
if [ "$species" == "mus_musculus" ] && [ ! -f ${species}.vcf.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/variation/vcf/${species}/${species}.vcf.gz
fi

## downloading the vep file
if [ ! -f ${species}_vep_${Ensemblv}_${genome_version}.tar.gz ]; then
    wget http://ftp.ensembl.org/pub/release-${Ensemblv}/variation/indexed_vep_cache/${species}_vep_${Ensemblv}_${genome_version}.tar.gz
fi

echo "Done downloading reference files."

## unzip ref files
echo "Unzipping reference files..."
if [ ! -f ref.gtf ]; then
    echo "Unzipping GTF file..."
    gunzip -c ${prefix}.${genome_version}.${Ensemblv}.chr.gtf.gz > ref.gtf
fi
if [ ! -f ref.fa ]; then
    echo "Unzipping FASTA file..."
    gunzip -c ${prefix}.${genome_version}.dna.primary_assembly.fa.gz > ref.fa
fi
if [ ! -f transcriptom.fa ]; then
    echo "Unzipping cDNA file..."
    gunzip -c ${prefix}.${genome_version}.cdna.all.fa.gz > transcriptom.fa
fi

## unzip VEP caches directory
if [ ! -d ${species} ]; then
    echo "Unzipping VEP cache..."
    tar -xf ${species}_vep_${Ensemblv}_${genome_version}.tar.gz
fi

##unzip the vcf file
if [ -f ${species}.vcf.gz ] && [ ! -f knowns_variants.vcf ]; then
    gunzip -c ${species}.vcf.gz > knowns_variants.vcf
elif [ -f 1000GENOMES-phase_3.vcf.gz ] && [ ! -f knowns_variants.vcf ]; then
    gunzip -c 1000GENOMES-phase_3.vcf.gz > knowns_variants.vcf
fi
echo "Done Unzipping reference files."

echo "Creating tool indexes ..."
##kallisto index
## run docker kallisto 
if [ ! -f kallisto_index ]; then
    echo "Creating Kallisto index..."
    apptainer exec kallisto_0.51.1--heb0cbe2_0.sif kallisto index -i kallisto_index transcriptom.fa
fi 
##create an  samtools index (needed for GATK4)
## run docker samtools
if [ ! -f ref.fa.fai ]; then
    echo "Creating samtools index..."
    apptainer exec genehetx-rnaseq_v1.6.0.sif samtools faidx ref.fa
fi

##create a dict (needed for GATK4)
##you should specify the path to your picard.jar file or run docker
if [ ! -f ref.dict ]; then
    echo "Creating Picard dict..."
    apptainer exec genehetx-rnaseq_v1.6.0.sif java -jar /usr/local/bin/picard.jar  CreateSequenceDictionary R= ref.fa O= ref.dict
fi

##indexinf the vcf file (needed for  haplotypeCaller)
##you should specify the path to your gatk jar file of the specific region you use
if [ ! -f knowns_variants.vcf.idx ]; then
    echo "Creating GATK index for VCF file..."
    apptainer exec genehetx-rnaseq_v1.6.0.sif java -jar /usr/local/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar IndexFeatureFile -I knowns_variants.vcf
fi

echo "Done creating tool indexes."

echo "Parsing GTF file to extract exon information..."
##parsing a gtf file in order to get exon informations (needed for featureCount output analysis)
awk -F "\t" '$3 == "exon" { print $4"\t"$5"\t"$7"\t"$9 }' ref.gtf |awk '{for(i=5;i<=NF;i++){if($i~/^"ENSE/){a=$i}} print a, $1,$2,$3,$5,$15}'| sed 's/\"//g'|sed 's/\;//g'| sort -d | awk 'BEGIN {print "exon_id\tstart\tend\tstrand\tgene_id\tgene_name"} { print }' >Exon_gtf_info.tab

grep -P "\tgene\t" ref.gtf > ref.GeneLvlOnly.gtf
Rscript ${RNApipe}/PrePostScripts/procGTF.R ref.GeneLvlOnly.gtf refGeneID_ensembl_v${Ensemblv}
echo "Done parsing GTF file."

## Copy files in the new directory
cp transcriptom.fa $newDir/transcriptom.fa
cp ref.gtf $newDir/ref.gtf
cp ref.fa $newDir/ref.fa
cp kallisto_index $newDir/kallisto_index
cp ref.fa.fai $newDir/ref.fa.fai
cp ref.dict $newDir/ref.dict
cp knowns_variants.vcf  $newDir/knowns_variants.vcf
cp knowns_variants.vcf.idx $newDir/knowns_variants.vcf.idx
cp Exon_gtf_info.tab $newDir/Exon_gtf_info.tab

#VEP files
mkdir -p $newDir/VEP/
cp -r ${species} $newDir/VEP/.


##generate a STAR index
## run docker star
echo "Creating STAR index..."
apptainer exec genehetx-rnaseq_v1.6.0.sif STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $newDir --genomeFastaFiles $newDir/ref.fa  --sjdbOverhang 100 --sjdbGTFfile $newDir/ref.gtf  --genomeSAindexNbases 11

echo "Removing temporary files..."
## remove all files except the new directory
rm -f ${prefix}.${genome_version}.${Ensemblv}.chr.gtf.gz
rm -f ${prefix}.${genome_version}.dna.primary_assembly.fa.gz
rm -f ${prefix}.${genome_version}.cdna.all.fa.gz
rm -f ${species}_vep_${Ensemblv}_${genome_version}.tar.gz
rm -f ${species}.vcf.gz
rm -f 1000GENOMES-phase_3.vcf.gz
rm -f ref.gtf
rm -f ref.fa
rm -f transcriptom.fa
rm -f kallisto_index
rm -f ref.fa.fai
rm -f ref.dict
rm -f knowns_variants.vcf
rm -f knowns_variants.vcf.idx
rm -f Exon_gtf_info.tab
rm -f ref.GeneLvlOnly.gtf
rm -f refGeneID_ensembl_v${Ensemblv}.*
rm -rf ${species}

echo "DONE GENERATING REFERENCE FILES FOR ${species} in ${Ensemblv} - ${genome_version}. All files are in $newDir"
#copy slurm job
#cp ${0} $newDir/build_ref${Ensemblv}_${species}_${genome_version}.slurm
