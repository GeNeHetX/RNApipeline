#!/bin/bash
## each user can ch{
mkdir ensembl_v$Ensembel_ver_GRCh38_p13} &&\
chmod +rwx ensembl_v${Ensembel_ver_GRCh38_p13}

##downoalding the gtf,fasta and cdna files 
wget http://ftp.ensembl.org/pub/release-${Ensembel_ver}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${Ensembel_ver}.chr.gtf.gz
wget http://ftp.ensembl.org/pub/release-${Ensembel_ver}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-${Ensembel_ver}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

##downoalding the vcf file 
wget http://ftp.ensembl.org/pub/release-${Ensembel_ver}/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz

	 
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ref.fa
	
gunzip -c Homo_sapiens.GRCh38.${Ensembel_ver}.chr.gtf.gz >ref.gtf
gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz >transcriptom.fa
	
	
##create an  samtools index (needed for GATK4)
samtools faidx ref.fa

##kallisto index 
kallisto index -i kalliso_index transcriptom.fa

	
##create a dict (needed for GATK4)
##you should specify the path to your picard.jar file 
java -jar picard.jar  CreateSequenceDictionary R= ref.fa O= ref.dict
	
##unzip the vcf file 
	
gunzip -c 1000GENOMES-phase_3.vcf.gz > knowns_variants.vcf
##indexinf the vcf file (needed for  haplotypeCaller) 
##you should specify the path to your  gatk-package-4.2.5.0-local.jar file 
java -jar gatk-package-4.2.5.0-local.jar IndexFeatureFile -I knowns_variants.vcf 

##parsing a gtf file in order to get exon informations (needed for featureCount output analysis)
awk -F "\t" '$3 == "exon" { print $4"\t"$5"\t"$7"\t"$9 }' ref.gtf |awk '{for(i=5;i<=NF;i++){if($i~/^"ENSE/){a=$i}} print a, $1,$2,$3,$5,$15}'| sed 's/\"//g'|sed 's/\;//g'| sort -d | awk 'BEGIN {print "exon_id\tstart\tend\tstrand\tgene_id\tgene_name"} { print }' >Exon_gtf_info.tab

mv ref.gtf ensembl_v${Ensembel_ver}_GRCh38_p13/ref.gtf
mv ref.fa.fai ensembl_v${Ensembel_ver}_GRCh38_p13/ref.fa.fai
mv ref.dict ensembl_v${Ensembel_ver}_GRCh38_p13/ref.dict
mv knowns_variants.vcf  ensembl_v${Ensembel_ver}_GRCh38_p13/knowns_variants.vcf 
mv knowns_variants.vcf.idx ensembl_v${Ensembel_ver}_GRCh38_p13/knowns_variants.vcf.idx
mv ref.fa ensembl_v${Ensembel_ver}_GRCh38_p13/ref.fa
mv Exon_gtf_info.tab ensembl_v${Ensembel_ver}_GRCh38_p13/Exon_gtf_info.tab
mv kalliso_index ensembl_v${Ensembel_ver}_GRCh38_p13/kalliso_index
cp ref_build.sh ensembl_v${Ensembel_ver}_GRCh38_p13/ref_data.sh


##generate a STAR index 
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ensembl_v${Ensembel_ver}_GRCh38_p13 --genomeFastaFiles ensembl_v${Ensembel_ver}_GRCh38_p13/ref.fa  --sjdbOverhang 100 --sjdbGTFfile ensembl_v${Ensembel_ver}_GRCh38_p13/ref.gtf  --genomeSAindexNbases 11 
		
