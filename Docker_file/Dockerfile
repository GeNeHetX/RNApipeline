FROM ubuntu:latest
MAINTAINER NASSIMA IMARAZENE  <nassima.imarazene@inserm.fr>
### Install jdk8
RUN apt-get update && apt-get install --yes openjdk-8-jre-headless \
	unzip\
	&& apt-get install -y libfindbin-libs-perl
	
###INSTALL FASTQC 
ENV FASTQC_url http://www.bioinformatics.babraham.ac.uk/projects/fastqc

WORKDIR /usr/local
ADD ${FASTQC_url}/fastqc_v0.11.9.zip /tmp/
RUN unzip /tmp/fastqc_v0.11.9.zip && sed -i 's/Xmx250m/Xmx2048m/' FastQC/fastqc && chmod 755 FastQC/fastqc
ENV PATH /usr/local/FastQC/:$PATH
CMD fastqc

###INSTALLING MULTIQC 

RUN apt-get update && apt-get install -y \
  wget \
  python3 \
  python3-dev \
  python3-pip 
  
RUN pip3 install "multiqc==1.12"

###installing star 
RUN apt-get install -y rna-star
##installing SUBREAD 
RUN apt-get install -y subread 

# Set WORKDIR to /data
RUN mkdir /data
WORKDIR /data
###install picard tools 
RUN wget -O picard.jar https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar
###install gatk
RUN wget  https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip
RUN unzip gatk-4.2.5.0.zip	
ARG DEBIAN_FRONTEND=noninteractive
ENV SA=Europe/Paris
RUN apt-get install -y samtools
RUN sed -i -e '/^assistive_technologies=/s/^/#/' /etc/java-*-openjdk/accessibility.properties
#install kallisto 
RUN apt-get update &&  apt-get -y install kallisto
# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*
