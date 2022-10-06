FROM ubuntu:latest
MAINTAINER REMY NICOLLE  <remy.nicolle@inserm.fr>
### Install jdk8
RUN apt-get update && \
 	apt-get install --yes openjdk-8-jre-headless \
	unzip libfindbin-libs-perl default-jre \
	wget \
	python3 \
	python3-dev \
	python3-pip \
	rna-star \
	subread \
	samtools \
	kallisto \
  curl \
	&& 	pip3 install "multiqc==1.12" \
	&&  curl -s https://get.nextflow.io | bash \
	&& 	mkdir /data


RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-404.0.0-linux-x86_64.tar.gz \
  && tar -xf google-cloud-cli-404.0.0-linux-x86_64.tar.gz \
  && ./google-cloud-sdk/install.sh -q   --rc-path /root/.bashrc

###INSTALL FASTQC
ENV SA=Europe/Paris

ENV FASTQC_url http://www.bioinformatics.babraham.ac.uk/projects/fastqc
ARG DEBIAN_FRONTEND=noninteractive


WORKDIR /usr/local
ADD ${FASTQC_url}/fastqc_v0.11.9.zip /tmp/
RUN unzip /tmp/fastqc_v0.11.9.zip && sed -i 's/Xmx250m/Xmx2048m/' FastQC/fastqc && chmod 755 FastQC/fastqc
ENV PATH /usr/local/FastQC/:$PATH



WORKDIR /data
###install picard tools
RUN wget -O picard.jar https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar \
   && wget  https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip  \
   && unzip gatk-4.2.5.0.zip \
   && sed -i -e '/^assistive_technologies=/s/^/#/' /etc/java-*-openjdk/accessibility.properties

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*


# docker build --rm . -t genehetx/genehetx-rnaseq:latest
# docker push genehetx/genehetx-rnaseq:latest
