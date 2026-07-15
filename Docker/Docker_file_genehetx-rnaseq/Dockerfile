FROM ubuntu:22.04
MAINTAINER REMY NICOLLE  <remy.nicolle@inserm.fr>
ARG TARGETARCH
ARG NEXTFLOW_VERSION=24.04.4
### Install jdk8
RUN apt-get update && \
 	apt-get install --yes openjdk-8-jre-headless \
	openjdk-17-jre-headless \
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
	&& case "${TARGETARCH}" in amd64) JAVA_ARCH=amd64 ;; arm64) JAVA_ARCH=arm64 ;; *) echo "Unsupported TARGETARCH: ${TARGETARCH}" >&2; exit 1 ;; esac \
	&& update-alternatives --set java "/usr/lib/jvm/java-17-openjdk-${JAVA_ARCH}/bin/java" \
	&& curl -fL "https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow" \
	    -o /usr/local/bin/nextflow \
	&& chmod 0755 /usr/local/bin/nextflow \
	&& 	mkdir /data


RUN case "${TARGETARCH}" in \
      amd64) GCLOUD_ARCH=x86_64 ;; \
      arm64) GCLOUD_ARCH=arm ;; \
      *) echo "Unsupported TARGETARCH: ${TARGETARCH}" >&2; exit 1 ;; \
    esac \
  && curl -fL "https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-${GCLOUD_ARCH}.tar.gz" \
       -o /tmp/google-cloud-cli.tar.gz \
  && tar -xf /tmp/google-cloud-cli.tar.gz \
  && ./google-cloud-sdk/install.sh -q   --rc-path /root/.bashrc

###INSTALL FASTQC
ENV SA=Europe/Paris

ENV FASTQC_URL=https://www.bioinformatics.babraham.ac.uk/projects/fastqc
ARG DEBIAN_FRONTEND=noninteractive


WORKDIR /usr/local
RUN curl -fL --retry 3 \
    "${FASTQC_URL}/fastqc_v0.11.9.zip" \
    -o /tmp/fastqc_v0.11.9.zip
RUN unzip /tmp/fastqc_v0.11.9.zip && sed -i 's/Xmx250m/Xmx2048m/' FastQC/fastqc && chmod 755 FastQC/fastqc
ENV PATH=/usr/local/FastQC/:$PATH



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
