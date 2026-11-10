# Build with --build-arg ENSEMBL_RELEASE=107 and tag genehetx/genehetx-vep:107.
# No cache or FASTA lives in this image; buildReference.nf supplies both.
FROM ubuntu:22.04
ARG ENSEMBL_RELEASE=107
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    ca-certificates curl git unzip perl python3 tabix \
    libdbi-perl libdbd-mysql-perl libjson-perl libtext-csv-perl \
    libset-intervaltree-perl libperlio-gzip-perl libwww-perl libarchive-zip-perl \
    libconfig-inifiles-perl libbio-perl-perl libbio-db-hts-perl \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/vep
RUN git clone --depth 1 --branch release/${ENSEMBL_RELEASE} https://github.com/Ensembl/ensembl-vep.git . \
    && perl INSTALL.pl --AUTO a --NO_UPDATE --NO_HTSLIB --NO_BIOPERL --NO_TEST \
    && git rev-parse HEAD > source_revision.txt
ENV PERL5LIB=/opt/vep:/opt/vep/modules:/opt/vep/Bio
ENV PATH=/opt/vep:$PATH
RUN vep --help | grep -F 'ensembl-vep'
LABEL org.opencontainers.image.title="RNApipeline VEP"
WORKDIR /data
