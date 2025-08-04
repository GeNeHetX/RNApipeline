FROM ubuntu:latest
MAINTAINER LUCIE GOMES  <lucie.gomes@inserm.fr>

# Installer les dépendances et outils nécessaires
RUN apt-get update && \
    apt-get install -y \
    build-essential \
	openjdk-21-jdk \
    wget \
    zlib1g-dev \
    libbz2-dev \ 
    libcurl4-openssl-dev \ 
    liblzma-dev \
    g++ \
    git \
    perl \
    unzip \
    libncurses5-dev \
    libreadline-dev \
    vim-common \
    && rm -rf /var/lib/apt/lists/*

# Télécharger, extraire, compiler et installer STAR en une seule commande RUN
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz -O /tmp/STAR-2.7.11b.tar.gz && \
    tar -xzvf /tmp/STAR-2.7.11b.tar.gz -C /tmp && \
    cd /tmp/STAR-2.7.11b/source && \
    echo "Starting compilation..." && \
    make STAR && \
    if [ $? -ne 0 ]; then echo "STAR Make failed"; exit 1; fi && \
    cp /tmp/STAR-2.7.11b/source/STAR /usr/local/bin/ && \
    rm -rf /tmp/STAR-2.7.11b.tar.gz /tmp/STAR-2.7.11b && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Télécharger, extraire, compiler et installer Samtools en une seule commande RUN
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 -O /tmp/samtools-1.21.tar.bz2 && \
    tar -xjf /tmp/samtools-1.21.tar.bz2 -C /tmp && \
    cd /tmp/samtools-1.21 && \
    echo "Starting Samtools compilation..." && \
    make && \
    if [ $? -ne 0 ]; then echo "Samtools Make failed"; exit 1; fi && \
    cp /tmp/samtools-1.21/samtools /usr/local/bin/ && \
    rm -rf /tmp/samtools-1.21.tar.bz2 /tmp/samtools-1.21 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Télécharger et installer Picard
RUN wget -O /usr/local/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar && \
    chmod +x /usr/local/bin/picard.jar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Télécharger et installer GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip -O /tmp/gatk-4.6.1.0.zip && \
    unzip /tmp/gatk-4.6.1.0.zip -d /usr/local/ && \
    rm /tmp/gatk-4.6.1.0.zip && \
    ln -s /usr/local/gatk-4.6.1.0/gatk /usr/local/bin/gatk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Définir le répertoire de travail par défaut
WORKDIR /data