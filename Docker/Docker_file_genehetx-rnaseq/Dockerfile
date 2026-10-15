# Native ARM64 and x86_64 builds, with the same versions on both platforms.
FROM genehetx/genehetx-rnaseq:v1.6.1@sha256:8454e0ff2394009e84f424876128a9b19f5ab21b7bce2781032052fa8920d007 AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build autoconf automake libhdf5-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt
RUN curl -fL --retry 3 https://github.com/pachterlab/kallisto/archive/refs/tags/v0.51.0.tar.gz -o kallisto.tar.gz \
    && echo 'efeb0191c1a6a0d6de69111fb66f4bda51ff31fb40c513280f072bd44556f80d  kallisto.tar.gz' | sha256sum -c - \
    && tar -xzf kallisto.tar.gz \
    && cd kallisto-0.51.0/ext/htslib && autoheader && autoconf \
    && cd ../.. \
    && sed -i 's/-mno-avx2//g' ext/bifrost/CMakeLists.txt \
    && cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_HDF5=ON \
         -DCOMPILATION_ARCH=OFF -DENABLE_AVX2=OFF \
    && cmake --build build --target bifrost --parallel 4 \
    && cmake --build build --parallel 4
RUN curl -fL --retry 3 https://github.com/ShiLab-Bioinformatics/subread/archive/refs/tags/2.0.2.tar.gz -o subread.tar.gz \
    && echo '02a54526cc77265e2c805382a8b4fdb20402b496c3d0b5726aed82ea0427ff10  subread.tar.gz' | sha256sum -c - \
    && tar -xzf subread.tar.gz \
    && sed -i 's/-mtune=core2 //g' subread-2.0.2/src/Makefile.Linux \
    && make -C subread-2.0.2/src -f Makefile.Linux -j4 featureCounts
RUN curl -fL --retry 3 https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip -o fastqc.zip \
    && echo '5f4dba8780231a25a6b8e11ab2c238601920c9704caa5458d9de559575d58aa7  fastqc.zip' | sha256sum -c - \
    && unzip -q fastqc.zip && chmod +x FastQC/fastqc

FROM genehetx/genehetx-rnaseq:v1.6.1@sha256:8454e0ff2394009e84f424876128a9b19f5ab21b7bce2781032052fa8920d007
COPY --from=build /opt/kallisto-0.51.0/build/src/kallisto /usr/local/bin/kallisto
COPY --from=build /opt/subread-2.0.2/src/featureCounts /usr/local/bin/featureCounts
COPY --from=build /opt/FastQC /usr/local/FastQC
RUN apt-get update && apt-get install -y --no-install-recommends libhdf5-103-1 \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir multiqc==1.25.1 pandas==2.2.3 numpy==1.26.4 tqdm==4.67.1
# A successful build verifies the binaries, not just the manifest architecture.
RUN STAR --version | grep -Fx '2.7.10a' \
    && kallisto version | grep -F 'version 0.51.0' \
    && featureCounts -v 2>&1 | grep -F 'v2.0.2' \
    && fastqc --version | grep -F 'v0.12.1' \
    && multiqc --version | grep -F 'version 1.25.1' \
    && samtools --version | head -1 | grep -Fx 'samtools 1.13'
LABEL org.opencontainers.image.version="1.7.1"
WORKDIR /data
