
FROM ensemblorg/ensembl-vep:latest

 
USER vep

RUN perl INSTALL.pl \
        --AUTO cfp \
        --NO_UPDATE \
        --ASSEMBLY GRCh38 \
        --PLUGINSDIR /opt/vep/.vep/Plugins/ \
        --CACHEDIR /opt/vep/.vep/ \
        --PLUGINS all \
        --SPECIES homo_sapiens  
