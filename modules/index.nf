process DownloadReference {
    input:
    val sources

    output:
    tuple path('ref.fa'), path('ref.gtf'), path('cdna.fa'), path('VEP'), path('sites.vcf.gz'), path('sources.json'), emit: files

    script:
    def source_json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(sources))
    """
    set -euo pipefail
    wget --tries=3 --timeout=60 -nv -O ref.fa.gz '${sources.fasta}'
    wget --tries=3 --timeout=60 -nv -O ref.gtf.gz '${sources.gtf}'
    wget --tries=3 --timeout=60 -nv -O cdna.fa.gz '${sources.cdna}'
    wget --tries=3 --timeout=60 -nv -O vep_cache.tar.gz '${sources.vep}'
    wget --tries=3 --timeout=60 -nv -O sites.vcf.gz '${sources.somalier}'
    echo '${params.somalier_sites_sha256}  sites.vcf.gz' | sha256sum -c -
    for asset in ref.fa ref.gtf cdna.fa; do
        gzip -t "\${asset}.gz"
        gzip -cd "\${asset}.gz" > "\${asset}"
        test -s "\${asset}"
    done
    mkdir VEP
    tar -xzf vep_cache.tar.gz -C VEP
    cat > sources.json <<'SOURCES'
${source_json}
SOURCES
    python3 - <<'PY'
import hashlib,json
from pathlib import Path
p=Path('sources.json'); sources=json.loads(p.read_text())
for key,name in [('fasta','ref.fa.gz'),('gtf','ref.gtf.gz'),('cdna','cdna.fa.gz'),('vep','vep_cache.tar.gz'),('somalier','sites.vcf.gz')]:
    h=hashlib.sha256()
    with open(name,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    sources[key]={'url':sources[key],'download_sha256':h.hexdigest()}
p.write_text(json.dumps(sources,indent=2)+'\\n')
PY
    rm ref.fa.gz ref.gtf.gz cdna.fa.gz vep_cache.tar.gz
    """
}

process buildref {
    cpus { params.reference_build_threads }
    memory { params.reference_build_memory }

    input:
    tuple path(fasta), path(gtf), path(cdna), path(cache), path(panel), path(sources)
    path catalogue
    path overrides
    path metadata_script
    path hotspot_script

    output:
    path 'refdata', emit: reference

    script:
    def override_arg = overrides ? "--transcript-overrides ${overrides}" : ''
    """
    set -euo pipefail
    mkdir -p refdata/hotspots
    cp -L ${fasta} refdata/ref.fa
    cp -L ${gtf} refdata/ref.gtf
    cp -L ${cdna} refdata/cdna.fa
    cp -RL ${cache} refdata/VEP
    cp -L ${sources} refdata/sources.json
    cp -L ${catalogue} refdata/hotspots/catalogue.tsv
    samtools faidx refdata/ref.fa
    python3 ${metadata_script} prepare --root refdata --release ${params.ensembl_release} --panel ${panel}
    bcftools sort refdata/somalier/sites.vcf -Oz -o refdata/somalier/sites.vcf.gz
    bcftools index -t refdata/somalier/sites.vcf.gz
    rm refdata/somalier/sites.vcf
    python3 ${hotspot_script} build --catalogue refdata/hotspots/catalogue.tsv \
        --gtf refdata/ref.gtf --fasta refdata/ref.fa --outdir refdata/hotspots ${override_arg}
    kallisto index -i refdata/kallisto_index refdata/cdna.fa
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
        --genomeDir refdata --genomeFastaFiles refdata/ref.fa \
        --sjdbOverhang 100 --sjdbGTFfile refdata/ref.gtf --genomeSAindexNbases 11
    """
}

process ValidateHotspots {
    cpus 1
    memory { params.vep_memory }

    input:
    path reference
    path hotspot_script

    output:
    path 'validation.json', emit: validation

    script:
    """
    set -euo pipefail
    vep --help | grep -E 'ensembl-vep[[:space:]]*:[[:space:]]*${params.ensembl_release}[.]'
    vep --input_file ${reference}/hotspots/candidate_variants.vcf --format vcf \
        --output_file candidates.vep.json --json --offline --cache --no_stats \
        --cache_version ${params.ensembl_release} --species homo_sapiens --assembly GRCh38 \
        --dir_cache ${reference}/VEP --fasta ${reference}/ref.fa --hgvs --protein --transcript_version
    python3 ${hotspot_script} validate --alleles ${reference}/hotspots/candidate_alleles.json \
        --vep candidates.vep.json --output validation.json
    """
}

process PublishReference {
    publishDir path: params.reference_build_root, mode: 'copy', overwrite: false,
        saveAs: { name -> name == 'bundle' ? params.reference_name : null }

    input:
    path reference
    path validation
    path metadata_script

    output:
    path 'bundle'

    script:
    def containers_json = groovy.json.JsonOutput.toJson(params.containers)
    """
    set -euo pipefail
    python3 ${metadata_script} assemble --root bundle --release ${params.ensembl_release} \
        --reference ${reference} --validation ${validation}
    python3 ${metadata_script} manifest --root bundle --release ${params.ensembl_release} \
        --sources bundle/sources.json --pipeline-version ${params.pipeline_version} --containers '${containers_json}'
    """
}
