nextflow.enable.dsl=2

/*
 * Infrastructure-neutral Ensembl reference build.
 *
 * Site/project configuration supplies reference_root, containers, executor,
 * queue, and resources. Each process
 * publishes its own completed artifacts to the shared staging directory.
 * Finalization validates them and performs only a same-filesystem rename.
 */

params.reference_root = '/ref'
params.reference_id = 'ensembl_v107_GRCh38'
params.ensembl_release = 107
params.reference_build_owner = System.getenv('NF_RUN_NAME') ?: 'manual'
params.force = false
params.check_only = false
params.reference_kallisto_before_script = ''
params.genome_url = null
params.gtf_url = null
params.cdna_url = null
params.known_vcf_url = null
params.vep_cache_url = null
params.reference_picard_jar = '/data/picard.jar'
params.reference_gatk_jar = '/data/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar'
params.reference_containers = [
    utility: 'docker://ubuntu:22.04',
    core: 'docker://genehetx/genehetx-rnaseq:v1.6.1',
    r: 'docker://rocker/r-ver:4.3.3',
    kallisto: 'docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0'
]

def release = params.ensembl_release.toString()
def genome_url = params.genome_url ?: ('https://ftp.ensembl.org/pub/release-' + release + '/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
def gtf_url = params.gtf_url ?: ('https://ftp.ensembl.org/pub/release-' + release + '/gtf/homo_sapiens/Homo_sapiens.GRCh38.' + release + '.chr.gtf.gz')
def cdna_url = params.cdna_url ?: ('https://ftp.ensembl.org/pub/release-' + release + '/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz')
def known_vcf_url = params.known_vcf_url ?: ('https://ftp.ensembl.org/pub/release-' + release + '/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz')
def vep_cache_url = params.vep_cache_url ?: ('https://ftp.ensembl.org/pub/release-' + release + '/variation/indexed_vep_cache/homo_sapiens_vep_' + release + '_GRCh38.tar.gz')

process INIT_REFERENCE_STAGE {
    tag 'reference stage'
    label 'reference'
    executor 'local'
    container null
    scratch false

    input:
    val stage_root
    val final_root
    val force
    val owner

    output:
    path 'reference-stage.ready', emit: ready

    script:
    """
    set -euo pipefail
    if [[ -e '${final_root}/reference_manifest.json' && '${force}' != 'true' ]]; then
        echo "Reference already finalized: ${final_root}"
        touch reference-stage.ready
        exit 0
    fi
    if [[ '${force}' == 'true' ]]; then
        rm -rf -- '${stage_root}'
    elif [[ -e '${stage_root}/.build.active' && -s '${stage_root}/.build.owner' \
            && "\$(cat '${stage_root}/.build.owner')" != '${owner}' ]]; then
        echo "Reference staging directory is owned by another run: ${stage_root}" >&2
        echo "Resume that run or use --force only for an intentional replacement." >&2
        exit 1
    elif [[ -e '${stage_root}' && ! -e '${stage_root}/.build.active' ]]; then
        echo "Refusing to reuse an unmarked staging directory: ${stage_root}" >&2
        exit 1
    fi
    mkdir -p '${stage_root}'
    touch '${stage_root}/.build.active'
    printf '%s\n' '${owner}' > '${stage_root}/.build.owner'
    touch reference-stage.ready
    """
}

process PREPARE_REFERENCE_INPUTS {
    tag 'reference inputs'
    label 'reference'
    container params.reference_containers.utility
    cpus 4
    memory '16 GB'
    time '12h'
    scratch false

    input:
    path stage_ready
    val stage_root

    output:
    path 'ref.fa', emit: ref_fa
    path 'ref.gtf', emit: ref_gtf
    path 'transcriptom.fa', emit: transcriptome
    path 'cdna.fa', emit: cdna
    path 'knowns_variants.vcf', emit: known_vcf
    path 'VEP', emit: vep
    path 'reference-inputs.ready', emit: ready

    script:
    """
    set -euo pipefail
    mkdir -p sources VEP
    if [[ -e '${stage_root}/.inputs.complete' ]]; then
        cp -a '${stage_root}/sources' '${stage_root}/VEP' \
            '${stage_root}/ref.fa' '${stage_root}/ref.gtf' \
            '${stage_root}/transcriptom.fa' '${stage_root}/cdna.fa' \
            '${stage_root}/knowns_variants.vcf' .
    else
    fetch() {
        url="\$1"; destination="\$2"; partial="\${destination}.part"
        if [[ ! -s "\$destination" ]]; then
            curl --fail --location --retry 5 --retry-delay 5 --continue-at - \
                --output "\$partial" "\$url"
            mv -- "\$partial" "\$destination"
        fi
    }
    fetch '${genome_url}' sources/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    fetch '${gtf_url}' sources/Homo_sapiens.GRCh38.${release}.chr.gtf.gz
    fetch '${cdna_url}' sources/Homo_sapiens.GRCh38.cdna.all.fa.gz
    fetch '${known_vcf_url}' sources/1000GENOMES-phase_3.vcf.gz
    fetch '${vep_cache_url}' sources/homo_sapiens_vep_${release}_GRCh38.tar.gz
    gzip -t sources/*.gz
    tar -tzf sources/homo_sapiens_vep_${release}_GRCh38.tar.gz >/dev/null
    gzip -dc sources/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > ref.fa
    gzip -dc sources/Homo_sapiens.GRCh38.${release}.chr.gtf.gz > ref.gtf
    gzip -dc sources/Homo_sapiens.GRCh38.cdna.all.fa.gz > transcriptom.fa
    gzip -dc sources/1000GENOMES-phase_3.vcf.gz > knowns_variants.vcf
    tar -xzf sources/homo_sapiens_vep_${release}_GRCh38.tar.gz -C VEP --no-same-owner
    test -d VEP/homo_sapiens/${release}_GRCh38
    cp -f transcriptom.fa cdna.fa
    mkdir -p '${stage_root}/sources' '${stage_root}/VEP'
    fi
    cp -a sources VEP ref.fa ref.gtf transcriptom.fa cdna.fa knowns_variants.vcf \
        '${stage_root}/'
    touch '${stage_root}/.inputs.complete' reference-inputs.ready
    """
}

process BUILD_SEQUENCE_INDEXES {
    tag 'sequence indexes'
    label 'reference'
    container params.reference_containers.core
    cpus 16
    memory '32 GB'
    time '8h'
    scratch false

    input:
    path ref_fa
    path inputs_ready
    path stage_ready
    val stage_root

    output:
    path 'ref.fa.fai', emit: fai
    path 'ref.dict', emit: dict
    path 'sequence-indexes.ready', emit: ready

    script:
    """
    set -euo pipefail
    samtools faidx ref.fa
    java -jar '${params.reference_picard_jar}' CreateSequenceDictionary R=ref.fa O=ref.dict
    cp -f ref.fa.fai ref.dict '${stage_root}/'
    touch '${stage_root}/.sequence-indexes.complete' sequence-indexes.ready
    """
}

process BUILD_VARIANT_INDEX {
    tag 'variant index'
    label 'reference'
    container params.reference_containers.core
    cpus 8
    memory '16 GB'
    time '4h'
    scratch false

    input:
    path known_vcf
    path inputs_ready
    path stage_ready
    val stage_root

    output:
    path 'knowns_variants.vcf.idx', emit: index
    path 'variant-index.ready', emit: ready

    script:
    """
    set -euo pipefail
    java -jar '${params.reference_gatk_jar}' IndexFeatureFile -I knowns_variants.vcf
    cp -f knowns_variants.vcf.idx '${stage_root}/'
    touch '${stage_root}/.variant-index.complete' variant-index.ready
    """
}

process BUILD_KALLISTO_INDEX {
    tag 'Kallisto index'
    label 'reference'
    container params.reference_containers.kallisto
    beforeScript params.reference_kallisto_before_script
    cpus 16
    memory '32 GB'
    time '8h'
    scratch false

    input:
    path transcriptome
    path inputs_ready
    path stage_ready
    val stage_root

    output:
    path 'kalliso_index', emit: index
    path 'kallisto_index', emit: compatible_index
    path 'kallisto-index.ready', emit: ready

    script:
    """
    set -euo pipefail
    kallisto index -i kalliso_index transcriptom.fa
    ln kalliso_index kallisto_index
    cp -f kalliso_index kallisto_index '${stage_root}/'
    touch '${stage_root}/.kallisto-index.complete' kallisto-index.ready
    """
}

process BUILD_GTF_ARTIFACTS {
    tag 'GTF artifacts'
    label 'reference'
    container params.reference_containers.r
    cpus 8
    memory '16 GB'
    time '4h'
    scratch false

    input:
    path ref_gtf
    path inputs_ready
    path stage_ready
    val stage_root

    output:
    path 'gtf', emit: artifacts
    path 'gtf-artifacts.ready', emit: ready

    script:
    """
    set -euo pipefail
    mkdir -p gtf
    awk -F '\\t' '\$3 == "gene" { print }' ref.gtf > gtf/ref.GeneLvlOnly.gtf
    awk -F '\\t' '\$3 == "exon" { print \$1 "\\t" (\$4 - 1) "\\t" \$5 }' ref.gtf \
        | sort -k1,1 -k2,2n -k3,3n -u > gtf/ref.ExonsOnly.bed
    awk -F '\\t' '\$3 == "exon" { print \$4 "\\t" \$5 "\\t" \$7 "\\t" \$9 }' ref.gtf \
        | awk '{for(i=5;i<=NF;i++){if(\$i~/^"ENSE/){a=\$i}} print a, \$1,\$2,\$3,\$5,\$15}' \
        | sed 's/"//g;s/;//g' | sort -d \
        | awk 'BEGIN {print "exon_id\\tstart\\tend\\tstrand\\tgene_id\\tgene_name"} { print }' \
        > gtf/Exon_gtf_info.tab
    cat > gtf/procGTF.R <<'R_SCRIPT'
    args=commandArgs(trailingOnly=TRUE)
    library(parallel)
    gtf=read.delim(args[1],sep="\\t",as.is=T,header=F,comment.char="#")
    getmeta=function(x){strsplit(x[,9],";| ")}
    getvalue=function(x,field){unlist(mclapply(getmeta(x),function(y){
      if(field %in% y) y[which(y==field)+1] else NA
    }))}
    geneTab=gtf[which(gtf[,3]=="gene"),]
    geneTab$GeneID=getvalue(geneTab,"gene_id")
    geneTab$GeneName=getvalue(geneTab,"gene_name")
    geneTab$GeneName[is.na(geneTab$GeneName)]=geneTab$GeneID[is.na(geneTab$GeneName)]
    geneTab$biotype=getvalue(geneTab,"gene_biotype")
    geneTab=unique(geneTab); rownames(geneTab)=geneTab$GeneID
    geneTab=unique(geneTab[,-c(9,6,2,3,8)])
    colnames(geneTab)[1:4]=c("seqname","start","end","strand")
    saveRDS(geneTab,file=paste0(args[2],".rds"))
    write.table(geneTab,file=paste0(args[2],".tsv"),quote=F,sep="\\t")
    R_SCRIPT
    Rscript gtf/procGTF.R gtf/ref.GeneLvlOnly.gtf gtf/refGeneID_ensembl_v${release}
    cp -a gtf/. '${stage_root}/'
    touch '${stage_root}/.gtf-artifacts.complete' gtf-artifacts.ready
    """
}

process BUILD_STAR_INDEX {
    tag 'STAR index'
    label 'reference'
    container params.reference_containers.core
    cpus 16
    memory '64 GB'
    time '24h'
    scratch false

    input:
    path ref_fa
    path ref_gtf
    path inputs_ready
    path stage_ready
    val stage_root

    output:
    path 'star', emit: index
    path 'star-index.ready', emit: ready

    script:
    """
    set -euo pipefail
    mkdir -p star
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
        --genomeDir star --genomeFastaFiles ref.fa --sjdbOverhang 100 \
        --sjdbGTFfile ref.gtf --genomeSAindexNbases 11
    cp -a star/. '${stage_root}/'
    touch '${stage_root}/.star-indexes.complete' star-index.ready
    """
}

process FINALIZE_REFERENCE {
    tag 'publish reference'
    label 'reference'
    executor 'local'
    container null
    scratch false
    cache false

    input:
    path stage_ready
    path inputs_ready
    path sequence_ready
    path variant_ready
    path kallisto_ready
    path gtf_ready
    path star_ready
    path sequence_fai
    path sequence_dict
    path variant_index
    path kallisto_index
    path kallisto_compatible_index
    path gtf_artifacts
    path star_index
    path vep
    val stage_root
    val final_root
    val force

    output:
    path 'reference-build.done'

    script:
    """
    set -euo pipefail
    mkdir -p '${stage_root}'
    restore_file() {
        source="\$1"; destination="\$2"
        if [[ ! -s "\$destination" ]]; then
            cp -f "\$source" "\$destination"
        fi
    }
    restore_file sequence_fai '${stage_root}/ref.fa.fai'
    restore_file sequence_dict '${stage_root}/ref.dict'
    restore_file variant_index '${stage_root}/knowns_variants.vcf.idx'
    restore_file kallisto_index '${stage_root}/kalliso_index'
    restore_file kallisto_compatible_index '${stage_root}/kallisto_index'
    for file in gtf_artifacts/*; do
        name="\$(basename "\$file")"
        [[ -s '${stage_root}/'\$name ]] || cp -f "\$file" '${stage_root}/'\$name
    done
    for file in star_index/*; do
        name="\$(basename "\$file")"
        [[ -s '${stage_root}/'\$name ]] || cp -f "\$file" '${stage_root}/'\$name
    done
    if [[ ! -d '${stage_root}/VEP/homo_sapiens/${release}_GRCh38' ]]; then
        rm -rf '${stage_root}/VEP'
        cp -a vep '${stage_root}/VEP'
    fi
    required=(
      ref.fa ref.fa.fai ref.dict ref.gtf transcriptom.fa cdna.fa
      knowns_variants.vcf knowns_variants.vcf.idx kalliso_index kallisto_index
      ref.GeneLvlOnly.gtf ref.ExonsOnly.bed Exon_gtf_info.tab
      refGeneID_ensembl_v${release}.rds refGeneID_ensembl_v${release}.tsv
      Genome SA SAindex chrLength.txt chrName.txt chrNameLength.txt chrStart.txt
      genomeParameters.txt sjdbInfo.txt sjdbList.fromGTF.out.tab sjdbList.out.tab
      Log.out exonGeTrInfo.tab exonInfo.tab geneInfo.tab transcriptInfo.tab
    )
    for item in "\${required[@]}"; do
        test -s '${stage_root}/'\$item || {
            echo "Missing reference artifact: \$item" >&2
            exit 1
        }
    done
    test -d '${stage_root}/VEP/homo_sapiens/${release}_GRCh38'

    if [[ -e '${final_root}/reference_manifest.json' && '${force}' != 'true' ]]; then
        echo "Reference already finalized: ${final_root}"
        touch reference-build.done
        exit 0
    fi

    python3 - '${stage_root}' <<'PY'
import json, sys
from pathlib import Path
stage = Path(sys.argv[1])
required = [
    "ref.fa", "ref.fa.fai", "ref.dict", "ref.gtf", "transcriptom.fa",
    "cdna.fa", "knowns_variants.vcf", "knowns_variants.vcf.idx",
    "kalliso_index", "kallisto_index", "ref.GeneLvlOnly.gtf",
    "ref.ExonsOnly.bed", "Exon_gtf_info.tab",
    "refGeneID_ensembl_v${release}.rds", "refGeneID_ensembl_v${release}.tsv",
    "Genome", "SA", "SAindex",
]
manifest = {
    "manifest_schema_version": 2,
    "reference_id": "${params.reference_id}",
    "species": "homo_sapiens",
    "assembly": "GRCh38",
    "ensembl_release": ${release},
    "source_urls": {
        "genome": "${genome_url}", "gtf": "${gtf_url}",
        "cdna": "${cdna_url}", "known_vcf": "${known_vcf_url}",
        "vep_cache": "${vep_cache_url}",
    },
    "files": [{"path": n, "size_bytes": (stage / n).stat().st_size}
              for n in required],
    "complete": True,
}
(stage / "reference_manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\\n"
)
PY

    if [[ -e '${final_root}' ]]; then
        [[ '${force}' == 'true' ]] || {
            echo "Refusing to replace existing reference: ${final_root}" >&2
            exit 1
        }
        backup='${stage_root}.previous.\$(date -u +%Y%m%dT%H%M%SZ)'
        mv -- '${final_root}' "\$backup"
    fi
    rm -f '${stage_root}/.build.active' '${stage_root}/.build.owner' '${stage_root}/.inputs.complete' \
        '${stage_root}/.sequence-indexes.complete' '${stage_root}/.variant-index.complete' \
        '${stage_root}/.kallisto-index.complete' '${stage_root}/.gtf-artifacts.complete' \
        '${stage_root}/.star-indexes.complete'
    mv -- '${stage_root}' '${final_root}'
    test -s '${final_root}/reference_manifest.json'
    touch reference-build.done
    """
}

process CHECK_REFERENCE {
    tag 'check reference'
    label 'reference'
    executor 'local'
    container null
    scratch false
    cache false

    input:
    val final_root

    output:
    path 'reference-preflight.ready'

    script:
    """
    set -euo pipefail
    test -s '${final_root}/reference_manifest.json'
    python3 - '${final_root}/reference_manifest.json' <<'PY'
import json, sys
from pathlib import Path
manifest = json.loads(Path(sys.argv[1]).read_text())
assert manifest.get("complete") is True
assert manifest.get("ensembl_release") == ${release}
for item in manifest["files"]:
    path = Path(sys.argv[1]).parent / item["path"]
    assert path.is_file() and path.stat().st_size >= item["size_bytes"], path
PY
    test -d '${final_root}/VEP/homo_sapiens/${release}_GRCh38'
    touch reference-preflight.ready
    """
}

workflow {
    def run_name = System.getenv('NF_RUN_NAME') ?: 'manual'
    final_root = params.reference_root + '/' + params.reference_id
    stage_root = params.reference_root + '/.staging/' + run_name + '.active'

    if (params.check_only) {
        CHECK_REFERENCE(final_root)
    } else {
        stage = INIT_REFERENCE_STAGE(stage_root, final_root, params.force, params.reference_build_owner)
        inputs = PREPARE_REFERENCE_INPUTS(stage.ready, stage_root)
        sequence = BUILD_SEQUENCE_INDEXES(inputs.ref_fa, inputs.ready, stage.ready, stage_root)
        variants = BUILD_VARIANT_INDEX(inputs.known_vcf, inputs.ready, stage.ready, stage_root)
        kallisto = BUILD_KALLISTO_INDEX(inputs.transcriptome, inputs.ready, stage.ready, stage_root)
        gtf = BUILD_GTF_ARTIFACTS(inputs.ref_gtf, inputs.ready, stage.ready, stage_root)
        star = BUILD_STAR_INDEX(inputs.ref_fa, inputs.ref_gtf, inputs.ready, stage.ready, stage_root)
        FINALIZE_REFERENCE(
            stage.ready, inputs.ready, sequence.ready, variants.ready, kallisto.ready,
            gtf.ready, star.ready, sequence.fai, sequence.dict, variants.index,
            kallisto.index, kallisto.compatible_index, gtf.artifacts, star.index,
            inputs.vep, stage_root, final_root, params.force
        )
    }
}
