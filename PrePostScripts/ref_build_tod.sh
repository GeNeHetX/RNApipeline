#!/usr/bin/env bash
#SBATCH --job-name=rnapipeline-ref107
#SBATCH --partition=pam_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/tmp/rnapipeline-ref-build-%j.out
#SBATCH --error=/tmp/rnapipeline-ref-build-%j.err

set -Eeuo pipefail
IFS=$'\n\t'

# Build the complete precomputed reference used by RNApipeline on TOD/PAM.
# The build is staged on local Slurm scratch and published atomically to the
# direct CephFS reference mount.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REFERENCE_ID="ensembl_v107_GRCh38"
REF_ROOT="/ref"
WORK_ROOT="${SLURM_TMPDIR:-/srv/slurm/scratch}"
CHECK_ONLY=0
FORCE=0
KEEP_WORK=0
SKIP_MOUNT_CHECK=0

CORE_IMAGE_SOURCE="docker://genehetx/genehetx-rnaseq:v1.6.1"
KALLISTO_IMAGE_SOURCE="docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0"
R_IMAGE_SOURCE="docker://rocker/r-ver:4.3.3"
VEP_IMAGE_SOURCE="docker://quay.io/biocontainers/ensembl-vep:113.2--pl5321h2a3209d_0"

GENOME_URL="https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.chr.gtf.gz"
CDNA_URL="https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
KNOWN_VCF_URL="https://ftp.ensembl.org/pub/release-105/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
VEP_CACHE_URL="https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz"

usage() {
    cat <<'EOF'
Usage: ref_build_tod.sh [options]

Build the RNApipeline Ensembl 107 reference on a TOD/PAM Slurm worker.

Options:
  --check-only             Validate the mount, tools, and source URLs only.
  --force                  Replace an existing reference after preserving it
                           under /ref/.staging/.
  --keep-work              Keep the local scratch build directory on success.
  --ref-root PATH          Override /ref for testing or controlled staging.
  --work-root PATH         Override the local scratch build root.
  --skip-mount-check       Skip the CephFS mount check (tests only).
  -h, --help               Show this help.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

while (($# > 0)); do
    case "$1" in
        --check-only)
            CHECK_ONLY=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --keep-work)
            KEEP_WORK=1
            shift
            ;;
        --ref-root)
            (($# >= 2)) || die "--ref-root requires a path"
            REF_ROOT="$2"
            shift 2
            ;;
        --work-root)
            (($# >= 2)) || die "--work-root requires a path"
            WORK_ROOT="$2"
            shift 2
            ;;
        --skip-mount-check)
            SKIP_MOUNT_CHECK=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown option: $1"
            ;;
    esac
done

CURRENT_JOB="${SLURM_JOB_ID:-manual-$$}"
FINAL_DIR="${REF_ROOT%/}/${REFERENCE_ID}"
BUILD_DIR="${WORK_ROOT%/}/rnapipeline-ref-${REFERENCE_ID}"
LOCK_DIR="${BUILD_DIR}.lock"
STAGE_DIR="${REF_ROOT%/}/.staging/${REFERENCE_ID}.${CURRENT_JOB}"
IMAGE_DIR="${BUILD_DIR}/images"

CORE_SIF="${IMAGE_DIR}/genehetx-rnaseq-v1.6.1.sif"
KALLISTO_SIF="${IMAGE_DIR}/kallisto-0.51.1.sif"
R_SIF="${IMAGE_DIR}/r-ver-4.3.3.sif"
VEP_SIF="${IMAGE_DIR}/ensembl-vep-113.2.sif"

PICARD_JAR=""
GATK_JAR=""
APPTAINER_ARCH=""
LOCK_HELD=0

release_build_lock() {
    if ((LOCK_HELD == 1)); then
        rm -rf -- "$LOCK_DIR"
        LOCK_HELD=0
    fi
}

cleanup() {
    local rc=$?
    trap - EXIT

    if ((KEEP_WORK == 1 || rc != 0)); then
        if [[ -n "${BUILD_DIR:-}" && -d "${BUILD_DIR:-}" ]]; then
            log "Keeping scratch build directory: ${BUILD_DIR}"
        fi
        if ((rc != 0)) && [[ -n "${STAGE_DIR:-}" && -d "${STAGE_DIR:-}" ]]; then
            rm -rf -- "$STAGE_DIR"
        fi
    else
        if [[ -n "${BUILD_DIR:-}" && -d "${BUILD_DIR:-}" ]]; then
            rm -rf -- "$BUILD_DIR"
        fi
    fi

    release_build_lock
    exit "$rc"
}
trap cleanup EXIT

check_ref_mount() {
    [[ -d "$REF_ROOT" ]] || die "reference root does not exist: $REF_ROOT"

    if ((SKIP_MOUNT_CHECK == 0)); then
        require_cmd findmnt
        local fstype
        fstype="$(findmnt -T "$REF_ROOT" -n -o FSTYPE 2>/dev/null || true)"
        [[ "$fstype" == ceph ]] || die "$REF_ROOT is not mounted as CephFS (detected: ${fstype:-none})"
    fi

    [[ -w "$REF_ROOT" ]] || die "reference root is not writable: $REF_ROOT"
    local probe="${REF_ROOT%/}/.ref-build-write-test.${CURRENT_JOB}"
    touch "$probe"
    rm -f -- "$probe"
    mkdir -p -- "${REF_ROOT%/}/.staging"
}

check_source_urls() {
    local url
    for url in "$GENOME_URL" "$GTF_URL" "$CDNA_URL" "$KNOWN_VCF_URL" "$VEP_CACHE_URL"; do
        log "Checking source URL: $url"
        curl -fsSIL --retry 3 --connect-timeout 10 --max-time 60 "$url" >/dev/null \
            || die "source URL is not reachable: $url"
    done
}

acquire_build_lock() {
    if ! mkdir -- "$LOCK_DIR" 2>/dev/null; then
        local owner="unknown"
        if [[ -s "${LOCK_DIR}/job_id" ]]; then
            owner="$(<"${LOCK_DIR}/job_id")"
        fi
        die "scratch build is already locked by Slurm job ${owner}: ${BUILD_DIR}"
    fi
    LOCK_HELD=1
    printf '%s\n' "$CURRENT_JOB" >"${LOCK_DIR}/job_id"
}

host_apptainer_arch() {
    case "$(uname -m)" in
        x86_64)
            printf 'amd64\n'
            ;;
        aarch64|arm64)
            printf 'arm64\n'
            ;;
        *)
            die "unsupported host architecture: $(uname -m)"
            ;;
    esac
}

image_matches_arch() {
    local image="$1"
    local image_machine
    image_machine="$($APPTAINER_BIN exec "$image" uname -m 2>/dev/null)" || return 1

    case "${APPTAINER_ARCH}:${image_machine}" in
        amd64:x86_64|arm64:aarch64|arm64:arm64)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

ensure_image() {
    local target="$1"
    local source="$2"

    if [[ -s "$target" ]] && image_matches_arch "$target"; then
        log "Using cached image: $target"
        return
    fi

    if [[ -s "$target" ]]; then
        log "Removing cached image with the wrong or unusable architecture: $target"
        rm -f -- "$target"
    fi

    log "Pulling Apptainer image: $source"
    mkdir -p -- "$(dirname "$target")"
    "$APPTAINER_BIN" pull --arch "$APPTAINER_ARCH" "$target" "$source"
    [[ -s "$target" ]] || die "Apptainer pull did not produce: $target"
    image_matches_arch "$target" \
        || die "Apptainer image has the wrong or unusable architecture: $target (expected $APPTAINER_ARCH)"
}

apptainer_exec() {
    local image="$1"
    shift
    "$APPTAINER_BIN" exec \
        --bind "${BUILD_DIR}:/work" \
        --pwd /work \
        "$image" "$@"
}

core_exec() {
    apptainer_exec "$CORE_SIF" "$@"
}

kallisto_exec() {
    apptainer_exec "$KALLISTO_SIF" "$@"
}

r_exec() {
    "$APPTAINER_BIN" exec \
        --bind "${BUILD_DIR}:/work" \
        --bind "${SCRIPT_DIR}:/pipeline:ro" \
        --pwd /work \
        "$R_SIF" "$@"
}

vep_exec() {
    apptainer_exec "$VEP_SIF" "$@"
}

detect_core_path() {
    local candidate
    for candidate in "$@"; do
        [[ -n "$candidate" ]] || continue
        if core_exec test -f "$candidate" >/dev/null 2>&1; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

download_file() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"

    if [[ -s "$destination" ]]; then
        log "Using downloaded source: $destination"
        return
    fi

    log "Downloading: $url"
    curl --fail --location --retry 5 --retry-delay 5 --continue-at - \
        --output "$partial" "$url"
    mv -- "$partial" "$destination"
}

record_tool_versions() {
    core_exec STAR --version >"${BUILD_DIR}/tool-versions.core.txt" 2>&1
    core_exec samtools --version >"${BUILD_DIR}/tool-versions.samtools.txt" 2>&1
    core_exec java -version >"${BUILD_DIR}/tool-versions.java.txt" 2>&1
    kallisto_exec kallisto version >"${BUILD_DIR}/tool-versions.kallisto.txt" 2>&1
    r_exec R --version >"${BUILD_DIR}/tool-versions.r.txt" 2>&1
    vep_exec vep --version >"${BUILD_DIR}/tool-versions.vep.txt" 2>&1
}

build_reference() {
    mkdir -p -- "${BUILD_DIR}/sources" "${BUILD_DIR}/VEP" "${BUILD_DIR}/images" \
        "${BUILD_DIR}/apptainer-tmp"

    export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/srv/slurm/scratch/.apptainer-cache}"
    export APPTAINER_TMPDIR="${BUILD_DIR}/apptainer-tmp"
    mkdir -p -- "$APPTAINER_CACHEDIR"

    ensure_image "$CORE_SIF" "$CORE_IMAGE_SOURCE"
    ensure_image "$KALLISTO_SIF" "$KALLISTO_IMAGE_SOURCE"
    ensure_image "$R_SIF" "$R_IMAGE_SOURCE"
    ensure_image "$VEP_SIF" "$VEP_IMAGE_SOURCE"

    PICARD_JAR="$(detect_core_path \
        "${PICARD_JAR:-}" \
        /usr/local/bin/picard.jar \
        /data/picard.jar \
    )" || die "could not find Picard in $CORE_SIF"
    GATK_JAR="$(detect_core_path \
        "${GATK_JAR:-}" \
        /usr/local/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar \
        /data/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar \
    )" || die "could not find GATK in $CORE_SIF"

    local genome_gz="${BUILD_DIR}/sources/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    local gtf_gz="${BUILD_DIR}/sources/Homo_sapiens.GRCh38.107.chr.gtf.gz"
    local cdna_gz="${BUILD_DIR}/sources/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    local known_vcf_gz="${BUILD_DIR}/sources/1000GENOMES-phase_3.vcf.gz"
    local vep_tar="${BUILD_DIR}/sources/homo_sapiens_vep_113_GRCh38.tar.gz"

    download_file "$GENOME_URL" "$genome_gz"
    download_file "$GTF_URL" "$gtf_gz"
    download_file "$CDNA_URL" "$cdna_gz"
    download_file "$KNOWN_VCF_URL" "$known_vcf_gz"
    download_file "$VEP_CACHE_URL" "$vep_tar"

    gzip -t "$genome_gz"
    gzip -t "$gtf_gz"
    gzip -t "$cdna_gz"
    gzip -t "$known_vcf_gz"
    tar -tzf "$vep_tar" >/dev/null

    log "Unpacking reference inputs"
    gzip -dc "$genome_gz" >"${BUILD_DIR}/ref.fa"
    gzip -dc "$gtf_gz" >"${BUILD_DIR}/ref.gtf"
    gzip -dc "$cdna_gz" >"${BUILD_DIR}/transcriptom.fa"
    gzip -dc "$known_vcf_gz" >"${BUILD_DIR}/knowns_variants.vcf"
    tar -xzf "$vep_tar" -C "${BUILD_DIR}/VEP" --no-same-owner
    [[ -d "${BUILD_DIR}/VEP/homo_sapiens/113_GRCh38" ]] \
        || die "VEP cache did not contain homo_sapiens/113_GRCh38"

    cp -p -- "${BUILD_DIR}/transcriptom.fa" "${BUILD_DIR}/cdna.fa"

    log "Building samtools, Picard, GATK, and Kallisto indexes"
    core_exec samtools faidx /work/ref.fa
    core_exec java -jar "$PICARD_JAR" CreateSequenceDictionary \
        R=/work/ref.fa O=/work/ref.dict
    core_exec java -jar "$GATK_JAR" IndexFeatureFile \
        -I /work/knowns_variants.vcf
    kallisto_exec kallisto index \
        -i /work/kalliso_index /work/transcriptom.fa

    if ! ln -- "${BUILD_DIR}/kalliso_index" "${BUILD_DIR}/kallisto_index" 2>/dev/null; then
        cp -p -- "${BUILD_DIR}/kalliso_index" "${BUILD_DIR}/kallisto_index"
    fi

    log "Generating pipeline metadata artifacts"
    awk -F $'\t' '$3 == "gene" { print }' "${BUILD_DIR}/ref.gtf" \
        >"${BUILD_DIR}/ref.GeneLvlOnly.gtf"
    awk -F $'\t' '$3 == "exon" { print $1 "\t" ($4 - 1) "\t" $5 }' \
        "${BUILD_DIR}/ref.gtf" \
        | sort -k1,1 -k2,2n -k3,3n -u \
        >"${BUILD_DIR}/ref.ExonsOnly.bed"

    awk -F $'\t' '$3 == "exon" { print $4 "\t" $5 "\t" $7 "\t" $9 }' "${BUILD_DIR}/ref.gtf" \
        | awk '{for(i=5;i<=NF;i++){if($i~/^"ENSE/){a=$i}} print a, $1,$2,$3,$5,$15}' \
        | sed 's/"//g;s/;//g' \
        | sort -d \
        | awk 'BEGIN {print "exon_id\tstart\tend\tstrand\tgene_id\tgene_name"} { print }' \
        >"${BUILD_DIR}/Exon_gtf_info.tab"

    r_exec Rscript /pipeline/procGTF.R \
        /work/ref.GeneLvlOnly.gtf \
        /work/refGeneID_ensembl_v107

    log "Building STAR index"
    core_exec STAR \
        --runThreadN "${SLURM_CPUS_PER_TASK:-16}" \
        --runMode genomeGenerate \
        --genomeDir /work \
        --genomeFastaFiles /work/ref.fa \
        --sjdbOverhang 100 \
        --sjdbGTFfile /work/ref.gtf \
        --genomeSAindexNbases 11

    record_tool_versions
}

validate_reference() {
    local required
    for required in \
        ref.fa ref.fa.fai ref.dict ref.gtf transcriptom.fa cdna.fa \
        kalliso_index kallisto_index knowns_variants.vcf knowns_variants.vcf.idx \
        ref.ExonsOnly.bed ref.GeneLvlOnly.gtf Exon_gtf_info.tab \
        refGeneID_ensembl_v107.rds refGeneID_ensembl_v107.tsv \
        Genome SA SAindex; do
        [[ -s "${BUILD_DIR}/${required}" ]] || die "missing reference artifact: ${required}"
    done

    [[ -d "${BUILD_DIR}/VEP/homo_sapiens/113_GRCh38" ]] \
        || die "missing VEP cache directory"

    log "Validating Kallisto index"
    kallisto_exec kallisto inspect -i /work/kalliso_index >/dev/null

    log "Validating exon BED with GATK"
    core_exec java -jar "$GATK_JAR" BedToIntervalList \
        -I /work/ref.ExonsOnly.bed \
        -O /work/ref.ExonsOnly.interval_list \
        -R /work/ref.fa \
        --SEQUENCE_DICTIONARY /work/ref.dict

    log "Validating VEP cache offline"
    {
        printf '##fileformat=VCFv4.2\n'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        awk '!/^#/ { print; exit }' "${BUILD_DIR}/knowns_variants.vcf"
    } >"${BUILD_DIR}/vep_smoke.vcf"
    [[ -s "${BUILD_DIR}/vep_smoke.vcf" ]] || die "could not create VEP smoke VCF"

    vep_exec vep \
        -i /work/vep_smoke.vcf \
        -o /work/vep_smoke.out.vcf \
        --format vcf \
        --vcf \
        --offline \
        --cache \
        --dir_cache /work/VEP \
        --fasta /work/ref.fa \
        --species homo_sapiens \
        --assembly GRCh38 \
        --force_overwrite
    [[ -s "${BUILD_DIR}/vep_smoke.out.vcf" ]] || die "VEP smoke annotation produced no output"

    rm -f -- "${BUILD_DIR}/ref.ExonsOnly.interval_list" \
        "${BUILD_DIR}/vep_smoke.vcf" "${BUILD_DIR}/vep_smoke.out.vcf"
}

stage_reference() {
    rm -rf -- "$STAGE_DIR"
    mkdir -p -- "$STAGE_DIR"

    local path name
    while IFS= read -r -d '' path; do
        name="$(basename "$path")"
        case "$name" in
            *.gz|tool-versions.*|vep_smoke.*|ref.ExonsOnly.interval_list)
                continue
                ;;
        esac
        cp -a -- "$path" "$STAGE_DIR/"
    done < <(find "$BUILD_DIR" -maxdepth 1 -type f -print0)

    cp -a -- "${BUILD_DIR}/VEP" "$STAGE_DIR/VEP"
}

write_manifest() {
    export MANIFEST_BUILD_DIR="$BUILD_DIR"
    export MANIFEST_STAGE_DIR="$STAGE_DIR"
    export MANIFEST_REFERENCE_ID="$REFERENCE_ID"
    export MANIFEST_CORE_IMAGE_SOURCE="$CORE_IMAGE_SOURCE"
    export MANIFEST_KALLISTO_IMAGE_SOURCE="$KALLISTO_IMAGE_SOURCE"
    export MANIFEST_R_IMAGE_SOURCE="$R_IMAGE_SOURCE"
    export MANIFEST_VEP_IMAGE_SOURCE="$VEP_IMAGE_SOURCE"
    export MANIFEST_CORE_SIF="$CORE_SIF"
    export MANIFEST_KALLISTO_SIF="$KALLISTO_SIF"
    export MANIFEST_R_SIF="$R_SIF"
    export MANIFEST_VEP_SIF="$VEP_SIF"
    export MANIFEST_GENOME_URL="$GENOME_URL"
    export MANIFEST_GTF_URL="$GTF_URL"
    export MANIFEST_CDNA_URL="$CDNA_URL"
    export MANIFEST_KNOWN_VCF_URL="$KNOWN_VCF_URL"
    export MANIFEST_VEP_CACHE_URL="$VEP_CACHE_URL"
    export MANIFEST_PICARD_JAR="$PICARD_JAR"
    export MANIFEST_GATK_JAR="$GATK_JAR"

    python3 - <<'PY'
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

build = Path(os.environ["MANIFEST_BUILD_DIR"])
stage = Path(os.environ["MANIFEST_STAGE_DIR"])

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def text(path: Path) -> str:
    return path.read_text(errors="replace").strip()

source_paths = {
    "genome_fasta": build / "sources/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "annotation_gtf": build / "sources/Homo_sapiens.GRCh38.107.chr.gtf.gz",
    "transcriptome_cdna": build / "sources/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    "known_variants": build / "sources/1000GENOMES-phase_3.vcf.gz",
    "vep_cache_archive": build / "sources/homo_sapiens_vep_113_GRCh38.tar.gz",
}
source_urls = {
    "genome_fasta": os.environ["MANIFEST_GENOME_URL"],
    "annotation_gtf": os.environ["MANIFEST_GTF_URL"],
    "transcriptome_cdna": os.environ["MANIFEST_CDNA_URL"],
    "known_variants": os.environ["MANIFEST_KNOWN_VCF_URL"],
    "vep_cache_archive": os.environ["MANIFEST_VEP_CACHE_URL"],
}

artifacts = []
for path in sorted(stage.rglob("*")):
    if path.is_file() and path.name != "reference_manifest.json":
        artifacts.append({
            "path": str(path.relative_to(stage)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

containers = {}
for name, source_key, sif_key in [
    ("core", "MANIFEST_CORE_IMAGE_SOURCE", "MANIFEST_CORE_SIF"),
    ("kallisto", "MANIFEST_KALLISTO_IMAGE_SOURCE", "MANIFEST_KALLISTO_SIF"),
    ("r", "MANIFEST_R_IMAGE_SOURCE", "MANIFEST_R_SIF"),
    ("vep", "MANIFEST_VEP_IMAGE_SOURCE", "MANIFEST_VEP_SIF"),
]:
    sif = Path(os.environ[sif_key])
    containers[name] = {
        "source": os.environ[source_key],
        "sif_sha256": sha256(sif),
    }

versions = {}
for name in ("core", "samtools", "java", "kallisto", "r", "vep"):
    versions[name] = text(build / f"tool-versions.{name}.txt")

manifest = {
    "manifest_schema_version": 1,
    "reference_id": os.environ["MANIFEST_REFERENCE_ID"],
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "pipeline_compatibility": {
        "pipeline": "RNApipeline",
        "pipeline_version": "1.7.0",
        "star_sjdb_overhang": 100,
        "star_genome_sa_index_nbases": 11,
        "kallisto_index_names": ["kalliso_index", "kallisto_index"],
    },
    "genome": {
        "species": "homo_sapiens",
        "assembly": "GRCh38",
        "ensembl_release": 107,
    },
    "annotation": {"ensembl_release": 107},
    "transcriptome": {"ensembl_release": 107},
    "known_variants": {
        "ensembl_release": 105,
        "file": "knowns_variants.vcf",
    },
    "vep": {
        "cache_release": 113,
        "cache_path": "VEP",
        "species": "homo_sapiens",
        "assembly": "GRCh38",
    },
    "source_inputs": {
        name: {"url": source_urls[name], "sha256": sha256(path)}
        for name, path in source_paths.items()
    },
    "containers": containers,
    "tool_versions": versions,
    "tool_paths": {
        "picard_jar": os.environ["MANIFEST_PICARD_JAR"],
        "gatk_jar": os.environ["MANIFEST_GATK_JAR"],
    },
    "artifacts": artifacts,
}

(stage / "reference_manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
PY
}

validate_manifest() {
    python3 - "$STAGE_DIR/reference_manifest.json" "$STAGE_DIR" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
stage = Path(sys.argv[2])
manifest = json.loads(manifest_path.read_text())

assert manifest["reference_id"] == "ensembl_v107_GRCh38"
assert manifest["genome"]["ensembl_release"] == 107
assert manifest["annotation"]["ensembl_release"] == 107
assert manifest["transcriptome"]["ensembl_release"] == 107
assert manifest["known_variants"]["ensembl_release"] == 105
assert manifest["vep"]["cache_release"] == 113

for artifact in manifest["artifacts"]:
    path = stage / artifact["path"]
    assert path.is_file(), path
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    assert digest.hexdigest() == artifact["sha256"], path
PY
}

main() {
    require_cmd awk
    require_cmd curl
    require_cmd find
    require_cmd gzip
    require_cmd python3
    require_cmd sha256sum
    require_cmd sort
    require_cmd tar

    APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"
    require_cmd "$APPTAINER_BIN"
    APPTAINER_ARCH="${APPTAINER_PULL_ARCH:-$(host_apptainer_arch)}"
    case "$APPTAINER_ARCH" in
        amd64|arm64)
            ;;
        *)
            die "APPTAINER_PULL_ARCH must be amd64 or arm64 (got: $APPTAINER_ARCH)"
            ;;
    esac

    check_ref_mount
    check_source_urls

    if ((CHECK_ONLY == 1)); then
        if [[ -e "$FINAL_DIR" ]]; then
            log "Reference already exists: $FINAL_DIR"
        else
            log "Reference destination is available: $FINAL_DIR"
        fi
        log "Preflight checks passed"
        exit 0
    fi

    if [[ -e "$FINAL_DIR" || -L "$FINAL_DIR" ]] && ((FORCE == 0)); then
        die "reference already exists: $FINAL_DIR (use --force to preserve and replace it)"
    fi

    mkdir -p -- "$WORK_ROOT"
    [[ -w "$WORK_ROOT" ]] || die "work root is not writable: $WORK_ROOT"
    acquire_build_lock

    build_reference
    validate_reference
    stage_reference
    write_manifest
    validate_manifest

    local backup=""
    if [[ -e "$FINAL_DIR" || -L "$FINAL_DIR" ]]; then
        backup="${REF_ROOT%/}/.staging/${REFERENCE_ID}.previous.${CURRENT_JOB}"
        log "Preserving existing reference at: $backup"
        mv -- "$FINAL_DIR" "$backup"
    fi

    if ! mv -- "$STAGE_DIR" "$FINAL_DIR"; then
        if [[ -n "$backup" && -e "$backup" ]]; then
            mv -- "$backup" "$FINAL_DIR" || true
        fi
        die "could not publish reference to: $FINAL_DIR"
    fi

    [[ -s "${FINAL_DIR}/reference_manifest.json" ]] || die "published reference has no manifest"
    log "Published reference: $FINAL_DIR"
}

main "$@"
