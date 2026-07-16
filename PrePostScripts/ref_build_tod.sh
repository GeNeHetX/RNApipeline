#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# Build the complete precomputed reference used by RNApipeline on TOD/PAM.
# The build is staged on local Slurm scratch and published atomically to the
# direct CephFS reference mount.

REFERENCE_ID="ensembl_v107_GRCh38"
PIPELINE_VERSION="1.7.0"
KALLISTO_VERSION="0.51.1"
REF_ROOT="/ref"
WORK_ROOT="${SLURM_TMPDIR:-/srv/slurm/scratch}"
CHECK_ONLY=0
FORCE=0
KEEP_WORK=0
SKIP_MOUNT_CHECK=0
KALLISTO_ONLY=0

CORE_IMAGE_SOURCE="docker://genehetx/genehetx-rnaseq:v1.6.1"
KALLISTO_IMAGE_SOURCE="docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0"
KALLISTO_SOURCE_URL="https://github.com/pachterlab/kallisto/archive/refs/tags/v0.51.1.tar.gz"
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
  --kallisto-only          Build and publish the ARM64 and AMD64 Kallisto SIFs,
                           then exit without rebuilding the reference.
  --force                  Replace an existing reference after preserving it
                           under /ref/.staging/. Also replace the Kallisto SIFs.
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
        --kallisto-only)
            KALLISTO_ONLY=1
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
KALLISTO_TOOL_ROOT="${REF_ROOT%/}/tools/rnapipeline/v${PIPELINE_VERSION}/kallisto/${KALLISTO_VERSION}"
KALLISTO_ARM_SIF_PATH="${KALLISTO_TOOL_ROOT}/arm64/kallisto-${KALLISTO_VERSION}.sif"
KALLISTO_AMD64_SIF_PATH="${KALLISTO_TOOL_ROOT}/amd64/kallisto-${KALLISTO_VERSION}.sif"
KALLISTO_SOURCE_ARCHIVE="${BUILD_DIR}/sources/kallisto-v${KALLISTO_VERSION}.tar.gz"
KALLISTO_STAGE_DIR=""

CORE_SIF="${IMAGE_DIR}/genehetx-rnaseq-v1.6.1.sif"
KALLISTO_SIF="${IMAGE_DIR}/kallisto-${KALLISTO_VERSION}-arm64.sif"
KALLISTO_AMD64_SIF="${IMAGE_DIR}/kallisto-${KALLISTO_VERSION}-amd64.sif"
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
        if ((rc != 0)) && [[ -n "${KALLISTO_STAGE_DIR:-}" && -d "${KALLISTO_STAGE_DIR:-}" ]]; then
            rm -rf -- "$KALLISTO_STAGE_DIR"
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
    for url in "$GENOME_URL" "$GTF_URL" "$CDNA_URL" "$KNOWN_VCF_URL" "$VEP_CACHE_URL" "$KALLISTO_SOURCE_URL"; do
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

pull_arch_image() {
    local target="$1"
    local source="$2"
    local arch="$3"

    if [[ -s "$target" ]]; then
        log "Using cached ${arch} image: $target"
        return
    fi

    log "Pulling ${arch} Apptainer image: $source"
    mkdir -p -- "$(dirname "$target")"
    "$APPTAINER_BIN" pull --arch "$arch" "$target" "$source"
    [[ -s "$target" ]] || die "Apptainer pull did not produce: $target"
    "$APPTAINER_BIN" inspect "$target" >/dev/null \
        || die "Apptainer produced an invalid SIF: $target"
}

write_kallisto_definition() {
    local definition="$1"
    local build_cpus="${SLURM_CPUS_PER_TASK:-16}"

    cat >"$definition" <<EOF
Bootstrap: docker
From: ubuntu:22.04

%files
    ${KALLISTO_SOURCE_ARCHIVE} /tmp/kallisto-source.tar.gz

%post
    set -eux
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \\
        autoconf automake build-essential ca-certificates cmake \\
        libhdf5-dev procps zlib1g-dev
    tar -xzf /tmp/kallisto-source.tar.gz -C /opt
    cd /opt/kallisto-${KALLISTO_VERSION}/ext/htslib
    autoheader
    autoconf
    cd ../..
    cmake -S . -B build \\
        -DCMAKE_BUILD_TYPE=Release \\
        -DCMAKE_INSTALL_PREFIX=/usr/local \\
        -DUSE_HDF5=ON
    cmake --build build --parallel ${build_cpus}
    install -m 0755 build/src/kallisto /usr/local/bin/kallisto
    kallisto version | grep -F 'kallisto, version ${KALLISTO_VERSION}'
    rm -rf /opt/kallisto-${KALLISTO_VERSION} /tmp/kallisto-source.tar.gz
    apt-get clean
    rm -rf /var/lib/apt/lists/*

%environment
    export PATH="/usr/local/bin:\$PATH"

%test
    /usr/local/bin/kallisto version
EOF
}

build_arm_kallisto_sif() {
    local target="$1"
    local definition="${BUILD_DIR}/kallisto-${KALLISTO_VERSION}-arm64.def"

    [[ "$(host_apptainer_arch)" == "arm64" ]] \
        || die "the ARM64 Kallisto SIF must be built on an ARM64 Slurm worker"

    write_kallisto_definition "$definition"
    log "Building native ARM64 Kallisto ${KALLISTO_VERSION} SIF"
    "$APPTAINER_BIN" build --fakeroot "$target" "$definition"
    [[ -s "$target" ]] || die "ARM64 Kallisto build did not produce: $target"
    image_matches_arch "$target" \
        || die "built Kallisto SIF is not executable on ARM64: $target"
}

prepare_kallisto_sifs() {
    local arm_work="$KALLISTO_SIF"
    local amd64_work="$KALLISTO_AMD64_SIF"

    mkdir -p -- "${BUILD_DIR}/sources" "$IMAGE_DIR"
    download_file "$KALLISTO_SOURCE_URL" "$KALLISTO_SOURCE_ARCHIVE"
    tar -tzf "$KALLISTO_SOURCE_ARCHIVE" \
        | awk -v expected="kallisto-${KALLISTO_VERSION}/" '$0 == expected { found = 1 } END { exit !found }' \
        || die "Kallisto source archive does not contain kallisto-${KALLISTO_VERSION}/"

    if [[ -e "$KALLISTO_TOOL_ROOT" || -L "$KALLISTO_TOOL_ROOT" ]] && ((FORCE == 0)); then
        [[ -s "$KALLISTO_ARM_SIF_PATH" && -s "$KALLISTO_AMD64_SIF_PATH" \
            && -x "${KALLISTO_TOOL_ROOT}/bin/kallisto" \
            && -s "${KALLISTO_TOOL_ROOT}/kallisto_manifest.json" ]] \
            || die "existing Kallisto SIF bundle is incomplete: $KALLISTO_TOOL_ROOT (use --force)"
        if ! image_matches_arch "$KALLISTO_ARM_SIF_PATH"; then
            die "existing ARM64 Kallisto SIF is unusable: $KALLISTO_ARM_SIF_PATH (use --force)"
        fi
        cp -p -- "$KALLISTO_ARM_SIF_PATH" "$arm_work"
        cp -p -- "$KALLISTO_AMD64_SIF_PATH" "$amd64_work"
    else
        build_arm_kallisto_sif "$arm_work"
        pull_arch_image "$amd64_work" "$KALLISTO_IMAGE_SOURCE" amd64
    fi
}

validate_kallisto_sifs() {
    [[ -s "$KALLISTO_SIF" ]] || die "missing ARM64 Kallisto SIF: $KALLISTO_SIF"
    [[ -s "$KALLISTO_AMD64_SIF" ]] || die "missing AMD64 Kallisto SIF: $KALLISTO_AMD64_SIF"
    image_matches_arch "$KALLISTO_SIF" \
        || die "ARM64 Kallisto SIF is not executable on this worker: $KALLISTO_SIF"
    "$APPTAINER_BIN" inspect "$KALLISTO_AMD64_SIF" >/dev/null \
        || die "AMD64 Kallisto SIF is invalid: $KALLISTO_AMD64_SIF"

    local version
    version="$(kallisto_exec kallisto version 2>&1)"
    [[ "$version" == *"${KALLISTO_VERSION}"* ]] \
        || die "ARM64 Kallisto version mismatch: $version"
    log "Validated native ARM64 Kallisto ${KALLISTO_VERSION} and AMD64 SIF"
}

write_kallisto_manifest() {
    local stage="$1"
    export KALLISTO_MANIFEST_STAGE="$stage"
    export KALLISTO_MANIFEST_ARM_SIF="$KALLISTO_SIF"
    export KALLISTO_MANIFEST_AMD64_SIF="$KALLISTO_AMD64_SIF"
    export KALLISTO_MANIFEST_SOURCE_ARCHIVE="$KALLISTO_SOURCE_ARCHIVE"
    export KALLISTO_MANIFEST_SOURCE_URL="$KALLISTO_SOURCE_URL"
    export KALLISTO_MANIFEST_IMAGE_SOURCE="$KALLISTO_IMAGE_SOURCE"
    export KALLISTO_MANIFEST_TOOL_ROOT="$KALLISTO_TOOL_ROOT"

    python3 - <<'PY'
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

stage = Path(os.environ["KALLISTO_MANIFEST_STAGE"])

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

manifest = {
    "manifest_schema_version": 1,
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "pipeline": "RNApipeline",
    "pipeline_version": "1.7.0",
    "tool": "kallisto",
    "tool_version": "0.51.1",
    "source": {
        "arm64_build_url": os.environ["KALLISTO_MANIFEST_SOURCE_URL"],
        "arm64_source_archive_sha256": sha256(Path(os.environ["KALLISTO_MANIFEST_SOURCE_ARCHIVE"])),
        "amd64_image": os.environ["KALLISTO_MANIFEST_IMAGE_SOURCE"],
    },
    "containers": {
        "arm64": {
            "path": os.path.join(os.environ["KALLISTO_MANIFEST_TOOL_ROOT"], "arm64", "kallisto-0.51.1.sif"),
            "sha256": sha256(Path(os.environ["KALLISTO_MANIFEST_ARM_SIF"])),
            "build": "pachterlab/kallisto v0.51.1 built natively on ARM64",
        },
        "amd64": {
            "path": os.path.join(os.environ["KALLISTO_MANIFEST_TOOL_ROOT"], "amd64", "kallisto-0.51.1.sif"),
            "sha256": sha256(Path(os.environ["KALLISTO_MANIFEST_AMD64_SIF"])),
            "build": "quay.io/biocontainers pinned Kallisto 0.51.1 image",
        },
    },
    "wrapper": {
        "path": os.path.join(os.environ["KALLISTO_MANIFEST_TOOL_ROOT"], "bin", "kallisto"),
        "sha256": sha256(stage / "bin/kallisto"),
        "selects_by_host_architecture": True,
    },
}

(stage / "kallisto_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
PY
}

write_kallisto_wrapper() {
    local wrapper="$1"

    cat >"$wrapper" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "$(uname -m)" in
    aarch64|arm64)
        IMAGE="${ROOT_DIR}/arm64/kallisto-0.51.1.sif"
        ;;
    x86_64)
        IMAGE="${ROOT_DIR}/amd64/kallisto-0.51.1.sif"
        ;;
    *)
        echo "Unsupported host architecture for Kallisto: $(uname -m)" >&2
        exit 1
        ;;
esac

APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"
command -v "$APPTAINER_BIN" >/dev/null 2>&1 \
    || { echo "Apptainer is required to run Kallisto" >&2; exit 1; }
[[ -s "$IMAGE" ]] \
    || { echo "No Kallisto SIF for host architecture: $IMAGE" >&2; exit 1; }

exec "$APPTAINER_BIN" exec --no-home "$IMAGE" kallisto "$@"
EOF
    chmod 0755 "$wrapper"
}

publish_kallisto_sifs() {
    local final_parent="${REF_ROOT%/}/tools/rnapipeline/v${PIPELINE_VERSION}/kallisto"
    local final_dir="$KALLISTO_TOOL_ROOT"
    local backup=""

    if [[ -e "$final_dir" || -L "$final_dir" ]] && ((FORCE == 0)); then
        log "Kallisto SIF bundle already published: $final_dir"
        return
    fi

    KALLISTO_STAGE_DIR="${REF_ROOT%/}/.staging/rnapipeline-kallisto.${CURRENT_JOB}"
    rm -rf -- "$KALLISTO_STAGE_DIR"
    mkdir -p -- "${KALLISTO_STAGE_DIR}/arm64" "${KALLISTO_STAGE_DIR}/amd64"
    cp -p -- "$KALLISTO_SIF" \
        "${KALLISTO_STAGE_DIR}/arm64/kallisto-${KALLISTO_VERSION}.sif"
    cp -p -- "$KALLISTO_AMD64_SIF" \
        "${KALLISTO_STAGE_DIR}/amd64/kallisto-${KALLISTO_VERSION}.sif"
    mkdir -p -- "${KALLISTO_STAGE_DIR}/bin"
    write_kallisto_wrapper "${KALLISTO_STAGE_DIR}/bin/kallisto"
    write_kallisto_manifest "$KALLISTO_STAGE_DIR"

    mkdir -p -- "$final_parent"
    if [[ -e "$final_dir" || -L "$final_dir" ]]; then
        backup="${REF_ROOT%/}/.staging/rnapipeline-kallisto.previous.${CURRENT_JOB}"
        log "Preserving existing Kallisto SIF bundle at: $backup"
        mv -- "$final_dir" "$backup"
    fi

    if ! mv -- "$KALLISTO_STAGE_DIR" "$final_dir"; then
        if [[ -n "$backup" && -e "$backup" ]]; then
            mv -- "$backup" "$final_dir" || true
        fi
        die "could not publish Kallisto SIF bundle to: $final_dir"
    fi
    KALLISTO_STAGE_DIR=""
    log "Published ARM64 and AMD64 Kallisto SIFs: $final_dir"
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
        --pwd /work \
        "$R_SIF" "$@"
}

write_proc_gtf_script() {
    local target="${BUILD_DIR}/procGTF.R"
    [[ -s "$target" ]] && return
    cat >"$target" <<'R_SCRIPT'
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(parallel)
gtf=read.delim(args[1],sep="\t",as.is=T,header=F,comment.char="#")
gtf.getmeta=function(apieceofgtf){strsplit((apieceofgtf)[,9],";| ")}
gtf.getmetavalue=function(apieceofgtf,field){
   metadata=gtf.getmeta(apieceofgtf)
   unlist(mclapply(metadata,function(x){
       if(field%in% x){return(x[which(x==field)+1])}else{return(NA)}
   }))
}
geneTab = gtf[which(gtf[,3] == "gene"),]
geneTab$GeneID = gtf.getmetavalue(geneTab,"gene_id")
geneTab$GeneName = gtf.getmetavalue(geneTab,"gene_name")
inoname = which(is.na(geneTab$GeneName))
geneTab$GeneName[inoname] = geneTab$GeneID[inoname]
geneTab$biotype = gtf.getmetavalue(geneTab,"gene_biotype")
geneTab = unique(geneTab)
rownames(geneTab) = geneTab$GeneID
geneTab = unique(geneTab[,- c(9,6,2,3,8)])
colnames(geneTab)[1:4] = c("seqname","start","end", "strand")
saveRDS(geneTab,file=paste0(args[2],".rds"))
write.table(geneTab,file=paste0(args[2],".tsv"),quote=F,sep="\t")
R_SCRIPT
    chmod 0755 "$target"
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
    prepare_kallisto_sifs
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
    write_proc_gtf_script
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

    r_exec Rscript /work/procGTF.R \
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
    export MANIFEST_KALLISTO_SOURCE_URL="$KALLISTO_SOURCE_URL"
    export MANIFEST_KALLISTO_SOURCE_ARCHIVE="$KALLISTO_SOURCE_ARCHIVE"
    export MANIFEST_R_IMAGE_SOURCE="$R_IMAGE_SOURCE"
    export MANIFEST_VEP_IMAGE_SOURCE="$VEP_IMAGE_SOURCE"
    export MANIFEST_CORE_SIF="$CORE_SIF"
    export MANIFEST_KALLISTO_SIF="$KALLISTO_SIF"
    export MANIFEST_KALLISTO_AMD64_SIF="$KALLISTO_AMD64_SIF"
    export MANIFEST_KALLISTO_TOOL_ROOT="$KALLISTO_TOOL_ROOT"
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
    "kallisto_source": build / "sources/kallisto-v0.51.1.tar.gz",
}
source_urls = {
    "genome_fasta": os.environ["MANIFEST_GENOME_URL"],
    "annotation_gtf": os.environ["MANIFEST_GTF_URL"],
    "transcriptome_cdna": os.environ["MANIFEST_CDNA_URL"],
    "known_variants": os.environ["MANIFEST_KNOWN_VCF_URL"],
    "vep_cache_archive": os.environ["MANIFEST_VEP_CACHE_URL"],
    "kallisto_source": os.environ["MANIFEST_KALLISTO_SOURCE_URL"],
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

containers["kallisto"].update({
    "architecture": os.environ["APPTAINER_ARCH"],
    "source_build_url": os.environ["MANIFEST_KALLISTO_SOURCE_URL"],
    "source_build_archive_sha256": sha256(Path(os.environ["MANIFEST_KALLISTO_SOURCE_ARCHIVE"])),
    "architectures": {
        "arm64": {
            "path": os.path.join(os.environ["MANIFEST_KALLISTO_TOOL_ROOT"], "arm64", "kallisto-0.51.1.sif"),
            "sha256": sha256(Path(os.environ["MANIFEST_KALLISTO_SIF"])),
        },
        "amd64": {
            "path": os.path.join(os.environ["MANIFEST_KALLISTO_TOOL_ROOT"], "amd64", "kallisto-0.51.1.sif"),
            "sha256": sha256(Path(os.environ["MANIFEST_KALLISTO_AMD64_SIF"])),
        },
    },
})

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
    if ((KALLISTO_ONLY == 1)); then
        log "Checking source URL: $KALLISTO_SOURCE_URL"
        curl -fsSIL --retry 3 --connect-timeout 10 --max-time 60 "$KALLISTO_SOURCE_URL" >/dev/null \
            || die "source URL is not reachable: $KALLISTO_SOURCE_URL"
    else
        check_source_urls
    fi

    if ((CHECK_ONLY == 1)); then
        if ((KALLISTO_ONLY == 1)); then
            if [[ -e "$KALLISTO_TOOL_ROOT" ]]; then
                log "Kallisto SIF bundle already exists: $KALLISTO_TOOL_ROOT"
            else
                log "Kallisto SIF bundle will be published to: $KALLISTO_TOOL_ROOT"
            fi
        else
            if [[ -e "$FINAL_DIR" ]]; then
                log "Reference already exists: $FINAL_DIR"
            else
                log "Reference destination is available: $FINAL_DIR"
            fi
        fi
        log "Preflight checks passed"
        exit 0
    fi

    mkdir -p -- "$WORK_ROOT"
    [[ -w "$WORK_ROOT" ]] || die "work root is not writable: $WORK_ROOT"
    acquire_build_lock
    export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/srv/slurm/scratch/.apptainer-cache}"
    export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${BUILD_DIR}/apptainer-tmp}"
    mkdir -p -- "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

    if ((KALLISTO_ONLY == 1)); then
        prepare_kallisto_sifs
        validate_kallisto_sifs
        publish_kallisto_sifs
        log "Kallisto-only build completed"
        exit 0
    fi

    if [[ -e "$FINAL_DIR" || -L "$FINAL_DIR" ]] && ((FORCE == 0)); then
        die "reference already exists: $FINAL_DIR (use --force to preserve and replace it)"
    fi

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
    publish_kallisto_sifs
}

main "$@"
