#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

VERSION="0.51.1"
SOURCE_URL="https://github.com/pachterlab/kallisto/archive/refs/tags/v${VERSION}.tar.gz"
DEST_DIR="${1:-}"

[[ -n "$DEST_DIR" && "$#" -eq 1 ]] || {
    echo "Usage: $0 DESTINATION_DIRECTORY" >&2
    exit 2
}
[[ "$(uname -m)" == "aarch64" || "$(uname -m)" == "arm64" ]] || {
    echo "This script must run on an ARM64 host" >&2
    exit 1
}

command -v apptainer >/dev/null 2>&1 || { echo "apptainer is required" >&2; exit 1; }
command -v curl >/dev/null 2>&1 || { echo "curl is required" >&2; exit 1; }

CPUS="${KALLISTO_BUILD_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
if [[ -z "${KALLISTO_BUILD_CPUS:-}" && -z "${SLURM_CPUS_PER_TASK:-}" ]] \
    && command -v nproc >/dev/null 2>&1; then
    CPUS="$(nproc)"
fi

WORK_DIR="$(mktemp -d "${TMPDIR:-/srv/slurm/scratch}/kallisto-arm64.XXXXXX")"
OUTPUT="${DEST_DIR%/}/kallisto-${VERSION}.sif"
TEMP_SIF="${DEST_DIR%/}/.kallisto-${VERSION}.${BASHPID}.sif"
trap 'rm -rf "$WORK_DIR" "$TEMP_SIF"' EXIT

mkdir -p "$DEST_DIR"
curl --fail --location --retry 5 --retry-delay 5 \
    --output "${WORK_DIR}/kallisto-source.tar.gz" "$SOURCE_URL"

cat >"${WORK_DIR}/kallisto.def" <<EOF
Bootstrap: docker
From: ubuntu:22.04

%files
    ${WORK_DIR}/kallisto-source.tar.gz /opt/kallisto-source.tar.gz

%post
    set -eux
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \\
        autoconf automake build-essential ca-certificates cmake \\
        libhdf5-dev procps zlib1g-dev
    tar -xzf /opt/kallisto-source.tar.gz -C /opt
    cd /opt/kallisto-${VERSION}/ext/htslib
    autoheader
    autoconf
    cd ../..
    cmake -S . -B build \\
        -DCMAKE_BUILD_TYPE=Release \\
        -DCMAKE_INSTALL_PREFIX=/usr/local \\
        -DUSE_HDF5=ON
    # Kallisto's CMake graph can link before its external Bifrost target is
    # complete when both are built in one parallel invocation.
    cmake --build build --target bifrost --parallel ${CPUS}
    cmake --build build --parallel ${CPUS}
    install -m 0755 build/src/kallisto /usr/local/bin/kallisto
    kallisto version | grep -F 'kallisto, version ${VERSION}'
    rm -rf /opt/kallisto-${VERSION} /opt/kallisto-source.tar.gz
    apt-get clean
    rm -rf /var/lib/apt/lists/*

%environment
    export PATH="/usr/local/bin:\$PATH"

%test
    /usr/local/bin/kallisto version
EOF

apptainer build --fakeroot "$TEMP_SIF" "${WORK_DIR}/kallisto.def"
[[ "$(apptainer exec "$TEMP_SIF" uname -m)" == "aarch64" || \
   "$(apptainer exec "$TEMP_SIF" uname -m)" == "arm64" ]] || {
    echo "Built SIF is not ARM64" >&2
    exit 1
}
apptainer exec "$TEMP_SIF" kallisto version | grep -F "kallisto, version ${VERSION}" >/dev/null
mv -f "$TEMP_SIF" "$OUTPUT"
echo "Wrote $OUTPUT"
