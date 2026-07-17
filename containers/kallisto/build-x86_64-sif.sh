#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

VERSION="0.51.1"
SOURCE="docker://quay.io/biocontainers/kallisto:0.51.1--heb0cbe2_0"
DEST_DIR="${1:-}"

[[ -n "$DEST_DIR" && "$#" -eq 1 ]] || {
    echo "Usage: $0 DESTINATION_DIRECTORY" >&2
    exit 2
}
[[ "$(uname -m)" == "x86_64" ]] || {
    echo "This script must run on an x86_64 host" >&2
    exit 1
}

command -v apptainer >/dev/null 2>&1 || { echo "apptainer is required" >&2; exit 1; }

WORK_DIR="$(mktemp -d "${TMPDIR:-/srv/slurm/scratch}/kallisto-x86_64.XXXXXX")"
OUTPUT="${DEST_DIR%/}/kallisto-${VERSION}.sif"
TEMP_SIF="${DEST_DIR%/}/.kallisto-${VERSION}.${BASHPID}.sif"
trap 'rm -rf "$WORK_DIR" "$TEMP_SIF"' EXIT

mkdir -p "$DEST_DIR"
apptainer pull --arch amd64 "$TEMP_SIF" "$SOURCE"
[[ "$(apptainer exec "$TEMP_SIF" uname -m)" == "x86_64" ]] || {
    echo "Pulled SIF is not x86_64" >&2
    exit 1
}
apptainer exec "$TEMP_SIF" kallisto version | grep -F "kallisto, version ${VERSION}" >/dev/null
mv -f "$TEMP_SIF" "$OUTPUT"
echo "Wrote $OUTPUT"
