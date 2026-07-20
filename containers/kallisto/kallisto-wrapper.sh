#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "$(uname -m)" in
    aarch64|arm64) IMAGE="${ROOT_DIR}/arm64/kallisto-0.51.1.sif" ;;
    x86_64) IMAGE="${ROOT_DIR}/amd64/kallisto-0.51.1.sif" ;;
    *) echo "Unsupported host architecture: $(uname -m)" >&2; exit 1 ;;
esac

APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"
command -v "$APPTAINER_BIN" >/dev/null 2>&1 || {
    echo "Apptainer is required" >&2
    exit 1
}
[[ -s "$IMAGE" ]] || { echo "Missing Kallisto SIF: $IMAGE" >&2; exit 1; }

exec "$APPTAINER_BIN" exec --no-home \
    --bind "$PWD:$PWD" \
    --bind /ref:/ref \
    --bind /biojobs:/biojobs \
    "$IMAGE" kallisto "$@"
