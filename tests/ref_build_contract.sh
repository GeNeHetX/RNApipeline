#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

command -v nextflow >/dev/null
NEXTFLOW_HOME="${NXF_HOME:-$HOME/.nextflow}"

grep -q "final_manifest.isFile()" "$ROOT_DIR/PrePostScripts/ref_build.nf"
grep -q "params.check_only || (!params.force" "$ROOT_DIR/PrePostScripts/ref_build.nf"
grep -q 'params.reference_id ?: "ensembl_v\${params.ensembl_release}_GRCh38"' "$ROOT_DIR/PrePostScripts/ref_build.nf"
grep -q "rnapipeline-buildref-v107" "$ROOT_DIR/README.md"
grep -q "PREPARE_REFERENCE_INPUTS" "$ROOT_DIR/PrePostScripts/ref_build.nf"
grep -q "FINALIZE_REFERENCE" "$ROOT_DIR/PrePostScripts/ref_build.nf"

mkdir -p "$WORK_DIR/complete/VEP/homo_sapiens/107_GRCh38"
printf 'reference\n' > "$WORK_DIR/complete/ref.fa"

python3 - "$WORK_DIR/complete/reference_manifest.json" <<'PY'
import json
import pathlib
import sys

manifest = {
    "manifest_schema_version": 2,
    "reference_id": "complete",
    "species": "homo_sapiens",
    "assembly": "GRCh38",
    "ensembl_release": 107,
    "files": [{"path": "ref.fa", "size_bytes": pathlib.Path(sys.argv[1]).parent.joinpath("ref.fa").stat().st_size}],
    "complete": True,
}
pathlib.Path(sys.argv[1]).write_text(json.dumps(manifest) + "\n")
PY

complete_log="$WORK_DIR/complete.log"
NXF_OFFLINE=true NXF_HOME="$NEXTFLOW_HOME" nextflow run \
  "$ROOT_DIR/PrePostScripts/ref_build.nf" \
  -work-dir "$WORK_DIR/complete-work" \
  -ansi-log false \
  --reference_root "$WORK_DIR" \
  --reference_id complete \
  >"$complete_log" 2>&1

grep -q 'CHECK_REFERENCE' "$complete_log"
! grep -q 'PREPARE_REFERENCE_INPUTS' "$complete_log"

mkdir -p "$WORK_DIR/derived/ensembl_v108_GRCh38/VEP/homo_sapiens/108_GRCh38"
printf 'reference\n' > "$WORK_DIR/derived/ensembl_v108_GRCh38/ref.fa"
python3 - "$WORK_DIR/derived/ensembl_v108_GRCh38/reference_manifest.json" <<'PY'
import json
import pathlib
import sys

manifest_path = pathlib.Path(sys.argv[1])
manifest = {
    "manifest_schema_version": 2,
    "reference_id": "ensembl_v108_GRCh38",
    "species": "homo_sapiens",
    "assembly": "GRCh38",
    "complete": True,
    "ensembl_release": 108,
    "files": [{"path": "ref.fa", "size_bytes": manifest_path.parent.joinpath("ref.fa").stat().st_size}],
}
manifest_path.write_text(json.dumps(manifest) + "\n")
PY

derived_log="$WORK_DIR/derived.log"
NXF_OFFLINE=true NXF_HOME="$NEXTFLOW_HOME" nextflow run \
  "$ROOT_DIR/PrePostScripts/ref_build.nf" \
  -work-dir "$WORK_DIR/derived-work" \
  -ansi-log false \
  --reference_root "$WORK_DIR/derived" \
  --ensembl_release 108 \
  >"$derived_log" 2>&1

grep -q 'CHECK_REFERENCE' "$derived_log"
! grep -q 'PREPARE_REFERENCE_INPUTS' "$derived_log"

printf '%s\n' 'reference build contract tests passed'
