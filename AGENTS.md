# RNApipeline agent notes

## Scope

- This repository is the source of truth for the RNA-seq Nextflow pipeline.
- The IAC repository may be inspected for runtime contracts, but changes for this
  work belong only in RNApipeline.
- Preserve the existing IFB and Google Cloud execution paths unless a change is
  explicitly additive and does not alter pipeline outputs.

## Compatibility rules

- Keep pipeline version and reference version separate in run metadata.
- Do not refactor the workflow or rename existing output files while building
  the TOD reference.
- Existing precomputed-reference behavior uses the legacy `kalliso_index`
  filename; compatibility aliases are safer than changing the module contract.
- The current production reference tuple is Ensembl 107 genome/GTF/cDNA,
  Ensembl 105 known variants, and Ensembl 113 VEP cache.
- TOD-specific execution must use `/ref` for reference data, Slurm, and
  Apptainer. Do not add TOD paths to the shared IAC repository from here.

## Reference builds

- Use `PrePostScripts/ref_build_tod.sh` for the TOD/PAM reference build.
- Build in Slurm scratch and publish only a complete, validated directory to
  `/ref/ensembl_v107_GRCh38`.
- Never commit generated reference data, container images, Slurm logs, or
  Nextflow work directories.
- Keep `reference_manifest.json` with every published reference. It records
  source versions, tool/container versions, and artifact checksums.

## Working-tree safety

- Existing uncommitted changes belong to the user. Inspect diffs before editing
  overlapping files and do not reset, discard, or stage unrelated changes.
- Prefer a small additive patch over broad formatting or line-ending rewrites.
- Run `bash -n` on shell changes, validate Nextflow configuration, and use
  `git diff --check` before handing work back.
