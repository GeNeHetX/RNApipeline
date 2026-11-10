# RNApipeline 2.0 validation — 2026-11-10

Implemented on `pipeline-2.0`, based on `tod-nf` commit `13b78c1` (1.9 line).
The sample sheet, lane grouping, alignment defaults, counting semantics and
optional Kallisto behavior are preserved. The source is ready for review; this
record does not constitute a production release or scientific validation.

## Passed

- Real Nextflow 24.04.3 wiring checks: paired/single-end, multiple lanes, legacy
  CSV input, reference-root derivation, Kallisto disabled, opt-in variants,
  retained aggregate-switch compatibility, sample identity, missing input and
  deprecated-switch errors, strandedness preflight, complete reference-builder
  stages and publication, plus the retained aggregation script: 18 cases.
- Reference/report regressions: forward/reverse strand codons, wrong residue
  auditing, VEP consequence mismatch rejection, ALT counts, zero-inclusive depth,
  broad median, dual-category hotspot minimum even when coding median is lower,
  and omission of an unassessable dual-category no-call.
- Active pipeline and VEP107 Dockerfiles built successfully for ARM64 and x86_64.
- Real FastQC, FeatureCounts, Kallisto, Somalier and BCFtools ran on synthetic
  directly constructed FASTQ/BAM inputs in both pipeline images. No real reads
  were mapped. Somalier emitted `S.somalier`; metrics parsed the real FastQC ZIP;
  a normalized synthetic variant retained AD=40,40 and DP=80.
- VEP107 on ARM64 loaded the matching FASTA through Bio::DB::HTS::Faidx.
- All 46 source catalogue symbols matched approved HGNC symbols. Source IDs,
  records and KLF4 entries were preserved; evidence limitations remain audited.
- The default release-107 FASTA, GTF, cDNA and indexed-VEP-cache URLs returned
  HTTP 200. The VEP download alone is 22,711,311,448 compressed bytes.
- Whitespace checks and Python compilation passed. The final reference assembly
  stages links in work; a Nextflow publication check confirmed copied output
  files are ordinary files, avoiding an extra full reference copy in work.

## Local image identities

These are local image IDs, not published registry manifest digests.

| Image | Architecture | Image ID |
|---|---|---|
| `genehetx/genehetx-rnaseq:v2.0.0` | ARM64 | `sha256:6d19d3e062c9e20515277162015d4cd1c403dcabfc84882a745e8b34cd8d530f` |
| `genehetx/genehetx-rnaseq:v2.0.0-amd64` | x86_64 | `sha256:bccbafbb6fd178eab3b0106e702b3654dd54008dbe2c9b6011b48e497319302b` |
| `genehetx/genehetx-vep:107` | ARM64 | `sha256:9d93aeaa8c2f027c405e1b239456beaea504872bd0b67c003ee5895375f6a467` |
| `genehetx/genehetx-vep:107-amd64` | x86_64 | `sha256:4128567beb75e7e765105ae0d0e7ddf1df84eede8dc8a132fcdf97161e82b609` |

The x86_64 tool execution was tested under Docker's emulation on this ARM64 host;
native x86_64 production performance has not been measured. Architecture suffixes
here keep both local images available during checks; production tags must carry
the matching platform image or a multiarchitecture manifest.

## Still required before production

Complete a human release-107 reference build on a host with sufficient memory
and work disk, then run representative GI samples. The local Docker engine has
about 8 GB RAM, below the 36 GB human STAR-build request. No complete human
reference bundle, sensitivity measurement or GI fingerprint informativeness
assessment was performed here. No image or code branch has been published by
these checks.

From a checkout of this version on a Docker host with at least the requested
memory, run as the reference maintainer with write access to the destination:

```bash
mkdir -p /ref
nextflow run buildReference.nf -with-docker \
  --ensembl_release 107 \
  --reference_build_root /ref \
  --reference_name ensembl_v107_GRCh38_rna2
```

Install the images first using [the build commands](containers.md). For IAC, the
README gives the equivalent `nf-run` invocation; the source checkout and site
container settings must point to this version before invoking it.
