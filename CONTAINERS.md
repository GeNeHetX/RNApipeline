# Building RNApipeline 2.0 containers

Run these commands from the RNApipeline checkout on a Docker host. They build
local images; workers need the same images available through your registry or
container cache before using the pipeline. No registry image is published by
editing this repository.

```bash
docker buildx build --platform linux/arm64 --load \
  -t genehetx/genehetx-rnaseq:v2.0.0 \
  -f Docker/Dockerfile.rnaseq .

docker buildx build --platform linux/arm64 --load \
  -t genehetx/genehetx-vep:107 --build-arg ENSEMBL_RELEASE=107 \
  -f Docker/Dockerfile.vep .
```

Use `linux/amd64` on x86_64. The main image retains STAR 2.7.11b, Kallisto 0.51.0,
FeatureCounts 2.0.2, FastQC 0.12.1, samtools 1.13 and BCFtools 1.21. Somalier 0.3.7
is added using native, explicitly selected Bioconda builds. The checked-in architecture-specific lock files pin every package URL and checksum.
The specification is also recorded inside `/opt/somalier/explicit-spec.txt`; pysam 0.23.3
supports reference preparation.
Somalier's libraries are scoped to its wrapper and do not replace the existing
RNA tools' runtime libraries.

VEP remains separate, using the selected Ensembl release for its API and cache.
The old mutable `latest` Dockerfile is replaced by a parameterized native build.
For another release, change the VEP build argument and tag to the same number;
`--ensembl_release` selects that tag in both workflows. No separate cache-release
parameter is needed. Production deployments should pin the published images by
digest in the site config; the defaults name the new images and do not prove
those tags have been published to a registry.

A reference can also be built directly with Nextflow:

```bash
mkdir -p /absolute/path/to/references
nextflow run buildReference.nf -with-docker \
  --ensembl_release 107 \
  --reference_build_root /absolute/path/to/references \
  --reference_name ensembl_v107_GRCh38_rna2
```

Allow enough work disk for downloads, indexes and the final assembled bundle.
The final reference is staged using links and copied to its published destination only after VEP validation. A full human
STAR build uses the configured 36 GB memory request. This is a substantial build;
the repository's small fixture checks do not replace completing it on the target
compute environment before production mapping.

Local checks (Python needs pysam 0.23.3):

```bash
python3 .github/check_hotspot_outputs.py
python3 .github/check_sample_identity.py
```

The workflow checks use mock tools to exercise real Nextflow scheduling, sample
identity, output publication and all reference-builder stages. They do not
validate biological sensitivity or a complete human reference.

See [the validation record](VALIDATION.md) for the checks actually completed and the remaining full-reference/production work.
