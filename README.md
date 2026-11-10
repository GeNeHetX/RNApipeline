# RNApipeline

RNApipeline 2.0 keeps the 1.9 sample sheet, FASTQ lane grouping, STAR alignment
and FeatureCounts counting behavior. It adds per-sample QC metrics and Somalier
fingerprints. Kallisto remains optional; hotspot calling remains opt-in.
MultiQC and cohort matrices are replaced by metrics that can be combined later.

This development version is on `pipeline-2.0`; `tod-nf` remains the 1.9 baseline.
See [completed validation and remaining deployment work](docs/validation.md),
[the output and metric schema](docs/outputs.md) and [the hotspot catalogue](reference/hotspots/README.md).

## What you need

- FASTQ files readable from the compute workers.
- A reference directory with the indexes required by the enabled tools. On IAC,
  the site config supplies the shared `/ref` root and Ensembl 107 by default.
- A headerless, whitespace-separated sample sheet with four fields per row:

  ```text
  sample_id  read1.fastq.gz  read2.fastq.gz  input_directory
  ```

Use one row per read pair. Repeat the sample ID for each lane; lane rows are
grouped into one sample and kept in sheet order. FASTQ names are relative to
the input directory. Use `-` for read 2 with single-end data, and set
`params.single_end = true` in the analysis config. Blank lines and lines
starting with `#` are ignored. Do not put spaces inside a field.

For example, two paired-end lanes from one sample:

```text
sampleA sampleA_L001_R1.fastq.gz sampleA_L001_R2.fastq.gz /data/cohortA
sampleA sampleA_L002_R1.fastq.gz sampleA_L002_R2.fastq.gz /data/cohortA
sampleB sampleB_L001_R1.fastq.gz sampleB_L001_R2.fastq.gz /data/cohortA
```

## Cohort analysis config

Create a small Nextflow config beside the sample sheet. This example is for a
reverse-stranded paired-end library; select values for the actual library
protocol before running:

```groovy
params.single_end = false
params.FeatureCountStrand = 2
params.kallistoStrand = ' --rf-stranded '
params.ensembl_release = 107
params.kallisto = true
params.variant_calling = false
params.somalier = true
params.bcftools_max_depth = 100000
```

For inward paired-end reads or single-end reads, Salmon's `ISF`/`SF` evidence
corresponds to Kallisto `--fr-stranded` and FeatureCounts `-s 1`; `ISR`/`SR`
corresponds to `--rf-stranded` and `-s 2`; `IU`/`U` is unstranded (no Kallisto
strand flag and FeatureCounts `-s 0`). Check the per-sample and batch evidence
before choosing. STAR itself does not take a strandness setting. See the
optional preflight below when the protocol is unknown. This mapping follows
the tools' [Salmon library-type definitions](https://salmon.readthedocs.io/en/latest/library_type.html),
[Kallisto strand options](https://munfred.github.io/kallisto/manual), and
[FeatureCounts strand options](https://subread.sourceforge.net/featureCounts.html).

Tool and scientific defaults live in `nextflow.config`. On IAC, `/etc/nextflow/site.config`
sets the execution profile, shared paths, default resources, containers, and
run output location. You do not need to copy or edit an IAC config. The IAC
defaults are 6 threads for STAR, Kallisto, Salmon preflight, and reference-index
build, and 36 GB for STAR. Override these for one run by adding
`params.max_threads` or `params.max_star_memory` to the cohort config.
`mainRNAseq.nf` takes a prepared reference through `--ref` or the existing
`reference_root`/`ensembl_release` convention. It does not build references while
mapping reads. Build a fresh 2.0 human GRCh38 reference with
[`buildReference.nf`](buildReference.nf); existing 1.9 bundles lack the Somalier
panel, hotspot targets and reference manifest. To use an older reference for
alignment/quantification only, set `params.somalier = false` and keep
`params.variant_calling = false`.

The reference builder has one release parameter: `--ensembl_release` (default
107). It downloads that release's primary-assembly FASTA, chromosome GTF, cDNA
and indexed VEP cache. It builds the FASTA, STAR and Kallisto indexes, prepares
the Somalier RNA panel, maps/audits the existing hotspot catalogue and checks
candidate consequences using same-release VEP. It records source URLs, SHA256,
versions and STAR index settings in `reference_manifest.json`. It also retains
the generic merged exon BED for legacy reference consumers. GATK-only known
variants and dictionaries are not part of the 2.0 bundle.

Somalier sites and GI hotspot targets are separate assets. The builder normalizes
Somalier chromosome names to the FASTA and checks every REF base. Hotspot records
that cannot be resolved are retained in `hotspots/mapping_audit.tsv`; they do not
become reassuring no-call rows. This is a research catalogue, not a clinically
validated panel. See its [scope and evidence limitations](reference/hotspots/README.md).

The builder publishes a new directory only after all indexes and checks succeed.
Choose a name that does not already exist. On IAC:

```bash
nf-run reference-ensembl-v107-rna2 /biojobs/pipelines/RNApipeline/buildReference.nf \
  -profile pam_cpu \
  --ensembl_release 107 \
  --reference_build_root /ref \
  --reference_name ensembl_v107_GRCh38_rna2
```

The account needs write access to `/ref`. For another human GRCh38 release,
change `--ensembl_release` and the destination name. The selected release must
have matching assets available on Ensembl; it also selects the VEP image tag.
The pipeline and matching VEP images must be installed/published by the maintainer
before running on compute workers (see [container builds](docs/containers.md)).

The index uses `sjdbOverhang=100`, preserving the existing protocol: STAR defines
this as the maximum read length minus one for splice-junction database sequences
([STAR manual](https://github.com/alexdobin/STAR/blob/2.7.11b/doc/STARmanual.pdf)).
This is a fixed index setting, recorded in the reference manifest; the update
does not infer a new index setting from each cohort's FASTQs.

Paired-end FeatureCounts runs use `-p --countReadPairs -O`: each fragment is
counted as a pair, and a fragment overlapping multiple genes or exons is counted
for each overlap. This can make counts across features non-exclusive. The
pipeline does not set `-B`, so counting does not require both mates to map.

`params.variant_calling = true` enables BCFtools mpileup, multiallelic calling,
normalization, VEP and one six-column hotspot TSV per sample. Targets come from
the prepared reference. `params.bcftools_max_depth` is the per-position/per-BAM
cap (default 100000); it affects read counts when reached. The caller prior
remains `params.bcftools_prior = 1.1e-3`. Base quality defaults to 13 and mapping
quality to 0, as in the previous caller; they are exposed as
`variant_min_base_quality` and `variant_min_mapping_quality`. Depth and ALT counts
come from the same filtered pileup. BCFtools' default BAQ, overlap and duplicate
handling are retained ([BCFtools manual](https://samtools.github.io/bcftools/bcftools.html)).
Do not toggle mpileup and VEP separately. No RNA no-call establishes DNA wild type.

Set `params.kallisto = false` to skip quantification or `params.somalier = false`
for alignment-only work with a legacy/nonhuman reference. Single-end Kallisto
uses `params.fragment_length` (mean fragment length, default 120) and
`params.read_sd` (default 20); legacy `read_len` overrides still work.
Somalier publishes one fingerprint per sample; cohort relatedness is left for
project analysis. Its [RNA sites panel](https://github.com/brentp/somalier/releases)
uses expressed sites and still needs informativeness assessment on your GI data.

The old `aggregate` and `samtools_depth` switches are accepted with a warning;
2.0 does not invoke aggregation or publish raw depth. Raw tool reports and
intermediates remain in Nextflow work for troubleshooting/resume. Keep the work
area until you deliberately clean it after a successful run.

## Run

Pass absolute paths for the pipeline, sample sheet, and cohort config. These
files may live in any directory you can access:

```bash
nf-run cohortA /biojobs/pipelines/RNApipeline/mainRNAseq.nf \
  -profile pam_cpu \
  -c /home/you/cohortA/analysis.config \
  --sampleList /home/you/cohortA/samples.txt \
  --ref /ref/ensembl_v107_GRCh38_rna2
```

Choose one execution profile: `pam_cpu`, `pam_gpu`, or `burst`. Burst mode is
available only when an administrator has enabled it. See the IAC
[Nextflow run guide](https://wiki.genehetx.fr/compute/nf-run/) for profiles,
input visibility on workers, and infrastructure troubleshooting.

Each run name has its own directory:

```text
/biojobs/nextflow/<run-name>/{results,logs,metadata,work}
```

Results and logs are under `results/` and `logs/`; Nextflow's task cache is in
`work/`. Keep `work/` to resume. The first invocation starts fresh. Reusing a
run name without `-resume` is refused. To resume after an interruption, repeat
the same command and add `-resume`:

```bash
nf-run cohortA /biojobs/pipelines/RNApipeline/mainRNAseq.nf \
  -profile pam_cpu -resume \
  -c /home/you/cohortA/analysis.config \
  --sampleList /home/you/cohortA/samples.txt \
  --ref /ref/ensembl_v107_GRCh38_rna2
```

Use a new run name for a new cohort or a deliberate change to the sample sheet
or analysis settings. The launcher retains each invocation's command, config
and sample-sheet snapshots, and pipeline revision under `metadata/`.

## Optional strandedness preflight

Run this only when library orientation is unknown. It uses the same sample
sheet, runs Salmon's automatic library inference, preserves each sample's raw
Salmon output, and writes a per-sample plus batch summary to
`results/strandedness_summary.tsv`. It reports evidence; it does not choose a
protocol or alter the production config. Salmon's format histogram counts
observed mappings, and a fragment may contribute to more than one format; the
reported percentages are normalized within that histogram, not percentages of
all input fragments. Salmon documents the fields in its
[`lib_format_counts.json` output](https://salmon.readthedocs.io/en/latest/file_formats.html#observed-library-format-counts).

```bash
nf-run cohortA-strand /biojobs/pipelines/RNApipeline/strandPreflight.nf \
  -profile pam_cpu \
  --transcriptome /ref/ensembl_v107_GRCh38/cdna.fa \
  --sampleList /home/you/cohortA/samples.txt \
  --ref /ref/ensembl_v107_GRCh38_rna2
```

Read the summary, set the chosen `FeatureCountStrand` and `kallistoStrand` in
the cohort analysis config, then launch production under a separate name such
as `cohortA`. The two runs have separate work caches, logs, and results.

For operational details, cleanup, burst mode, and troubleshooting, use the
IAC [Nextflow run guide](https://wiki.genehetx.fr/compute/nf-run/) and its
[pipeline/runtime reference](https://wiki.genehetx.fr/reference/pipeline-and-nextflow-reference/).
