# RNApipeline

## Start here

For a normal cohort run, use **`mainRNAseq.nf`** and the included cohort config.
Build a reference only when the required release is not already available; that
uses **`buildReference.nf`** and its separate example config.

```bash
nf-run cohortA /biojobs/pipelines/RNApipeline/mainRNAseq.nf \
  -profile pam_cpu \
  -c /biojobs/pipelines/RNApipeline/configExamples/cohort.config \
  --sampleList /home/you/cohortA/samples.txt \
  --ref /ref/ensembl_v107_GRCh38_rna2
```

To build that reference once:

```bash
nf-run reference-ensembl-v107 /biojobs/pipelines/RNApipeline/buildReference.nf \
  -profile pam_cpu \
  -c /biojobs/pipelines/RNApipeline/configExamples/reference-build.config
```

For the example protocol, run with this config as shown. For a different
protocol, copy it and change the strand settings. Edit the reference-build
config if your reference destination differs.
The examples assume the pipeline checkout is at `/biojobs/pipelines/RNApipeline`;
replace that path if yours differs. See [outputs and metrics](#outputs-and-metrics),
[container builds](CONTAINERS.md) and the [hotspot reference files](reference/README.md).

## Where to change things

| Change | File or folder |
|---|---|
| Connect or configure the standard analysis workflow | `mainRNAseq.nf` |
| Add or change an RNA process | `modules/rna_seq_pipe.nf`; connect it in `mainRNAseq.nf` |
| Add or change Kallisto processes | `modules/kallisto.nf`; connect them in `mainRNAseq.nf` |
| Add or change variant processes | `modules/mpileup.nf` or `modules/variant_calling.nf`; connect them in `mainRNAseq.nf` |
| Add or change metrics, Somalier or hotspot-output processes | `modules/sample_outputs.nf`; connect them in `mainRNAseq.nf` |
| Add or change input checks or the run summary | `modules/check_prepost_pipeline.nf`; connect them in `mainRNAseq.nf` |
| Add or change a reference-build process | `modules/index.nf`, then connect it in `buildReference.nf` |
| Edit helper code called by a process | `scripts/` for analysis helpers; `reference/` for reference-build helpers |
| Change a container | `Docker/` |
| Select cohort or reference-build parameters | `configExamples/` |

`buildReference.nf` is the runnable reference workflow. The scripts beside it
under `reference/` are helper code used by its processes; they are not separate
pipelines. `scripts/` contains helpers used by the analysis and optional
strandedness preflight.

The pipeline folders are `Docker/`, `configExamples/`, `modules/`, `reference/`
and `scripts/`. Maintainer checks are kept separately under the hidden `.github/`.

## Outputs and metrics

| Path | Contents |
|---|---|
| `<run>_pipeline_summary.json` | Run parameters, samples/inputs, pipeline and reference identities, containers and completion checks |
| `FeatureCounts_output/` | Per-sample `exonscount.txt` and `genecount.txt` files in the existing table formats |
| `UnmappedReads_output/` | Optional unmapped R1/R2 FASTQs |
| `Metrics_output/<sample>.metrics.tsv` | One normalized, project-aggregatable metrics table per sample |
| `Kallisto_output/<sample>/abundance.tsv` | When Kallisto is enabled |
| `Somalier_output/<sample>.somalier` | When Somalier is enabled (default) |
| `Hotspot_output/<sample>.hotspots.tsv` | When variant calling is enabled (default off) |

FastQC HTML/ZIP reports, STAR logs, FeatureCounts assignment summaries,
Kallisto run info, BAM/BAI, BCF/VCF, VEP annotations and depth remain in the
external Nextflow task cache. The workflow publishes no cohort matrices or
standalone MD5, status or sample-check files.

### Metrics schema

Each per-sample file has one header and one row per metric, using these columns:

`sample_id`, `metric_scope`, `fastq_id`, `read_role`, `metric_source`, `metric`,
`value`, `unit`, `qc_status`.

Additional lanes add FASTQ rows with the same columns. `fastq_id` is the original
FASTQ basename; read roles and order come from the sample sheet. Sample-level
rows have empty FASTQ/read-role fields. Values are scalar numbers or strings;
timestamps, versions and statuses remain strings. No FASTQ quality distributions
are retained.

| Source | Collected metrics |
|---|---|
| FastQC, each FASTQ | Version; every module status; Basic Statistics summary fields except file name/type; minimum and mean of per-base mean quality; weighted mean per-sequence quality; minimum/maximum read lengths; maximum per-base N%; maximum base-content deviation from 25%; mean sequence GC%; maximum absolute tile quality deviation; deduplicated%; number and summed percentage of overrepresented sequences; maximum adapter%; number of reported kmers and maximum observed/expected ratio |
| STAR | Every populated scalar field from `Log.final.out`, including mapped/unmapped/splice/chimeric counts/rates and run timing |
| FeatureCounts | Every assignment category, count and percentage, separately for exon and gene summaries; number of gene-table rows with count > 0 |
| Kallisto | Enabled/disabled status; `n_targets`, `n_bootstraps`, `n_processed`, `n_pseudoaligned`, `n_unique`, `p_pseudoaligned`, `p_unique`, `kallisto_version`, `index_version` when emitted by the tool |

FastQC metrics retain pass/warn/fail status. STAR, FeatureCounts and Kallisto
metrics have no invented pass/fail threshold. Metrics join by sample ID, and a
completion check covers every input sample.

RNApipeline 2.0 keeps the 1.9 sample sheet, FASTQ lane grouping, STAR alignment
and FeatureCounts counting behavior. It adds per-sample QC metrics and Somalier
fingerprints. Kallisto remains optional; hotspot calling remains opt-in. The
sample-first metrics can be combined outside the pipeline for project analysis.

This development version is on `pipeline-2.0`; `tod-nf` remains the 1.9 baseline.

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

[`configExamples/cohort.config`](configExamples/cohort.config) is a
reverse-stranded paired-end example. Copy it for your cohort and set the
strand options to match the library protocol before running.

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
the generic merged exon BED for legacy reference consumers. The 2.0 bundle
contains only assets used by its supported alignment, quantification, fingerprint
and hotspot tasks.

Somalier sites and GI hotspot targets are separate assets. The builder normalizes
Somalier chromosome names to the FASTA and checks every REF base. Hotspot records
that cannot be resolved are retained in `hotspots/mapping_audit.tsv`; they do not
become reassuring no-call rows. This is a research catalogue, not a clinically
validated panel. See its [scope and evidence limitations](reference/README.md).

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
before running on compute workers (see [container builds](CONTAINERS.md)).

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

The 2.0 workflow has no cohort aggregation or standalone raw-depth task. Raw tool
reports and intermediates remain in the Nextflow task cache for troubleshooting
and resume. With `nf-run` on IAC, the launcher stores it under the run directory
shown below. For direct Nextflow runs, the default cache is outside the checkout
under `~/.cache/nextflow/RNApipeline`; set `NF_RUN_WORK_DIR` to use other storage.
Retain the cache while you may resume.

## Execution and resume

Choose one execution profile: `pam_cpu`, `pam_gpu`, or `burst`. Burst mode is
available only when an administrator has enabled it. See the IAC
[Nextflow run guide](https://wiki.genehetx.fr/compute/nf-run/) for profiles,
input visibility on workers, and infrastructure troubleshooting.

Each IAC `nf-run` name has its own directory:

```text
/biojobs/nextflow/<run-name>/{results,logs,metadata,work}
```

Results and logs are under `results/` and `logs/`; Nextflow's task cache is in
the external `work/` directory. This is runtime data, not part of the pipeline
checkout. Keep it to resume. The first invocation starts fresh. Reusing a
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
