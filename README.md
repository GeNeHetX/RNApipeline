# RNApipeline

RNApipeline is a Nextflow workflow for RNA-seq quality control, STAR alignment,
gene/exon counts, transcript quantification, and MultiQC reporting. Variant
calling and annotation are optional.

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
The STAR BAM sort cap defaults to 30 GB; if you lower STAR memory below that,
lower `params.limitBAMsortRAM` in the same config as well.
For a non-IAC run, set `params.ref` to the reference directory in the cohort
config.

## Run

Pass absolute paths for the pipeline, sample sheet, and cohort config. These
files may live in any directory you can access:

```bash
nf-run cohortA /biojobs/pipelines/RNApipeline/fullPairedEnd.nf \
  -profile pam_cpu \
  -c /home/you/cohortA/analysis.config \
  --sampleList /home/you/cohortA/samples.txt
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
nf-run cohortA /biojobs/pipelines/RNApipeline/fullPairedEnd.nf \
  -profile pam_cpu -resume \
  -c /home/you/cohortA/analysis.config \
  --sampleList /home/you/cohortA/samples.txt
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
  --sampleList /home/you/cohortA/samples.txt
```

Read the summary, set the chosen `FeatureCountStrand` and `kallistoStrand` in
the cohort analysis config, then launch production under a separate name such
as `cohortA`. The two runs have separate work caches, logs, and results.

For operational details, cleanup, burst mode, and troubleshooting, use the
IAC [Nextflow run guide](https://wiki.genehetx.fr/compute/nf-run/) and its
[pipeline/runtime reference](https://wiki.genehetx.fr/reference/pipeline-and-nextflow-reference/).
