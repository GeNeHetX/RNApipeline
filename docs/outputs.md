# RNApipeline 2.0 outputs

| Path | Contents |
|---|---|
| `<run>_pipeline_summary.json` | Run parameters, samples/inputs, pipeline and reference identities, containers and completion checks |
| `FeatureCounts_output/` | Existing per-sample `exonscount.txt` and `genecount.txt` filenames and table formats |
| `UnmappedReads_output/` | Existing optional unmapped R1/R2 FASTQs |
| `Metrics_output/<sample>.metrics.tsv` | One normalized, project-aggregatable metrics table per sample |
| `Kallisto_output/<sample>/abundance.tsv` | When Kallisto is enabled |
| `Somalier_output/<sample>.somalier` | When Somalier is enabled (default) |
| `Hotspot_output/<sample>.hotspots.tsv` | When variant calling is enabled (default off) |

FastQC HTML/ZIP, STAR logs, FeatureCounts assignment summary files, Kallisto run
info, BAM/BAI, BCF/VCF, VEP annotation and depth stay in Nextflow work. No MultiQC,
cohort matrices, standalone MD5/status/sample-check files are published.

## Metrics schema

`sample_id`, `metric_scope`, `fastq_id`, `read_role`, `metric_source`, `metric`,
`value`, `unit`, `qc_status`.

There is one header and one row per metric. A sample with additional lanes has
more FASTQ rows, with the same columns. `fastq_id` is the original FASTQ basename;
read roles and order come from the sample sheet. Run metadata retains the full
input association. Sample-level rows have empty FASTQ/read-role fields.
Concatenate files keeping one header, or import them as tables and join on
sample ID plus metric identity. Values are scalar numbers or strings; timestamps,
versions and statuses remain strings. No FASTQ quality distributions are stored.

| Source | Collected metrics |
|---|---|
| FastQC, each FASTQ | Version; every module status; Basic Statistics summary fields except file name/type; minimum and mean of per-base mean quality; weighted mean per-sequence quality; minimum/maximum read lengths; maximum per-base N%; maximum base-content deviation from 25%; mean sequence GC%; maximum absolute tile quality deviation; deduplicated%; number and summed percentage of overrepresented sequences; maximum adapter%; number of reported kmers and maximum observed/expected ratio |
| STAR | Every populated scalar field from `Log.final.out`, including mapped/unmapped/splice/chimeric counts/rates and run timing |
| FeatureCounts | Every assignment category, count and percentage, separately for exon and gene summaries; number of gene-table rows with count > 0 |
| Kallisto | Enabled/disabled status; `n_targets`, `n_bootstraps`, `n_processed`, `n_pseudoaligned`, `n_unique`, `p_pseudoaligned`, `p_unique`, `kallisto_version`, `index_version` when emitted by the tool |

FastQC summary metrics retain the source module's pass/warn/fail status. STAR,
FeatureCounts and Kallisto numeric metrics have no invented pass/fail threshold.
Metrics join by sample ID and a final completion check covers every input sample.
