# Hotspot catalogue

`catalogue.tsv` preserves the existing 193 GI entries and 46 genes, including
all original KLF4 IPMN entries. It has no genome coordinates assigned by this
update. `source/gi_hotspots_catalogue.tsv` is the unchanged source retrieved from
[the project handoff](https://drive.google.com/drive/folders/15t7JZVdUC7kS1SZxPPRu7YIFqQU1uAug).
Its SHA256 is recorded in `source/source.sha256`. The pan-cancer additions in
Worknotes remain proposed and have not been added to this source snapshot.

Prepare/audit this snapshot without a reference or reads:

```bash
python3 PrePostScripts/prepare_catalogue.py \
  --source reference/hotspots/source/gi_hotspots_catalogue.tsv \
  --outdir reference/hotspots
```

The preparation checks IDs, gene-symbol syntax and agreement between protein
substitution fields. All 46 gene symbols match approved HGNC symbols in the
dated `hgnc_symbol_audit.tsv`; `hgnc_source.json` records the official source and
its checksum. `catalogue_audit.tsv` records event limitations and the
original evidence labels. These checks do not establish that every candidate
allele has been observed biologically. 111 entries have no source URL. The shared KRAS ClinVar link identifies G12D in
nevus sebaceous and does not substantiate the other 13 entries sharing it or their
tumour context. `source_review.tsv` records the source-access review; unresolved
biological evidence remains explicit. Some citations support a gene or mechanism
rather than each listed allele.

`buildReference.nf` then selects one transcript per gene: MANE Select, Ensembl
canonical, APPRIS principal 1, basic, longest CDS, then transcript ID as a stable
tie-break. `--transcript_overrides` accepts a TSV with `gene` and `transcript_id`
(versioned IDs are preferred). A residue mismatch stays unresolved; the builder
does not search alternate isoforms just to find the requested amino acid.
Incomplete CDS, wrong phase, internal stops, missing transcripts and unsupported
events are audited. Candidate coding SNVs are enumerated from the reference
codon and independently checked by same-release VEP before publication.

The reference contains `hotspots/targets.json`, `targets.bed`, the source
catalogue, candidate alleles, mapping audit and VEP validation result. Possible
protein-derived alleles remain candidates; VEP agreement is a technical
coordinate/consequence check, not clinical validation. Codon-only requests
produce codon targets; unspecified indels/repeats, fusions, amplification and
TERT promoter DNA targets are outside this report. Broad-inactivation targets
use selected CDS for the source gene-mechanism entries, plus TP53 as specified
in Worknotes. Only frameshift or stop-gained variants qualify in broad regions.

Variant calling remains disabled by default. Enabling it reports protein-altering
variants in validated hotspot codons and the supported broad truncating events.
Unassessable genes stay in the reference audit and are omitted from the report.

The six report columns are:

`sample_id`, `genesymbol`, `mutated_reads`, `total_reads`, `cdna_mutation`, `protein_mutation`.

Detected rows use FORMAT/AD for ALT-supporting reads and FORMAT/DP for total
local depth. The same BCFtools pileup supplies coverage for no-call rows, with
missing target positions treated as zero. Hotspot and dual-category genes use
minimum hotspot coverage; broad-only genes use median coding coverage. All
include zero-depth bases. A no-call is not a DNA wild-type determination; RNA
expression, allele-specific expression and NMD can hide variants.
