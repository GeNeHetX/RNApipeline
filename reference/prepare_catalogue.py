#!/usr/bin/env python3
"""Prepare the source catalogue without assigning genome coordinates."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re


def prepare(source, outdir):
    source, out = Path(source), Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    with source.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        rows = list(reader)
    ids = set()
    audit = []
    for r in rows:
        if not r['id'] or r['id'] in ids:
            raise ValueError('Missing or duplicate catalogue ID: ' + r['id'])
        ids.add(r['id'])
        if not re.fullmatch(r'[A-Za-z0-9-]+', r['gene']):
            raise ValueError('Invalid gene symbol: ' + r['gene'])
        kind = r['variant_class']
        if kind in ('protein_snv', 'protein_nonsense'):
            m = re.fullmatch(r'([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY*])', r['feature'])
            if not m or m[1] == m[3] or int(m[2]) <= 0:
                raise ValueError('Invalid protein substitution: ' + r['id'])
            if (r['reference_aa'], r['aa_position'], r['alternate_aa']) != m.groups():
                raise ValueError('Inconsistent protein fields: ' + r['id'])
        status = 'requires_reference_mapping'
        if kind in ('fusion', 'amplification', 'promoter_snv', 'coding_repeat'):
            status = 'outside_standard_RNA_SNV_assay'
        elif kind in ('protein_indel', 'frameshift_unspecified'):
            status = 'requires_exact_nucleotide_allele'
        audit.append(dict(id=r['id'], genesymbol=r['gene'], feature=r['feature'],
                          status=status, evidence=r['source_allele_specificity'], source_url=r['source_url'],
                          evidence_review=('missing_source_url' if not r['source_url'] else
                              'citation_is_KRAS_G12D_only' if 'RCV000029214' in r['source_url'] and r['feature']!='G12D' else
                              'allele_specific_review_pending')))
    with (out/'catalogue.tsv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    with (out/'catalogue_audit.tsv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(audit[0]), delimiter='\t')
        w.writeheader()
        w.writerows(audit)
    summary = dict(entries=len(rows), genes=len({r['gene'] for r in rows}),
                   variant_classes=dict(Counter(r['variant_class'] for r in rows)),
                   source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   mapping_status='not_yet_mapped',
                   entries_without_source_url=sum(not r['source_url'] for r in rows),
                   clinical_validation=False,
                   report_columns=['sample_id', 'genesymbol', 'mutated_reads', 'total_reads', 'cdna_mutation', 'protein_mutation'],
                   coverage_rule=dict(hotspot='minimum_including_zeros', broad_inactivation='median_including_zeros',
                                      dual_category='hotspot_minimum_including_zeros'))
    (out/'catalogue_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True)
    p.add_argument('--outdir', required=True)
    a = p.parse_args()
    print(json.dumps(prepare(a.source, a.outdir), indent=2))
