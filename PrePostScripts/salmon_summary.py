#!/usr/bin/env python3
"""Summarize Salmon's raw library-format evidence without choosing a protocol."""
import csv
import json
from pathlib import Path
import sys

FORMATS = ('MSF', 'OSF', 'ISF', 'MSR', 'OSR', 'ISR', 'SF', 'SR', 'MU', 'OU', 'IU', 'U')
METRICS = ('compatible_fragment_ratio', 'num_compatible_fragments',
           'num_incompatible_fragments', 'num_assigned_fragments',
           'num_frags_with_concordant_consistent_mappings',
           'num_frags_with_inconsistent_or_orphan_mappings',
           'strand_mapping_bias')
ADDITIVE_METRICS = tuple(name for name in METRICS if name.startswith('num_'))
FIELDS = ('scope', 'sample', 'expected_format', *METRICS, 'format_mapping_total',
          *(f'{name}_mapping_count' for name in FORMATS),
          *(f'{name}_pct_of_format_mappings' for name in FORMATS))


def load(directory):
    sample = Path(directory).name
    evidence = json.loads((Path(directory) / 'lib_format_counts.json').read_text())
    counts = {name: int(evidence.get(name, 0) or 0) for name in FORMATS}
    return sample, evidence, counts


def as_row(scope, sample, evidence, counts):
    total = sum(counts.values())
    metrics = [evidence.get(name, '') if evidence is not None else '' for name in METRICS]
    values = [counts[name] for name in FORMATS]
    percents = [f'{100 * counts[name] / total:.3f}' if total else '' for name in FORMATS]
    return [scope, sample, evidence.get('expected_format', '') if evidence else '',
            *metrics, total, *values, *percents]


def main(directories):
    evidence = [load(directory) for directory in directories]
    if not evidence:
        raise SystemExit('No Salmon sample outputs were supplied')
    writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    writer.writerow(FIELDS)
    total_counts = {name: 0 for name in FORMATS}
    for sample, metrics, counts in sorted(evidence):
        writer.writerow(as_row('sample', sample, metrics, counts))
        for name in FORMATS:
            total_counts[name] += counts[name]
    totals = {
        name: sum(int(metrics.get(name, 0) or 0) for _, metrics, _ in evidence)
        for name in ADDITIVE_METRICS
    }
    writer.writerow(as_row('batch', 'ALL', totals, total_counts))


if __name__ == '__main__':
    main(sys.argv[1:])
