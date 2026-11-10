#!/usr/bin/env python3
"""Write concatenatable per-sample QC metrics; never copy FastQC distributions."""
import argparse
import csv
import json
import re
import zipfile
from pathlib import Path

FIELDS = ['sample_id', 'metric_scope', 'fastq_id', 'read_role', 'metric_source',
          'metric', 'value', 'unit', 'qc_status']
KALLISTO_FIELDS = ('n_targets', 'n_bootstraps', 'n_processed', 'n_pseudoaligned',
                   'n_unique', 'p_pseudoaligned', 'p_unique', 'kallisto_version',
                   'index_version')


def collect(sample, reads, fastqc, star, summaries, genes, kallisto, single_end, output):
    rows = []

    def add(source, metric, value, unit='', status='', fastq_id='', role=''):
        rows.append([sample, 'fastq' if fastq_id else 'sample', fastq_id, role,
                     source, metric, value, unit, status])

    # File order comes from the sample sheet, including all lanes.
    for i, read in enumerate(reads):
        name = Path(read).name
        stem = re.sub(r'\.(fastq|fq)(\.gz)?$', '', name)
        archive = next(p for p in fastqc if Path(p).name == stem + '_fastqc.zip')
        role = 'R1' if single_end or i % 2 == 0 else 'R2'
        with zipfile.ZipFile(archive) as z:
            data = z.read(next(n for n in z.namelist() if n.endswith('/fastqc_data.txt'))).decode()
        modules = {}
        version = re.search(r'^##FastQC\s+([^\n]+)', data)
        if version:
            add('FastQC', 'version', version[1], fastq_id=name, role=role)
        for block in data.split('>>END_MODULE'):
            lines = block.strip().splitlines()
            start = next((i for i,line in enumerate(lines) if line.startswith('>>')), None)
            if start is None:
                continue
            lines = lines[start:]
            module, status = lines[0][2:].split('\t')
            values = [line.split('\t') for line in lines[1:] if line and not line.startswith('#')]
            modules[module] = (status, values, lines[1:])
            add('FastQC', module + '.status', status, status=status, fastq_id=name, role=role)
        def metric(module, key, value, unit=''):
            add('FastQC', module + '.' + key, value, unit, modules[module][0], name, role)
        for key, value, *_ in modules['Basic Statistics'][1]:
            if key not in ('Filename', 'File type'):
                metric('Basic Statistics', key, value, 'percent' if key == '%GC' else '')
        # Scalar summaries only. No per-base, per-quality, GC or length histogram rows.
        for module, column, label, reducer in [
            ('Per base sequence quality', 1, 'mean_quality_min', min),
            ('Per base sequence quality', 1, 'mean_of_reported_base_bin_means', lambda v: sum(v)/len(v)),
            ('Per base N content', 1, 'max_N_percent', max),
            ('Adapter Content', None, 'max_adapter_percent', max),
            ('Kmer Content', 3, 'max_observed_expected', max),
            ('Per tile sequence quality', 2, 'max_absolute_quality_deviation', lambda v: max(abs(x) for x in v)),
        ]:
            if module in modules:
                values = [float(x) for row in modules[module][1]
                          for x in (row[1:] if column is None else row[column:column+1])]
                if values:
                    metric(module, label, reducer(values))
        if 'Per base sequence content' in modules:
            values = modules['Per base sequence content'][1]
            if values:
                metric('Per base sequence content', 'max_base_bias_from_25_percent',
                       max(abs(float(x)-25) for r in values for x in r[1:]), 'percent')
        if 'Per sequence GC content' in modules:
            values = modules['Per sequence GC content'][1]
            total = sum(float(r[1]) for r in values)
            if total:
                metric('Per sequence GC content', 'mean_GC_percent',
                       sum(float(r[0])*float(r[1]) for r in values)/total, 'percent')
        if 'Kmer Content' in modules:
            metric('Kmer Content', 'n_kmers', len(modules['Kmer Content'][1]), 'count')
        if 'Sequence Duplication Levels' in modules:
            for line in modules['Sequence Duplication Levels'][2]:
                if line.startswith('#Total Deduplicated Percentage'):
                    metric('Sequence Duplication Levels', 'deduplicated_percent', line.split('\t')[1], 'percent')
        if 'Overrepresented sequences' in modules:
            values = modules['Overrepresented sequences'][1]
            metric('Overrepresented sequences', 'n_sequences', len(values), 'count')
            metric('Overrepresented sequences', 'total_percent', sum(float(r[2]) for r in values), 'percent')
        if 'Per sequence quality scores' in modules:
            values = modules['Per sequence quality scores'][1]
            n = sum(float(r[1]) for r in values)
            if n:
                metric('Per sequence quality scores', 'mean_quality', sum(float(r[0])*float(r[1]) for r in values)/n)
        if 'Sequence Length Distribution' in modules:
            values = modules['Sequence Length Distribution'][1]
            if values:
                metric('Sequence Length Distribution', 'minimum_length', min(int(r[0].split('-')[0]) for r in values), 'bases')
                metric('Sequence Length Distribution', 'maximum_length', max(int(r[0].split('-')[-1]) for r in values), 'bases')

    for line in Path(star).read_text().splitlines():
        if '|' in line:
            key, value = (p.strip() for p in line.split('|', 1))
            if key and value:
                add('STAR', key, value.rstrip('%'), 'percent' if value.endswith('%') else '')
    for summary in summaries:
        scope = 'exon' if 'exonscount' in str(summary) else 'gene'
        with open(summary) as f:
            table = csv.reader(f, delimiter='\t')
            next(table)
            counts = [(r[0], int(r[1])) for r in table]
        total = sum(v for _, v in counts)
        for category, count in counts:
            add('FeatureCounts', scope + '.' + category, count, 'count')
            add('FeatureCounts', scope + '.' + category + '_percent', 100*count/total if total else 0, 'percent')
    with open(genes) as f:
        nonzero = sum(float(r.rstrip().split('\t')[-1]) > 0 for r in f
                      if r.strip() and not r.startswith(('#', 'Geneid\t')))
    add('FeatureCounts', 'nonzero_genes', nonzero, 'count')
    add('Kallisto', 'status', 'enabled' if kallisto else 'disabled', status='enabled' if kallisto else 'disabled')
    if kallisto:
        info = json.loads(Path(kallisto).read_text())
        for key in KALLISTO_FIELDS:
            if key in info:
                add('Kallisto', key, info[key], 'percent' if key.startswith('p_') else '')
    with open(output, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(FIELDS)
        w.writerows(rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--reads', nargs='+', required=True)
    p.add_argument('--fastqc', nargs='+', required=True)
    p.add_argument('--star', required=True)
    p.add_argument('--summaries', nargs='+', required=True)
    p.add_argument('--genes', required=True)
    p.add_argument('--kallisto')
    p.add_argument('--single-end', action='store_true')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    collect(a.sample, a.reads, a.fastqc, a.star, a.summaries, a.genes,
            a.kallisto, a.single_end, a.output)
