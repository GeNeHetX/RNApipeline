#!/usr/bin/env python3
"""Build candidate hotspot and coding targets from one Ensembl FASTA/GTF.

Coordinates are reference checked here and independently consequence checked by
VEP during buildReference.nf. Unresolved entries remain in the build audit.
"""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import re

CODONS = dict(zip((a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),
                  'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMP = str.maketrans('ACGTN', 'TGCAN')


def intervals(positions):
    result = []
    for p in sorted(set(positions)):
        if result and result[-1][1] == p-1:
            result[-1][1] = p
        else:
            result.append([p-1, p])
    return result


def build(catalogue, gtf, fasta, outdir, overrides=None):
    import pysam
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(catalogue) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    wanted = {r['gene'] for r in rows}
    transcripts = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip().split('\t')
            if c[2] not in ('transcript', 'CDS'):
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)"', c[8]))
            gene, tid = a.get('gene_name'), a.get('transcript_id')
            if gene not in wanted or not tid:
                continue
            tx = transcripts.setdefault(tid, dict(gene=gene, transcript=tid,
                version=a.get('transcript_version', ''), chrom=c[0], strand=c[6], tags=set(), cds=[]))
            tx['tags'].update(re.findall(r'tag "([^"]+)"', c[8]))
            if c[2] == 'CDS':
                tx['cds'].append((int(c[3]), int(c[4]), int(c[7])))
    by_gene = defaultdict(list)
    for tx in transcripts.values():
        if tx['cds']:
            by_gene[tx['gene']].append(tx)
    selected = {}
    overrides_by_gene = {}
    if overrides:
        with open(overrides) as f:
            overrides_by_gene = {r['gene']: r['transcript_id'] for r in csv.DictReader(f, delimiter='\t')}
    def rank(t):
        return ('MANE_Select' in t['tags'], 'Ensembl_canonical' in t['tags'],
                'appris_principal_1' in t['tags'], 'basic' in t['tags'],
                sum(b-a+1 for a,b,_ in t['cds']))
    fa = pysam.FastaFile(fasta)
    for gene, txs in by_gene.items():
        txs.sort(key=lambda t: (tuple(-int(v) for v in rank(t)), t['transcript']))
        if gene in overrides_by_gene:
            txs = [t for t in txs if overrides_by_gene[gene] in
                   (t['transcript'], t['transcript']+'.'+t['version'])]
            if not txs:
                raise ValueError('Transcript override absent from this GTF: ' + gene)
        tx = txs[0]
        tx['transcript'] += '.'+tx['version'] if tx['version'] else ''
        tx['positions'], seq = [], ''
        try:
            for start, end, phase in sorted(set(tx['cds']), reverse=tx['strand']=='-'):
                if phase != (3-len(seq)%3)%3:
                    raise ValueError('Incomplete/inconsistent CDS phase')
                piece = fa.fetch(tx['chrom'], start-1, end).upper()
                ps = list(range(start, end+1))
                if tx['strand'] == '-':
                    piece, ps = piece.translate(COMP)[::-1], ps[::-1]
                seq += piece
                tx['positions'].extend(ps)
            if not tx['version'] or len(seq)%3 or seq[:3] != 'ATG':
                raise ValueError('Unversioned or incomplete CDS')
            if any(CODONS.get(seq[i:i+3]) in (None, '*') for i in range(0,len(seq),3)):
                raise ValueError('Ambiguous CDS or internal stop')
            tx['sequence'] = seq
        except (ValueError, KeyError) as e:
            tx['error'] = str(e)
        selected[gene] = tx
    broad = {r['gene'] for r in rows if r['variant_class']=='gene_mechanism'} | {'TP53'}
    hotspot_genes = {r['gene'] for r in rows if r['variant_class'] in ('protein_snv','protein_nonsense','residue_unspecified')}
    genes, audit, variants = {}, [], []
    for r in rows:
        gene, kind, feature = r['gene'], r['variant_class'], r['feature']
        tx = selected.get(gene)
        entry = dict(id=r['id'], genesymbol=gene, feature=feature, status='', transcript='', reason='')
        try:
            if kind not in ('protein_snv','protein_nonsense','residue_unspecified','gene_mechanism'):
                raise ValueError('Unsupported event or exact nucleotide allele missing: '+kind)
            if not tx or tx.get('error'):
                raise ValueError(tx.get('error') if tx else 'No coding transcript')
            entry['transcript'] = tx['transcript']
            g = genes.setdefault(gene, dict(transcript=tx['transcript'], chrom=tx['chrom'],
                 hotspot=[], coding=intervals(tx['positions']) if gene in broad else [], codons=[]))
            if kind == 'gene_mechanism':
                entry['status'] = 'coding_target_reference_checked'
            else:
                m = re.fullmatch(r'([A-Z])(\d+)([A-Z*]?)', feature)
                if not m:
                    raise ValueError('Invalid protein target')
                offset = 3*(int(m[2])-1)
                codon = tx['sequence'][offset:offset+3]
                if CODONS.get(codon) != m[1]:
                    raise ValueError('Reference amino acid mismatch on selected transcript')
                positions = tx['positions'][offset:offset+3]
                g['hotspot'].extend(positions)
                g['codons'].append(int(m[2]))
                entry['status'] = 'codon_reference_checked'
                # These are possible one-base alleles, not evidence of clinical occurrence.
                for i in range(3):
                    for alt in 'ACGT':
                        alternate = codon[:i]+alt+codon[i+1:]
                        aa = CODONS[alternate]
                        if alt == codon[i] or (m[3] and aa != m[3]) or (not m[3] and aa == m[1]):
                            continue
                        refbase, altbase = codon[i], alt
                        if tx['strand']=='-':
                            refbase, altbase = refbase.translate(COMP), altbase.translate(COMP)
                        variants.append(dict(id=f"{r['id']}.a{len(variants)+1}", gene=gene,
                            chrom=tx['chrom'], pos=positions[i], ref=refbase, alt=altbase,
                            transcript=tx['transcript'], aa_ref=m[1], aa_alt=aa, aa_pos=int(m[2])))
                if m[3] and not any(v['id'].startswith(r['id']+'.') for v in variants):
                    entry['reason'] = 'No single-base solution; codon retained, exact allele unresolved'
        except ValueError as e:
            entry.update(status='unresolved', reason=str(e))
        audit.append(entry)
    fa.close()
    # A gene with no valid targets is omitted from reporting, retained in audit.
    genes = {g: t for g,t in genes.items() if t['hotspot'] or t['coding']}
    for g in genes.values():
        g['hotspot'] = intervals(g['hotspot'])
        g['codons'] = sorted(set(g['codons']))
    for name,g in genes.items():
        g['category'] = 'dual' if name in broad and name in hotspot_genes else 'broad' if name in broad else 'hotspot'
    (out/'targets.json').write_text(json.dumps(genes, indent=2)+'\n')
    (out/'candidate_alleles.json').write_text(json.dumps(variants, indent=2)+'\n')
    with (out/'mapping_audit.tsv').open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(audit[0]), delimiter='\t'); w.writeheader(); w.writerows(audit)
    regions = defaultdict(list)
    for g in genes.values():
        regions[g['chrom']].extend(g['hotspot']+g['coding'])
    with (out/'targets.bed').open('w') as f:
        for chrom, spans in sorted(regions.items()):
            merged=[]
            for start,end in sorted(spans):
                if merged and start <= merged[-1][1]: merged[-1][1]=max(merged[-1][1],end)
                else: merged.append([start,end])
            # Include a left anchor for indels; coverage statistics retain exact targets.
            for start,end in merged: f.write(f'{chrom}\t{max(0,start-1)}\t{end}\n')
    with (out/'candidate_variants.vcf').open('w') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for v in sorted(variants, key=lambda v:(v['chrom'],v['pos'],v['id'])):
            f.write(f"{v['chrom']}\t{v['pos']}\t{v['id']}\t{v['ref']}\t{v['alt']}\t.\t.\t.\n")
    if not genes or not variants:
        raise ValueError('No valid hotspot targets were built; inspect mapping_audit.tsv')


def validate(alleles, vep, output):
    expected={v['id']:v for v in json.loads(Path(alleles).read_text())}
    seen=set()
    with open(vep) as f:
        for line in f:
            data=json.loads(line)
            original=data['input'].split('\t')
            key=original[2]
            v=expected[key]
            if (original[0],int(original[1]),original[3],original[4]) != (v['chrom'],v['pos'],v['ref'],v['alt']):
                raise ValueError('VEP input allele mismatch: '+key)
            if key in seen: raise ValueError('Duplicate VEP allele: '+key)
            seen.add(key)
            matches=[t for t in data.get('transcript_consequences', [])
                     if t.get('transcript_id','').split('.')[0]==v['transcript'].split('.')[0]
                     and str(t.get('protein_start'))==str(v['aa_pos'])
                     and t.get('amino_acids')==v['aa_ref']+'/'+v['aa_alt']]
            if not matches: raise ValueError('VEP consequence mismatch: '+key)
    if seen != set(expected): raise ValueError('Missing VEP allele validation')
    Path(output).write_text(json.dumps(dict(candidate_alleles=len(seen), status='passed',
        clinical_validation=False), indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command', required=True)
    b=sub.add_parser('build')
    for field in ('catalogue','gtf','fasta','outdir'): b.add_argument('--'+field, required=True)
    b.add_argument('--transcript-overrides')
    v=sub.add_parser('validate')
    for field in ('alleles','vep','output'): v.add_argument('--'+field, required=True)
    a=p.parse_args()
    if a.command=='build': build(a.catalogue,a.gtf,a.fasta,a.outdir,a.transcript_overrides)
    else: validate(a.alleles,a.vep,a.output)
