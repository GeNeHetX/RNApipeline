#!/usr/bin/env python3
"""Six-column RNA hotspot report using the reference builder's selected transcripts."""
import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import re

FIELDS = ['sample_id','genesymbol','mutated_reads','total_reads','cdna_mutation','protein_mutation']
PROTEIN_ALTERING = {'missense_variant','stop_gained','stop_lost','start_lost',
                    'frameshift_variant','inframe_insertion','inframe_deletion','protein_altering_variant'}


def contains(spans, position):
    i = bisect_left(spans, [position]) - 1
    return i >= 0 and position <= spans[i][1]


def overlaps(spans, position, ref, alt):
    end = position + max(len(ref), 2 if len(alt)>len(ref) else 1) - 1
    i = bisect_left(spans, [end]) - 1
    return i >= 0 and spans[i][1] >= position


def histogram_median(hist):
    n = sum(hist.values())
    indices = ((n-1)//2, n//2)
    values=[]
    cumulative=0
    for depth, count in sorted(hist.items()):
        values.extend(depth for i in indices[len(values):] if cumulative <= i < cumulative+count)
        cumulative += count
    return sum(values)/2


def report(sample, targets, depths, vcf, output):
    genes=json.loads(Path(targets).read_text())
    by_chrom=defaultdict(list)
    hist={}
    for gene,t in genes.items():
        by_chrom[t['chrom']].append((gene,t))
        hist[gene]={kind:Counter({0:sum(e-s for s,e in t[kind])}) for kind in ('hotspot','coding')}
    with open(depths) as f:
        for line in f:
            chrom,pos,depth=line.rstrip().split('\t')
            pos,depth=int(pos),0 if depth=='.' else int(depth)
            for gene,t in by_chrom[chrom]:
                for kind in ('hotspot','coding'):
                    if contains(t[kind],pos):
                        hist[gene][kind][0]-=1
                        hist[gene][kind][depth]+=1
    rows,detected=[],set()
    opener=gzip.open if str(vcf).endswith('.gz') else open
    with opener(vcf,'rt') as f:
        for line in f:
            if not line.strip(): continue
            if line.startswith('##INFO=<ID=CSQ,'):
                csq_fields=re.search(r'Format: ([^"]+)',line)[1].strip().split('|')
            if line.startswith('#'): continue
            c=line.rstrip().split('\t')
            if c[6] not in ('.','PASS'): continue
            chrom,pos,ref,alt=c[0],int(c[1]),c[3],c[4]
            if ',' in alt: raise ValueError('Report requires normalized biallelic VCF records')
            fmt=dict(zip(c[8].split(':'),c[9].split(':')))
            ad=[int(x) for x in fmt['AD'].split(',')]
            depth=int(fmt['DP']) if fmt.get('DP','.')!='.' else sum(ad)
            if len(ad)!=2 or ad[1]==0 or '1' not in re.split(r'[/|]',fmt.get('GT','')): continue
            info=dict(x.split('=',1) for x in c[7].split(';') if '=' in x)
            consequences=[dict(zip(csq_fields,x.split('|'))) for x in info.get('CSQ','').split(',')]
            for gene,t in by_chrom[chrom]:
                in_hotspot = overlaps(t['hotspot'],pos,ref,alt)
                in_coding = overlaps(t['coding'],pos,ref,alt)
                if not (in_hotspot or in_coding): continue
                for consequence in consequences:
                    if consequence.get('Feature','').split('.')[0]!=t['transcript'].split('.')[0]: continue
                    terms=set(consequence.get('Consequence','').split('&'))
                    hotspot=in_hotspot and bool(terms & PROTEIN_ALTERING)
                    broad=in_coding and bool(terms & {'frameshift_variant','stop_gained'})
                    if not (hotspot or broad): continue
                    rows.append([sample,gene,ad[1],depth,consequence.get('HGVSc',''),consequence.get('HGVSp','')])
                    detected.add(gene)
                    break
    for gene,t in sorted(genes.items()):
        if gene in detected: continue
        if t.get('category') in ('dual','hotspot') and not t['hotspot']:
            continue  # Unresolved hotspot coverage cannot produce a no-call row.
        # For dual-category genes the selected rule is HOTSPOT MINIMUM.
        if t['hotspot']:
            depth=min(d for d,n in hist[gene]['hotspot'].items() if n>0)
        else:
            depth=histogram_median(hist[gene]['coding'])
        rows.append([sample,gene,0,depth,'',''])
    with open(output,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(FIELDS);w.writerows(sorted(rows,key=lambda r:(r[1],r[4],r[5])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for field in ('sample','targets','depths','vcf','output'): p.add_argument('--'+field,required=True)
    a=p.parse_args()
    report(a.sample,a.targets,a.depths,a.vcf,a.output)
