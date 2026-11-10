#!/usr/bin/env python3
"""Reference contig checks, Somalier panel preparation and build manifest."""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()


def assemble(reference, validation, root):
    # Stage links, then let publishDir copy the validated bundle once. Avoid an
    # extra full copy of the human STAR index and VEP cache in the work area.
    def link(source, destination):
        Path(destination).symlink_to(Path(source).resolve())
        return destination
    shutil.copytree(reference, root, copy_function=link)
    link(validation, Path(root)/'hotspots'/'validation.json')


def prepare(root, panel, release):
    import pysam
    root=Path(root)
    fa=pysam.FastaFile(str(root/'ref.fa'))
    lengths=dict(zip(fa.references,fa.lengths))
    exons=defaultdict(list)
    header=[]
    transcript_versions={}
    with (root/'ref.gtf').open() as f:
        for line in f:
            if line.startswith('#'):
                header.append(line.strip())
                continue
            c=line.rstrip().split('\t')
            if len(c)!=9: raise ValueError('Invalid GTF line')
            if c[0] not in lengths or not 1<=int(c[3])<=int(c[4])<=lengths[c[0]]:
                raise ValueError('GTF interval incompatible with FASTA: '+line[:150])
            if c[2]=='exon': exons[c[0]].append((int(c[3])-1,int(c[4])))
            if c[2]=='transcript':
                attrs=dict(re.findall(r'(\w+) "([^"]+)"',c[8]))
                transcript_versions[attrs['transcript_id']]=attrs.get('transcript_version','')
    if not any(line in ('#!genome-build GRCh38','#!genome-version GRCh38') for line in header):
        raise ValueError('This builder requires the human GRCh38 GTF')
    if not exons: raise ValueError('GTF has no exons')
    shared=0
    with (root/'cdna.fa').open() as f:
        for line in f:
            if not line.startswith('>'): continue
            tid,_,version=line[1:].split()[0].partition('.')
            if tid in transcript_versions:
                if version!=transcript_versions[tid]:
                    raise ValueError('cDNA/GTF transcript version mismatch: '+tid)
                shared+=1
    if not shared: raise ValueError('cDNA and GTF have no versioned transcripts in common')
    with (root/'ref.ExonsOnly.bed').open('w') as f:
        for chrom in sorted(exons):
            merged=[]
            for start,end in sorted(exons[chrom]):
                if merged and start<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],end)
                else: merged.append([start,end])
            for start,end in merged: f.write(f'{chrom}\t{start}\t{end}\n')
    cache=root/'VEP'/'homo_sapiens'/f'{release}_GRCh38'/'info.txt'
    if not cache.is_file(): raise ValueError('Expected same-release GRCh38 VEP cache: '+str(cache))
    (root/'somalier').mkdir(exist_ok=True)
    n=0
    with gzip.open(panel,'rt') as source,(root/'somalier'/'sites.vcf').open('w') as target:
        for line in source:
            if line.startswith('##contig='): continue
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    for chrom,length in lengths.items(): target.write(f'##contig=<ID={chrom},length={length}>\n')
                target.write(line)
                continue
            c=line.rstrip().split('\t')
            aliases=[c[0],c[0].removeprefix('chr'),'chr'+c[0]]
            chrom=next((a for a in aliases if a in lengths),None)
            if not chrom: raise ValueError('Somalier contig absent in FASTA: '+c[0])
            pos=int(c[1]);ref=c[3].upper()
            if fa.fetch(chrom,pos-1,pos-1+len(ref)).upper()!=ref:
                raise ValueError(f'Somalier REF mismatch at {chrom}:{pos}')
            c[0]=chrom;target.write('\t'.join(c)+'\n');n+=1
    fa.close()
    if not n: raise ValueError('Empty Somalier panel')
    (root/'somalier'/'panel.json').write_text(json.dumps(dict(sites=n,
        original_sha256=sha256(panel), contigs='matched_to_reference', assembly='GRCh38'),indent=2)+'\n')


def manifest(root, release, sources, pipeline_version, containers):
    root=Path(root)
    required=['ref.fa','ref.fa.fai','ref.gtf','cdna.fa','kallisto_index','Genome',
              'SA','SAindex','genomeParameters.txt','ref.ExonsOnly.bed',
              'somalier/sites.vcf.gz','somalier/sites.vcf.gz.tbi',
              'hotspots/catalogue.tsv','hotspots/targets.bed','hotspots/targets.json',
              'hotspots/mapping_audit.tsv','hotspots/validation.json',
              f'VEP/homo_sapiens/{release}_GRCh38/info.txt']
    hashes={}
    for name in required:
        path=root/name
        if not path.is_file() or not path.stat().st_size: raise ValueError('Missing reference asset: '+name)
        hashes[name]=sha256(path)
    validation=json.loads((root/'hotspots/validation.json').read_text())
    if validation['status']!='passed': raise ValueError('Hotspot VEP validation did not pass')
    versions={}
    for tool,args in [('STAR',['--version']),('kallisto',['version']),('samtools',['--version']),
                      ('bcftools',['--version']),('somalier',[])]:
        result=subprocess.run([tool]+args,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=True)
        versions[tool]=result.stdout.splitlines()[0]
    data=dict(pipeline_version=pipeline_version,species='homo_sapiens',assembly='GRCh38',
              ensembl_release=int(release),vep_cache_release=int(release),
              sources=json.loads(Path(sources).read_text()),sha256=hashes,tools=versions,
              containers=json.loads(containers),star=dict(sjdbOverhang=100,genomeSAindexNbases=11),
              hotspot_validation=validation,
              hotspot_coverage_rule='minimum for hotspot and dual-category; median for broad-only; includes zero bases')
    (root/'reference_manifest.json').write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['prepare','manifest','assemble'])
    p.add_argument('--root',required=True)
    p.add_argument('--release',required=True)
    p.add_argument('--panel')
    p.add_argument('--sources')
    p.add_argument('--pipeline-version')
    p.add_argument('--containers')
    p.add_argument('--reference')
    p.add_argument('--validation')
    a=p.parse_args()
    if a.command=='prepare': prepare(a.root,a.panel,a.release)
    elif a.command=='assemble': assemble(a.reference,a.validation,a.root)
    else: manifest(a.root,a.release,a.sources,a.pipeline_version,a.containers)
