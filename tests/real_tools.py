#!/usr/bin/env python3
"""Exercise the image tools and metrics on synthetic, directly constructed reads/BAM.

Run inside the 2.0 image, with this checkout at /pipeline and a writable cwd.
No real samples, human reference or read mapping is used.
"""
import csv
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys

import pysam

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'PrePostScripts'))
from collect_metrics import collect


def run(*args):
    subprocess.run(args,check=True)


random.seed(10)
sequence=''.join(random.choices('ACGT',k=1000))
Path('ref.fa').write_text('>1\n'+sequence+'\n')
pysam.faidx('ref.fa')
Path('cdna.fa').write_text('>TX.1\n'+sequence+'\n')
Path('ref.gtf').write_text('1\tfixture\texon\t1\t1000\t.\t+\t.\tgene_id "GENE"; exon_id "EXON"; transcript_id "TX";\n')
header=dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN='1',LN=1000)],RG=[dict(ID='S',SM='S')])
with pysam.AlignmentFile('S.bam','wb',header=header) as bam, gzip.open('S.fastq.gz','wt') as fq:
    for i in range(80):
        read=pysam.AlignedSegment()
        read.query_name=f'read{i}'
        bases=list(sequence[100:175])
        if i%2: bases[20]=next(x for x in 'ACGT' if x!=bases[20])
        read.query_sequence=''.join(bases)
        read.flag=0;read.reference_id=0;read.reference_start=100;read.mapping_quality=60
        read.cigar=[(0,75)];read.query_qualities=pysam.qualitystring_to_array('I'*75)
        read.set_tag('RG','S');bam.write(read)
        fq.write(f'@read{i}\n{read.query_sequence}\n+\n'+('I'*75)+'\n')
pysam.index('S.bam')
ref=sequence[120];alt=next(x for x in 'ACGT' if x!=ref)
Path('sites.vcf').write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+f'1\t121\t.\t{ref}\t{alt}\t.\tPASS\t.\n')
run('somalier','extract','--sites','sites.vcf','-f','ref.fa','-d','fingerprints','S.bam')
assert Path('fingerprints/S.somalier').stat().st_size>0
run('fastqc','-q','S.fastq.gz')
for name,extra in [('Sgenecount.txt',[]),('Sexonscount.txt',['-f','-g','exon_id'])]:
    run('featureCounts','-a','ref.gtf','-o',name,*extra,'S.bam')
run('kallisto','index','-i','kallisto_index','cdna.fa')
run('kallisto','quant','--single','-l','75','-s','5','-i','kallisto_index','-o','quant','S.fastq.gz')
Path('STAR.log').write_text('Number of input reads | 80\nUniquely mapped reads % | 100.00%\n')
collect('S',['S.fastq.gz'],['S_fastqc.zip'],'STAR.log',
        ['Sgenecount.txt.summary','Sexonscount.txt.summary'],'Sgenecount.txt',
        'quant/run_info.json',True,'S.metrics.tsv')
with open('S.metrics.tsv') as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert any(r['metric']=='Basic Statistics.Total Sequences' and r['value']=='80' for r in rows)
assert any(r['metric']=='nonzero_genes' and r['value']=='1' for r in rows)
assert any(r['metric']=='n_processed' and r['value']=='80' for r in rows)
assert any(r['metric']=='version' and r['value']=='0.12.1' for r in rows)
run('bcftools','mpileup','-f','ref.fa','-d','100000','-q','0','-Q','13',
    '-a','FORMAT/AD,FORMAT/DP','-Ob','-o','pileup.bcf','S.bam')
run('bcftools','call','pileup.bcf','-m','-v','-P','0.0011','-Ov','-o','calls.vcf')
run('bcftools','norm','-f','ref.fa','-m','-any','-Oz','-o','calls.norm.vcf.gz','calls.vcf')
result=subprocess.run(['bcftools','query','-f','%POS[\t%AD\t%DP]\n','calls.norm.vcf.gz'],
                      check=True,stdout=subprocess.PIPE,text=True)
called=[r.split('\t') for r in result.stdout.splitlines()]
assert any(r==['121','40,40','80'] for r in called),called
Path('smoke_result.json').write_text(json.dumps(dict(status='passed',metrics=len(rows),
    inputs='synthetic directly constructed BAM/FASTQ; no mapping'),indent=2)+'\n')
print('PASS real FastQC, FeatureCounts, Kallisto, Somalier, BCFtools and metrics')
