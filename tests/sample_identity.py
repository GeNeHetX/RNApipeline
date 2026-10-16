#!/usr/bin/env python3
"""Exercise real Nextflow wiring with fake tools that reject mislabeled inputs."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
RUN = Path(tempfile.mkdtemp(prefix='rna-identity-'))
BIN = RUN / 'bin'
BIN.mkdir()
FAKE = r'''#!/usr/bin/env python3
import gzip, json, sys, time
from pathlib import Path
name, a = Path(sys.argv[0]).name, sys.argv[1:]
def opt(flag): return a[a.index(flag)+1]
def marker(path):
    p=Path(path)
    return (gzip.open(p,'rt').read() if p.suffix=='.gz' else p.read_text()).strip().split()[0]
def write(path, text): Path(path).write_text(text+'\n')
if name=='STAR':
    if '--runMode' in a:
        write(Path(opt('--genomeDir'))/'Genome', 'index')
    else:
        read_pos=a.index('--readFilesIn')
        mate1=a[read_pos+1].split(',')
        mate2=[]
        if read_pos+2<len(a) and not a[read_pos+2].startswith('--'):
            mate2=a[read_pos+2].split(',')
        reads=mate1+mate2
        sample=marker(mate1[0]).split('_L')[0].split('_R')[0]
        assert all(marker(x).split('_L')[0].split('_R')[0]==sample for x in reads)
        assert not mate2 or len(mate1)==len(mate2)
        assert not mate2 or all(marker(r1)==marker(r2) for r1,r2 in zip(mate1,mate2))
        assert opt('--outFileNamePrefix')==sample+'StarOut'
        assert opt('--outFilterMultimapNmax')=='10'
        time.sleep((6-int(sample[-1]))*.15)  # Finish in a different order.
        prefix=opt('--outFileNamePrefix')
        write(prefix+'Aligned.sortedByCoord.out.bam',sample)
        write(prefix+'Log.final.out',sample)
elif name=='featureCounts':
    sample=marker(a[-1]); out=opt('-o')
    assert out in [sample+'exonscount.txt',sample+'genecount.txt'], (sample,out)
    write(out,sample); write(out+'.summary',sample)
elif name=='kallisto':
    if a[0]=='index': write(opt('-i'),'index')
    else:
        sample=opt('-o'); reads=[x for x in a if x.endswith('.fastq.gz')]
        markers=[marker(x) for x in reads]
        assert reads and all(x.split('_L')[0].split('_R')[0]==sample for x in markers)
        assert '--single' in a or all(markers[i]==markers[i+1] for i in range(0,len(markers),2))
        Path(sample).mkdir(); write(sample+'/abundance.tsv',sample)
        write(sample+'/run_info.json',json.dumps({'sample':sample}))
elif name=='salmon':
    if a[0]=='index':
        Path(opt('-i')).mkdir()
        write(Path(opt('-i'))/'versionInfo.json','{}')
    else:
        mate1=[]; mate2=[]; mode=None
        for token in a:
            if token=='-1': mode='r1'
            elif token=='-2': mode='r2'
            elif token in ['-p','-o']: mode=None
            elif mode=='r1': mate1.append(token)
            elif mode=='r2': mate2.append(token)
        sample=opt('-o')
        reads=mate1+mate2
        assert mate1 and len(mate1)==len(mate2)
        assert all(marker(x).split('_L')[0].split('_R')[0]==sample for x in reads)
        lane1=[marker(x) for x in mate1]; lane2=[marker(x) for x in mate2]
        assert lane1==lane2
        Path(sample).mkdir()
        isf=80 if sample=='S1' else 800
        isr=20 if sample=='S1' else 200
        write(Path(sample)/'lib_format_counts.json',json.dumps({
            'expected_format':'ISR' if sample=='S1' else 'ISF',
            'compatible_fragment_ratio':0.99,
            'num_compatible_fragments':1000,
            'num_incompatible_fragments':10,
            'num_assigned_fragments':900,
            'num_frags_with_concordant_consistent_mappings':850,
            'num_frags_with_inconsistent_or_orphan_mappings':20,
            'strand_mapping_bias':0.8,
            'ISF':isf,'ISR':isr,'IU':10}))
        write(Path(sample)/'meta_info.json','{}')
        write(Path(sample)/'quant.sf','Name\tLength\tEffectiveLength\tTPM\tNumReads')
elif name=='fastqc':
    for fq in [x for x in a if x.endswith('.fastq.gz')]:
        for ext in ['zip','html']: write(fq[:-9]+'_fastqc.'+ext,marker(fq))
elif name=='multiqc':
    write('multiqc_report.html','QC'); Path('multiqc_data').mkdir(); write('multiqc_data/data.txt','QC')
elif name=='samtools':
    if a[0]=='index': write(a[-1],marker(opt('-b')))
    elif a[0]=='faidx': write(a[1]+'.fai','index')
    else: write(opt('-o'),marker(next(x for x in a if x.endswith('.bam'))))
elif name=='bcftools':
    if a[0]=='mpileup': print(marker(next(x for x in a if x.endswith('.bam'))))
    else: write(opt('-o'),marker(a[1]))
elif name=='vep':
    assert opt('--cache_version')=='107'
    sample=marker(opt('-i')); assert Path(opt('-o')).name.startswith(sample+'_')
    write(opt('-o'),sample)
elif name=='bgzip':
    p=Path(a[-1])
    with gzip.open(str(p)+'.gz','wt') as f: f.write(p.read_text())
    p.unlink()
elif name=='java':
    if 'IndexFeatureFile' in a: write(opt('-I')+'.idx','index')
    elif 'CreateSequenceDictionary' in a: write(next(x[2:] for x in a if x.startswith('O=')),'index')
    else:
        sample=marker(opt('-I'))
        if '-RGSM' in a: assert opt('-RGSM')==sample
        if '-M' in a: write(opt('-M'),sample)
        write(opt('-O'),sample)
elif name=='wget':
    with gzip.open(opt('-O'),'wt') as f: f.write('reference\n')
else: raise AssertionError(name)
'''
(BIN / 'fake').write_text(FAKE)
(BIN / 'fake').chmod(0o755)
for name in ['STAR','featureCounts','kallisto','fastqc','multiqc','samtools',
             'bcftools','vep','bgzip','java','wget','salmon']:
    (BIN / name).symlink_to('fake')
REF = RUN / 'ensembl_v107_GRCh38'
REF.mkdir()
for name in ['ref.gtf','ref.fa','ref.ExonsOnly.bed','kallisto_index']:
    (REF / name).write_text('reference\n')
READS = RUN / 'reads'
READS.mkdir()
SAMPLES = [f'S{i}' for i in range(1,6)]
for sample in SAMPLES:
    for mate in [1,2]:
        import gzip
        with gzip.open(READS / f'{sample}_R{mate}.fastq.gz','wt') as f: f.write(sample+'\n')
SHEET = RUN / 'samples.txt'
SHEET.write_text(''.join(f'{s} {s}_R1.fastq.gz {s}_R2.fastq.gz {READS}\n' for s in SAMPLES))
TRANSCRIPTOME = RUN / 'cdna.fa'
TRANSCRIPTOME.write_text('>transcript\nACGT\n')
LANES = RUN / 'lanes.txt'
lane_lines=[]
for sample in SAMPLES:
    for lane in (['L001','L002'] if sample=='S1' else ['L001']):
        for mate in (1,2):
            import gzip
            with gzip.open(READS / f'{sample}_{lane}_R{mate}.fastq.gz','wt') as f:
                f.write(sample+'_'+lane+'\n')
        lane_lines.append(f'{sample} {sample}_{lane}_R1.fastq.gz {sample}_{lane}_R2.fastq.gz {READS}\n')
LANES.write_text(''.join(lane_lines))
CSV = RUN / 'samples.csv'
CSV.write_text('ID_Sample,suffix1,suffix2\n'+''.join(f'{s},_R1,_R2\n' for s in SAMPLES))
CFG = RUN / 'test.config'
CFG.write_text("""process {
 executor = 'local'
 cpus = 1
 memory = 1.GB
 container = null
 withName: '.*' { container = null }
 withName: 'doSTAR|doOnlySTARnCount' { cpus = 1; memory = 1.GB }
 withName: buildref { cpus = 1; memory = 1.GB }
 withName: FCounts { cpus = 1; memory = 1.GB }
 withName: 'KallistoPE|Kallisto_single_end' { cpus = 1; memory = 1.GB }
 withName: fastqc { cpus = 1; memory = 1.GB }
 withName: multiqc { cpus = 1; memory = 1.GB }
 withName: gatk_vc { cpus = 1; memory = 1.GB }
 withName: samtools_index { cpus = 1; memory = 1.GB }
 withName: samtools_depth { cpus = 1; memory = 1.GB }
 withName: 'Vep.*' { cpus = 1; memory = 1.GB }
 beforeScript = 'export PATH=%s:$PATH'
}
params.max_threads = 1
params.max_star_memory = 1.GB
params.limitBAMsortRAM = 512000000
docker.enabled = false
apptainer.enabled = false
""" % BIN)
NF = os.environ.get('NEXTFLOW', 'nextflow')

def run(case, extra, success=True, reference=REF, sample_sheet=SHEET, pipeline=None):
    cwd=RUN/case; cwd.mkdir()
    pipeline = pipeline or ROOT/'fullPairedEnd.nf'
    cmd=[NF,'run',str(pipeline),'-ansi-log','false','-c',str(CFG),
         '--outputdir',str(cwd/'results'),'--runNumber',case]
    if reference is not None: cmd += ['--ref',str(reference)]
    cmd += [arg.replace('__SHEET__',str(sample_sheet)) for arg in extra]
    result=subprocess.run(cmd,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,
                          env={**os.environ,'NXF_DISABLE_CHECK_LATEST':'true'},timeout=180)
    (cwd/'engine.log').write_text(result.stdout)
    if (result.returncode==0)!=success:
        print(result.stdout); raise AssertionError(f'{case}: unexpected exit {result.returncode}; {cwd}')
    if success and pipeline.name=='fullPairedEnd.nf':
        report=json.loads((cwd/'results'/f'{case}_pipeline_summary.json').read_text())
        assert report['pipeline']=='RNApipeline NF - V1.8.0'
        assert report['version'].startswith('RNAv1.8.0_')
        assert sorted(report['samples'])==SAMPLES
        for sample in SAMPLES:
            if report['run_fcounts']:
                for suffix in ['exonscount.txt','genecount.txt']:
                    assert (cwd/'results/FeatureCounts_output'/f'{sample}{suffix}').read_text().strip()==sample
            if report['run_kallisto']:
                assert (cwd/'results/Kallisto_output'/sample/'abundance.tsv').read_text().strip()==sample
    elif success:
        summary=(cwd/'results/strandedness_summary.tsv').read_text()
        assert 'sample\tS1\tISR' in summary
        assert 'batch\tALL\t' in summary
        assert 'ISF_pct_of_format_mappings' in summary
        assert 'ISR_pct_of_format_mappings' in summary
        assert 'num_compatible_fragments' in summary
        assert all(f'sample\t{sample}\t' in summary for sample in SAMPLES)
    print(f'PASS {case}')
    return result.stdout

try:
    run('paired',['--sampleList','__SHEET__','--single_end','false'])
    run('multi_lane',['--sampleList','__SHEET__','--single_end','false'],sample_sheet=LANES)
    run('derivedref',['--sampleList','__SHEET__','--reference_root',str(RUN)],reference=None)
    run('csv',['--csvSample',str(CSV),'--sampleInputDir',str(READS)])
    single=RUN/'single.txt'
    single.write_text(''.join(f'{s} {s}_R1.fastq.gz - {READS}\n' for s in SAMPLES))
    run('single',['--sampleList',str(single),'--single_end','true'])
    run('variants',['--sampleList','__SHEET__','--mpileup','true',
        '--gatk4','true','--vep','true','--samtools_depth','true','--bed',str(REF/'ref.ExonsOnly.bed')])
    run('builtref',['--sampleList','__SHEET__','--ref','no_ref'])
    invalid=RUN/'invalid.txt'; invalid.write_text(SHEET.read_text().replace('S1 ', '!bad '))
    assert 'Invalid sample ID' in run('invalid',['--sampleList',str(invalid)],False)
    reused=RUN/'reused.txt'; reused.write_text(SHEET.read_text().replace('S2_R1','S1_R1'))
    assert 'assigned more than once' in run('reused',['--sampleList',str(reused)],False)
    missing=RUN/'missing.txt'; missing.write_text(SHEET.read_text().replace('S1_R1','missing_R1'))
    assert 'FASTQ missing' in run('missing',['--sampleList',str(missing)],False)
    run('kallisto_only',['--sampleList','__SHEET__','--star','false','--fcounts','false',
                         '--fastqc','false','--multiqc','false'])
    run('strand_preflight',['--sampleList','__SHEET__','--transcriptome',str(TRANSCRIPTOME)],
        sample_sheet=LANES,pipeline=ROOT/'strandPreflight.nf')
    print('All sample-identity checks passed')
except Exception:
    print(f'Diagnostic files retained: {RUN}')
    raise
else:
    shutil.rmtree(RUN)
