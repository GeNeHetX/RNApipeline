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
        reads=[x for x in a if x.endswith('.fastq.gz')]
        sample=marker(reads[0]); assert all(marker(x)==sample for x in reads)
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
        assert reads and all(marker(x)==sample for x in reads)
        Path(sample).mkdir(); write(sample+'/abundance.tsv',sample)
        write(sample+'/run_info.json',json.dumps({'sample':sample}))
elif name=='fastqc':
    for fq in [x for x in a if x.endswith('.fastq.gz')]:
        for ext in ['zip','html']: write(fq[:-9]+'_fastqc.'+ext,marker(fq))
elif name=='multiqc':
    write('multiqc_report.html','QC'); Path('multiqc_data').mkdir(); write('multiqc_data/data.txt','QC')
elif name=='samtools':
    if a[0]=='index': write(opt('-o'),marker(opt('-b')))
    elif a[0]=='faidx': write(a[1]+'.fai','index')
    else: write(opt('-o'),marker(next(x for x in a if x.endswith('.bam'))))
elif name=='mosdepth':
    prefix,bam=a[-2:]; sample=marker(bam)
    assert prefix==sample+'_coverage' and marker(bam+'.bai')==sample
    with gzip.open(prefix+'.per-base.bed.gz','wt') as f: f.write(sample+'\t0\t100\t5\n')
elif name=='bedtools': print(sys.stdin.read(),end='')
elif name=='run_deepvariant':
    opts=dict(x[2:].split('=',1) for x in a if x.startswith('--') and '=' in x)
    sample=marker(opts['reads']); assert marker(opts['regions'])==sample
    assert Path(opts['output_vcf']).name==sample+'.output.vcf.gz'
    with gzip.open(opts['output_vcf'],'wt') as f: f.write(sample+'\n')
elif name=='bcftools':
    if a[0]=='mpileup': print(marker(next(x for x in a if x.endswith('.bam'))))
    else: write(opt('-o'),marker(a[1]))
elif name=='vep':
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
for name in ['STAR','featureCounts','kallisto','fastqc','multiqc','samtools','mosdepth',
             'bedtools','run_deepvariant','bcftools','vep','bgzip','java','wget']:
    (BIN / name).symlink_to('fake')
REF = RUN / 'ensembl_v107_GRCh38'
REF.mkdir()
for name in ['ref.gtf','ref.fa','ref.ExonsOnly.bed','kallisto_index']:
    (REF / name).write_text('reference\n')
MODEL = RUN / 'model'
MODEL.mkdir()
(MODEL / 'model.ckpt').write_text('model')
READS = RUN / 'reads'
READS.mkdir()
SAMPLES = [f'S{i}' for i in range(1,6)]
for sample in SAMPLES:
    for mate in [1,2]:
        import gzip
        with gzip.open(READS / f'{sample}_R{mate}.fastq.gz','wt') as f: f.write(sample+'\n')
SHEET = RUN / 'samples.txt'
SHEET.write_text(''.join(f'{s} {s}_R1.fastq.gz {s}_R2.fastq.gz {READS}\n' for s in SAMPLES))
CSV = RUN / 'samples.csv'
CSV.write_text('ID_Sample,suffix1,suffix2\n'+''.join(f'{s},_R1,_R2\n' for s in SAMPLES))
CFG = RUN / 'test.config'
CFG.write_text("""process {
 executor = 'local'
 cpus = 1
 memory = 1.GB
 container = null
 withName: '.*' { container = null }
 beforeScript = 'export PATH=%s:$PATH'
}
docker.enabled = false
apptainer.enabled = false
""" % BIN)
NF = os.environ.get('NEXTFLOW', 'nextflow')

def run(case, extra, success=True, reference=REF):
    cwd=RUN/case; cwd.mkdir()
    cmd=[NF,'run',str(ROOT/'fullPairedEnd.nf'),'-ansi-log','false','-c',str(CFG),
         '--outputdir',str(cwd/'results'),'--runNumber',case]
    if reference is not None: cmd += ['--ref',str(reference)]
    cmd += extra
    result=subprocess.run(cmd,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,
                          env={**os.environ,'NXF_DISABLE_CHECK_LATEST':'true'},timeout=180)
    (cwd/'engine.log').write_text(result.stdout)
    if (result.returncode==0)!=success:
        print(result.stdout); raise AssertionError(f'{case}: unexpected exit {result.returncode}; {cwd}')
    if success:
        report=json.loads((cwd/'results'/f'{case}_pipeline_summary.json').read_text())
        assert sorted(report['samples'])==SAMPLES
        for sample in SAMPLES:
            if report['run_fcounts']:
                for suffix in ['exonscount.txt','genecount.txt']:
                    assert (cwd/'results/FeatureCounts_output'/f'{sample}{suffix}').read_text().strip()==sample
            if report['run_kallisto']:
                assert (cwd/'results/Kallisto_output'/sample/'abundance.tsv').read_text().strip()==sample
    print(f'PASS {case}')
    return result.stdout

try:
    run('paired',['--sampleList',str(SHEET),'--single_end','false'])
    run('derivedref',['--sampleList',str(SHEET),'--reference_root',str(RUN)],reference=None)
    run('csv',['--csvSample',str(CSV),'--sampleInputDir',str(READS)])
    single=RUN/'single.txt'
    single.write_text(''.join(f'{s} {s}_R1.fastq.gz - {READS}\n' for s in SAMPLES))
    run('single',['--sampleList',str(single),'--single_end','true'])
    run('variants',['--sampleList',str(SHEET),'--deepvariant','true','--mpileup','true',
        '--gatk4','true','--vep','true','--samtools_depth','true','--bed',str(REF/'ref.ExonsOnly.bed'),
        '--modelckptdeepar',str(MODEL)])
    run('builtref',['--sampleList',str(SHEET),'--ref','no_ref'])
    duplicate=RUN/'duplicate.txt'; duplicate.write_text(SHEET.read_text()+SHEET.read_text().splitlines()[0]+'\n')
    assert 'Duplicate sample IDs' in run('duplicate',['--sampleList',str(duplicate)],False)
    reused=RUN/'reused.txt'; reused.write_text(SHEET.read_text().replace('S2_R1','S1_R1'))
    assert 'assigned more than once' in run('reused',['--sampleList',str(reused)],False)
    missing=RUN/'missing.txt'; missing.write_text(SHEET.read_text().replace('S1_R1','missing_R1'))
    assert 'FASTQ missing' in run('missing',['--sampleList',str(missing)],False)
    run('kallisto_only',['--sampleList',str(SHEET),'--star','false','--fcounts','false',
                         '--fastqc','false','--multiqc','false'])
    print('All sample-identity checks passed')
except Exception:
    print(f'Diagnostic files retained: {RUN}')
    raise
else:
    shutil.rmtree(RUN)
