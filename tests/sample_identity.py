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
        write(Path(opt('--genomeDir'))/'genomeParameters.txt', 'parameters')
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
    assert '-O' in a
    assert ('-p' in a)==('--countReadPairs' in a)
    assert out in [sample+'exonscount.txt',sample+'genecount.txt'], (sample,out)
    write(out,sample); write(out+'.summary',sample)
elif name=='kallisto':
    if a[0]=='index': write(opt('-i'),'index')
    else:
        sample=opt('-o'); reads=[x for x in a if x.endswith('.fastq.gz')]
        markers=[marker(x) for x in reads]
        assert reads and all(x.split('_L')[0].split('_R')[0]==sample for x in markers)
        if '--single' in a:
            assert len(reads)==1
            with gzip.open(reads[0], 'rt') as f:
                tokens=f.read().split()
                assert tokens and all(token.split('_L')[0]==sample for token in tokens)
        else:
            assert len(reads)==2 and reads[0].endswith('_R1.fastq.gz') and reads[1].endswith('_R2.fastq.gz')
            with gzip.open(reads[0], 'rt') as r1, gzip.open(reads[1], 'rt') as r2:
                assert r1.read().split()==r2.read().split()
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
elif name=='Rscript':
    import gzip
    for path in ['StarFCGeneCount_s5.tsv.gz','StarFCExonCount_s5.tsv.gz','StarMetrics_s5.tsv']:
        if path.endswith('.gz'):
            with gzip.open(path,'wt') as f: f.write('matrix\n')
        else: write(path,'metrics')
elif name=='samtools':
    if a[0]=='index': write(a[-1],marker(opt('-b')))
    elif a[0]=='faidx': write(a[1]+'.fai','index')
    else: write(opt('-o'),marker(next(x for x in a if x.endswith('.bam'))))
elif name=='bcftools':
    if a[0]=='mpileup': print(marker(next(x for x in a if x.endswith('.bam'))))
    else:
        assert '-m' in a and '-v' in a and '-P' in a and '-Ov' in a
        write(opt('-o'),marker(a[1]))
elif name=='vep':
    assert opt('--cache_version')=='113'
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
    import tarfile, tempfile
    output=Path(opt('-O'))
    if output.name=='vep_cache.tar.gz':
        with tempfile.TemporaryDirectory() as tmp:
            cache=Path(tmp)/'homo_sapiens'/'113_GRCh38'
            cache.mkdir(parents=True)
            write(cache/'info.txt','cache')
            with tarfile.open(output,'w:gz') as archive:
                archive.add(Path(tmp)/'homo_sapiens',arcname='homo_sapiens')
    else:
        content='reference\n'
        if output.name=='ref.gtf.gz':
            content='1\ttest\texon\t2\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        with gzip.open(output,'wt') as f: f.write(content)
else: raise AssertionError(name)
'''
(BIN / 'fake').write_text(FAKE)
(BIN / 'fake').chmod(0o755)
for name in ['STAR','featureCounts','kallisto','fastqc','multiqc','Rscript','samtools',
             'bcftools','vep','bgzip','java','wget','salmon']:
    (BIN / name).symlink_to('fake')
REF = RUN / 'ensembl_v107_GRCh38'
REF.mkdir()
for name in ['ref.gtf','ref.fa','ref.ExonsOnly.bed','kallisto_index']:
    (REF / name).write_text('reference\n')
(REF / 'VEP' / 'homo_sapiens' / '113_GRCh38').mkdir(parents=True)
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
SINGLE_LANES = RUN/'single_lanes.txt'
single_lane_lines=[]
for sample in SAMPLES:
    for lane in (['L001','L002'] if sample=='S1' else ['L001']):
        import gzip
        with gzip.open(READS/f'{sample}_{lane}_R1.fastq.gz','wt') as f: f.write(sample+'_'+lane+'\n')
        single_lane_lines.append(f'{sample} {sample}_{lane}_R1.fastq.gz - {READS}\n')
SINGLE_LANES.write_text(''.join(single_lane_lines))
CSV = RUN / 'samples.csv'
CSV.write_text('ID_Sample,suffix1,suffix2\n'+''.join(f'{s},_R1,_R2\n' for s in SAMPLES))
CFG = RUN / 'test.config'
CFG.write_text("""process {
 executor = 'local'
 cpus = 1
 memory = 1.GB
 container = null
 withName: '.*' { container = null }
 withName: doSTAR { cpus = 1; memory = 1.GB }
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
    pipeline = pipeline or ROOT/'mainRNAseq.nf'
    cmd=[NF,'run',str(pipeline),'-ansi-log','false','-c',str(CFG),
         '--outputdir',str(cwd/'results'),'--runNumber',case]
    if reference is not None: cmd += ['--ref',str(reference)]
    cmd += [arg.replace('__SHEET__',str(sample_sheet)) for arg in extra]
    result=subprocess.run(cmd,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,
                          env={**os.environ,'NXF_DISABLE_CHECK_LATEST':'true'},timeout=180)
    (cwd/'engine.log').write_text(result.stdout)
    if (result.returncode==0)!=success:
        print(result.stdout); raise AssertionError(f'{case}: unexpected exit {result.returncode}; {cwd}')
    if success and pipeline.name=='mainRNAseq.nf':
        report=json.loads((cwd/'results'/f'{case}_pipeline_summary.json').read_text())
        assert report['pipeline']=='RNApipeline NF - V1.9.0'
        assert report['version'].startswith('RNAv1.9.0_')
        assert sorted(report['samples'])==SAMPLES
        assert report['run_star'] and report['run_fastqc'] and report['run_fcounts'] and report['run_multiqc']
        for sample in SAMPLES:
            for suffix in ['exonscount.txt','genecount.txt']:
                assert (cwd/'results/FeatureCounts_output'/f'{sample}{suffix}').read_text().strip()==sample
            if report['run_kallisto']:
                assert (cwd/'results/Kallisto_output'/sample/'abundance.tsv').read_text().strip()==sample
            else:
                assert not (cwd/'results/Kallisto_output'/sample).exists()
        if report['run_aggregate']:
            assert (cwd/'results/Aggregated_output/StarFCGeneCount_s5.tsv.gz').is_file()
        if report['run_variant_calling']:
            assert report['run_mpileup'] and report['run_vep']
            assert all((cwd/'results/Vep_output'/f'{sample}_mpileup_annot.vcf.gz').is_file() for sample in SAMPLES)
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

def process_scripts(case, command):
    return [p.read_text() for p in (RUN/case/'work').rglob('.command.sh')
            if command in p.read_text()]

try:
    run('paired',['--sampleList','__SHEET__','--single_end','false'])
    paired_counts=process_scripts('paired','featureCounts -T')
    assert paired_counts and all('-p --countReadPairs' in s and ' -O ' in s for s in paired_counts)
    run('multi_lane',['--sampleList','__SHEET__','--single_end','false'],sample_sheet=LANES)
    run('derivedref',['--sampleList','__SHEET__','--reference_root',str(RUN)],reference=None)
    run('csv',['--csvSample',str(CSV),'--sampleInputDir',str(READS)])
    single=RUN/'single.txt'
    single.write_text(''.join(f'{s} {s}_R1.fastq.gz - {READS}\n' for s in SAMPLES))
    run('single',['--sampleList',str(single),'--single_end','true'])
    single_counts=process_scripts('single','featureCounts -T')
    assert single_counts and all(' -p ' not in s and ' -O ' in s for s in single_counts)
    run('single_multi_lane',['--sampleList',str(SINGLE_LANES),'--single_end','true'])
    run('variants',['--sampleList','__SHEET__','--variant_calling','true'])
    run('aggregate',['--sampleList','__SHEET__','--aggregate','true'])
    run('no_kallisto',['--sampleList','__SHEET__','--kallisto','false'])
    assert '`--ref no_ref` is not supported' in run('builtref',['--sampleList','__SHEET__','--ref','no_ref'],False)
    assert 'always on' in run('disable_core',['--sampleList','__SHEET__','--star','false'],False)
    assert 'params.variant_calling' in run('legacy_variant_switch',['--sampleList','__SHEET__','--mpileup','true'],False)
    invalid=RUN/'invalid.txt'; invalid.write_text(SHEET.read_text().replace('S1 ', '!bad '))
    assert 'Invalid sample ID' in run('invalid',['--sampleList',str(invalid)],False)
    reused=RUN/'reused.txt'; reused.write_text(SHEET.read_text().replace('S2_R1','S1_R1'))
    assert 'assigned more than once' in run('reused',['--sampleList',str(reused)],False)
    missing=RUN/'missing.txt'; missing.write_text(SHEET.read_text().replace('S1_R1','missing_R1'))
    assert 'FASTQ missing' in run('missing',['--sampleList',str(missing)],False)
    run('strand_preflight',['--sampleList','__SHEET__','--transcriptome',str(TRANSCRIPTOME)],
        sample_sheet=LANES,pipeline=ROOT/'strandPreflight.nf')
    ref_root=RUN/'ref_root'; ref_root.mkdir()
    build_dir=RUN/'reference_build'; build_dir.mkdir()
    build_cmd=[NF,'run',str(ROOT/'buildReference.nf'),'-ansi-log','false','-c',str(CFG),
        '--reference_build_root',str(ref_root),'--reference_name','mock_ref',
        '--fasta_ref','https://example/ref.fa.gz','--GTF','https://example/ref.gtf.gz',
        '--cdna','https://example/cdna.fa.gz','--known_vcf','https://example/known.vcf.gz']
    built=subprocess.run(build_cmd,cwd=build_dir,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
        text=True,env={**os.environ,'NXF_DISABLE_CHECK_LATEST':'true'},timeout=180)
    (build_dir/'engine.log').write_text(built.stdout)
    if built.returncode:
        print(built.stdout); raise AssertionError(f'reference build failed; {build_dir}')
    built_ref=ref_root/'mock_ref'
    for required in ['Genome','genomeParameters.txt','ref.fa','ref.gtf','ref.ExonsOnly.bed',
                     'kallisto_index','VEP/homo_sapiens/113_GRCh38/info.txt']:
        assert (built_ref/required).exists(), f'missing reference output {required}'
    print('PASS reference_build')

    rscript=shutil.which('Rscript')
    if rscript:
        import gzip
        aggregate_results=RUN/'aggregate_fixture'/'results'
        counts_dir=aggregate_results/'FeatureCounts_output'; logs_dir=aggregate_results/'StarLog_output'
        counts_dir.mkdir(parents=True); logs_dir.mkdir(parents=True)
        for sample, values in [('S1',(11,21)),('S2',(12,22))]:
            header='Geneid\tChr\tStart\tEnd\tStrand\tLength\t'+sample+'StarOutAligned.sortedByCoord.out.bam\n'
            for suffix, counts in [('genecount.txt',values),('exonscount.txt',(v*10 for v in values))]:
                with (counts_dir/f'{sample}{suffix}').open('w') as f:
                    f.write('# featureCounts fixture\n'+header)
                    for gene, count in zip(['g1','g2'],counts):
                        f.write(f'{gene}\tchr1\t1\t10\t+\t10\t{count}\n')
            (logs_dir/f'{sample}StarOutLog.final.out').write_text(
                'Number of input reads | 100\nUniquely mapped reads number | 80\n'
                'Uniquely mapped reads % | 80.00%\n')
        aggregation=subprocess.run([rscript,'--vanilla',str(ROOT/'PrePostScripts/aggregSTAR.R'),
            str(aggregate_results),'2'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60)
        if aggregation.returncode:
            print(aggregation.stdout); raise AssertionError('R aggregation failed')
        with gzip.open(aggregate_results/'StarFCGeneCount_s2.tsv.gz','rt') as f:
            gene_matrix=f.read()
        assert gene_matrix.splitlines()[0] == 'S1\tS2'
        assert 'g1\t11\t12' in gene_matrix and 'g2\t21\t22' in gene_matrix
        with gzip.open(aggregate_results/'StarFCExonCount_s2.tsv.gz','rt') as f:
            exon_matrix=f.read()
        assert 'g1\t110\t120' in exon_matrix and 'g2\t210\t220' in exon_matrix
        metrics=(aggregate_results/'StarMetrics_s2.tsv').read_text()
        assert 'Uniquelymappedreads' in metrics and '80.00' in metrics
        print('PASS aggregation_script')
    else:
        print('SKIP aggregation_script (Rscript is not installed)')
    print('All sample-identity checks passed')
except Exception:
    print(f'Diagnostic files retained: {RUN}')
    raise
else:
    shutil.rmtree(RUN)
