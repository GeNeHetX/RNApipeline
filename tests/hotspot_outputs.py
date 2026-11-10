#!/usr/bin/env python3
"""Check reference strand/phase handling and the agreed reporting rules."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'PrePostScripts'))
from hotspot_reference import build, validate
from hotspot_report import report, histogram_median


class HotspotOutputs(unittest.TestCase):
    def test_reference_strands_and_unresolved(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'ref.fa').write_text('>1\nATGGTGNNCACCAT\n')
            pysam.faidx(str(p/'ref.fa'))
            (p/'ref.gtf').write_text(''.join(
                f'1\tfixture\tCDS\t{start}\t{end}\t.\t{strand}\t0\tgene_name "{gene}"; transcript_id "{tx}"; transcript_version "1"; tag "Ensembl_canonical";\n'
                for gene,tx,start,end,strand in [('POS','TP',1,6,'+'),('NEG','TN',9,14,'-')]))
            (p/'catalogue.tsv').write_text('id\tgene\tfeature\tvariant_class\nP\tPOS\tM1V\tprotein_snv\nN\tNEG\tM1V\tprotein_snv\nBAD\tPOS\tA2V\tprotein_snv\n')
            build(p/'catalogue.tsv',p/'ref.gtf',p/'ref.fa',p/'out')
            alleles=json.loads((p/'out/candidate_alleles.json').read_text())
            self.assertEqual([(v['chrom'],v['pos'],v['ref'],v['alt']) for v in alleles],
                             [('1',1,'A','G'),('1',14,'T','C')])
            with (p/'out/mapping_audit.tsv').open() as f:
                audit=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(audit[-1]['status'],'unresolved')
            targets=json.loads((p/'out/targets.json').read_text())
            self.assertEqual(targets['NEG']['hotspot'],[[11,14]])
            results=[dict(input=f"{v['chrom']}\t{v['pos']}\t{v['id']}\t{v['ref']}\t{v['alt']}",
                          transcript_consequences=[dict(transcript_id=v['transcript'],protein_start=1,amino_acids='M/V')]) for v in alleles]
            (p/'vep.json').write_text(''.join(json.dumps(v)+'\n' for v in results))
            validate(p/'out/candidate_alleles.json',p/'vep.json',p/'validated.json')
            results[0]['transcript_consequences'][0]['amino_acids']='M/L'
            (p/'vep.json').write_text(''.join(json.dumps(v)+'\n' for v in results))
            with self.assertRaises(ValueError):
                validate(p/'out/candidate_alleles.json',p/'vep.json',p/'validated.json')

    def test_coverage_rules_and_alt_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            targets={
                'HOT':dict(chrom='1',transcript='TH.1',hotspot=[[0,3]],coding=[]),
                'DUAL':dict(chrom='1',transcript='TD.1',hotspot=[[0,3]],coding=[[3,7]]),
                'DUAL_HIGH':dict(chrom='1',transcript='TX.1',hotspot=[[0,2]],coding=[[7,11]]),
                'BROAD':dict(chrom='1',transcript='TB.1',hotspot=[],coding=[[3,7]]),
                'UNASSESSABLE':dict(chrom='1',transcript='TU.1',category='dual',hotspot=[],coding=[[3,7]]),
                'CALL':dict(chrom='2',transcript='TC.1',hotspot=[[0,3]],coding=[])}
            (p/'targets.json').write_text(json.dumps(targets))
            # Missing hotspot position 3 and coding position 7 both mean zero.
            (p/'depth.tsv').write_text('1\t1\t9\n1\t2\t7\n1\t4\t10\n1\t5\t20\n1\t6\t30\n1\t8\t1\n1\t9\t2\n1\t10\t3\n')
            (p/'input.vcf').write_text('##fileformat=VCFv4.2\n'
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Feature|Consequence|HGVSc|HGVSp">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
                '2\t1\t.\tA\tG\t50\tPASS\tCSQ=TC.1|missense_variant|TC.1:c.1A>G|PC:p.Met1Val\tGT:AD:DP\t0/1:16,4:20\n')
            report('S',p/'targets.json',p/'depth.tsv',p/'input.vcf',p/'out.tsv')
            with (p/'out.tsv').open() as f:
                reader=csv.DictReader(f,delimiter='\t'); rows={r['genesymbol']:r for r in reader}
                self.assertEqual(reader.fieldnames,['sample_id','genesymbol','mutated_reads','total_reads','cdna_mutation','protein_mutation'])
            self.assertEqual(rows['HOT']['total_reads'],'0')
            self.assertEqual(rows['DUAL']['total_reads'],'0')
            self.assertEqual(rows['DUAL_HIGH']['total_reads'],'7')
            self.assertEqual(rows['BROAD']['total_reads'],'15.0')
            self.assertEqual((rows['CALL']['mutated_reads'],rows['CALL']['total_reads']),('4','20'))
            self.assertEqual(rows['HOT']['protein_mutation'],'')
            self.assertNotIn('UNASSESSABLE',rows)


if __name__=='__main__': unittest.main()
