"""Unit and optional real-Bowtie2 regression tests for RBC post-analysis."""
import copy
import csv
import json
import random
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
import os
import sqlite3

import csrna_rbc_diagnostics as d


def args(root):
    return d.parser().parse_args(['--analysis-report',str(root/'main.json'),'--manifest',str(root/'samples.csv'),
        '--reference',str(root/'ref.fa'),'--output-dir',str(root/'out'),'--index-cache',str(root/'indexes'),
        '--threads','1','--pilot-pairs','0'])


def write_fastq(path, seqs, mate=1, name='read'):
    path.write_text(''.join(f'@{name}{i}/{mate}\n{s}\n+\n'+('I'*len(s))+'\n' for i,s in enumerate(seqs)))


class DiagnosticTests(unittest.TestCase):
    def test_complexity(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=args(Path(tmp))
            self.assertEqual(d.complexity('A'*30)['entropy_bits'],0)
            self.assertEqual(d.complexity('ACGT'*20)['entropy_bits'],2)
            self.assertIn('low_3mer_complexity',d.complexity_labels(d.complexity('AC'*75),a))
            self.assertNotIn('homopolymer_rich',d.complexity_labels(d.complexity('AC'*75),a))
            self.assertIn('terminal_poly_G',d.complexity_labels(d.complexity('ACGT'+'G'*20),a))
            self.assertEqual(d.complexity('NNNN')['n_fraction'],1)

    def test_clip_orientation_and_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            a=args(Path(tmp));s='A'*100
            fields=['q1','0','ref','1','0','10S80M10S','*','0','0',s,'*','NM:i:2','AS:i:140']
            h=d.parse_hit(fields,s,a);self.assertEqual(h[-3:],(10,10,1))
            fields[1]='16';fields[5]='5S80M15S'
            h=d.parse_hit(fields,s,a);self.assertEqual(h[-3:],(15,5,1))
            fields[5]='5S30M65S';self.assertEqual(d.parse_hit(fields,s,a)[-1],0)
            fields[5]='100M';fields[-2]='NM:i:11';self.assertEqual(d.parse_hit(fields,s,a)[-1],0)
            fields[5]='100S1M'
            with self.assertRaises(ValueError):d.parse_hit(fields,s,a)

    def test_mode_partition(self):
        no={'accepted_best_loci':0,'parent':None};a={'accepted_best_loci':1,'parent':'a'};b={'accepted_best_loci':1,'parent':'b'}
        amb={'accepted_best_loci':2,'parent':None}
        self.assertEqual(d.mode_state(no,no),'neither_accepted')
        self.assertEqual(d.mode_state(a,no),'end_to_end_only')
        self.assertEqual(d.mode_state(no,a),'local_only')
        self.assertEqual(d.mode_state(a,a),'both_same_parent')
        self.assertEqual(d.mode_state(a,b),'both_different_parent')
        self.assertEqual(d.mode_state(a,amb),'both_parent_ambiguous')

    def test_pair_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);r1=root/'r1.fq';r2=root/'r2.fq'
            write_fastq(r1,['ACGT']*2,1);write_fastq(r2,['TGCA']*2,2)
            self.assertEqual(len(list(d.pair_records(r1,r2))),2)
            write_fastq(r2,['TGCA'],2)
            with self.assertRaises(ValueError):list(d.pair_records(r1,r2))
            write_fastq(r2,['TGCA']*2,1)
            with self.assertRaises(ValueError):list(d.pair_records(r1,r2))

    def test_sampling_resume_and_multiplicity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=args(root);a.pilot_pairs=17
            write_fastq(root/'r1.fq',['ACGT']*100,1);write_fastq(root/'r2.fq',['TGCA']*100,2)
            row={'sample_id':'RBC_1','r1':str(root/'r1.fq'),'r2':str(root/'r2.fq')}
            db=d.connect(root/'db.sqlite');v=d.ingest(db,row,a,root)
            self.assertEqual(v['validated_pairs'],100);self.assertEqual(v['analyzed_pairs'],17)
            self.assertEqual(db.execute('SELECT n FROM counts WHERE mate=1').fetchone()[0],17)
            ordinals=[r[0] for r in db.execute('SELECT ordinal FROM pairs ORDER BY ordinal')]
            self.assertGreater(max(ordinals),17)
            self.assertEqual(d.ingest(db,row,a,root),v)
            db.close()
            db=d.connect(root/'other.sqlite');d.ingest(db,row,a,root)
            self.assertEqual(ordinals,[r[0] for r in db.execute('SELECT ordinal FROM pairs ORDER BY ordinal')]);db.close()

    def test_manifest_rbc_only_and_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);write_fastq(root/'a.fq',['ACGT'],1);write_fastq(root/'b.fq',['ACGT'],2)
            p=root/'samples.csv';p.write_text('sample_id,condition,r1,r2\nRBC_1,RBC,a.fq,b.fq\nM1_1,M1,missing1,missing2\n')
            r={'samples':{'RBC_1':{}},'sample_metadata':[{'sample_id':'RBC_1','condition':'RBC'}]}
            self.assertEqual(len(d.selected_manifest(p,r)),1)
            r['sample_metadata'][0]['condition']='M1'
            with self.assertRaises(ValueError):d.selected_manifest(p,r)

    def test_missing_pilot_blocks_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=args(root);a.pilot_pairs=1
            write_fastq(root/'r1.fq',['ACGT'],1);write_fastq(root/'r2.fq',['TGCA'],2)
            row={'sample_id':'RBC_1','r1':str(root/'r1.fq'),'r2':str(root/'r2.fq')}
            db=d.connect(root/'db.sqlite');v=d.ingest(db,row,a,root)
            Path(v['retained_pilot_fastqs'][0]['path']).write_text('corrupted')
            with self.assertRaises(ValueError):d.ingest(db,row,a,root)
            db.close()

    def test_partial_ingestion_restarts_without_double_counting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=args(root)
            write_fastq(root/'r1.fq',['ACGT']*4,1);write_fastq(root/'r2.fq',['TGCA']*4,2)
            row={'sample_id':'RBC_1','r1':str(root/'r1.fq'),'r2':str(root/'r2.fq')}
            db=d.connect(root/'db.sqlite');d.ingest(db,row,a,root)
            db.execute('DELETE FROM done');db.commit()
            d.ingest(db,row,a,root)
            self.assertEqual(db.execute('SELECT sum(n) FROM counts').fetchone()[0],8)
            self.assertEqual(db.execute('SELECT count(*) FROM pairs').fetchone()[0],4);db.close()

    def test_locks_and_atomic_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with d.lock(p):
                with self.assertRaises(ValueError):
                    with d.lock(p): pass
            self.assertFalse((p/'.lock').exists())
            d.atomic_json(p/'r.json',{'a':1})
            with self.assertRaises(ValueError):d.atomic_json(p/'r.json',{'a':float('nan')})
            self.assertEqual(json.loads((p/'r.json').read_text()),{'a':1})


@unittest.skipUnless(os.environ.get('CSRNA_TEST_TOOLS')=='1','Set CSRNA_TEST_TOOLS=1 for actual Bowtie2 integration')
class BowtieTests(unittest.TestCase):
    def test_real_alignment_resume_and_screen(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);rng=random.Random(18)
            ref=''.join(rng.choice('ACGT') for _ in range(1000));unknown=''.join(rng.choice('ACGT') for _ in range(100))
            (root/'ref.fa').write_text('>ribo|test rRNA\n'+ref+'\n')
            (root/'screen.fa').write_text('>technical_test\n'+unknown+'\n')
            seqs=[ref[100:200],ref[100:200],('T'*25)+ref[300:375],unknown,'A'*100]
            write_fastq(root/'r1.fq',seqs,1);write_fastq(root/'r2.fq',[d.reverse_complement(s) for s in seqs],2)
            (root/'samples.csv').write_text('sample_id,condition,r1,r2,strandedness,pcr_cycles,pcr_cycles_document_reported,pcr_cycles_review_status\nRBC_1,RBC,r1.fq,r2.fq,strand_specific_orientation_unverified,15,24,conflicting_sources\n')
            (root/'main.json').write_text(json.dumps({'schema_version':'2.1.1','status':'complete','samples':{'RBC_1':{}},
                 'sample_metadata':[{'sample_id':'RBC_1','condition':'RBC'}],'reference':{'sha256':d.sha256(root/'ref.fa')}}))
            a=args(root);a.screen_reference=['technical='+str(root/'screen.fa')]
            r=d.run(a);s=r['samples']['RBC_1/R1']
            self.assertEqual(r['sample_metadata'][0]['pcr_cycles'],'15')
            self.assertEqual(r['sample_metadata'][0]['pcr_cycles_document_reported'],'24')
            self.assertEqual(r['sample_metadata'][0]['pcr_cycles_review_status'],'conflicting_sources')
            self.assertNotIn('pcr_cycles',r['main_analysis_sample_metadata'][0])
            self.assertNotIn('pcr_cycles',json.loads((root/'main.json').read_text())['sample_metadata'][0])
            self.assertEqual(s['reads'],5);self.assertEqual(s['distinct_sequences'],4)
            self.assertEqual(s['states']['both_same_parent'],2)
            self.assertEqual(s['states']['local_only'],1)
            self.assertEqual(s['states']['neither_accepted'],2)
            self.assertEqual(s['origin_screen_matches_neither']['technical'],1)
            self.assertEqual(sum(r['pair_agreement']['RBC_1']['local'].values()),5)
            self.assertEqual(s['labels']['homopolymer_rich']['reads'],2)
            with patch.object(d,'checked',side_effect=AssertionError('index rebuilt')),patch.object(d.subprocess,'Popen',wraps=d.subprocess.Popen):
                resumed=d.run(a)
            self.assertEqual(r['samples'],resumed['samples'])
            # A second raw reference with an identical sequence must not create unique-parent counts.
            (root/'duplicate.fa').write_text('>ribo|a\n'+ref+'\n>ribo|b\n'+ref+'\n')
            main=json.loads((root/'main.json').read_text());main['reference']['sha256']=d.sha256(root/'duplicate.fa')
            (root/'duplicate_main.json').write_text(json.dumps(main))
            alternate=copy.deepcopy(a);alternate.reference=str(root/'duplicate.fa');alternate.analysis_report=str(root/'duplicate_main.json');alternate.output_dir=str(root/'duplicate_out')
            dup=d.run(alternate)
            self.assertEqual(dup['samples']['RBC_1/R1']['states']['both_parent_ambiguous'],2)
            self.assertEqual(len((root/'duplicate_out/parent_counts.tsv').read_text().splitlines()),1)
            # Original source identity remains verifiable, even after moving paths.
            main=json.loads((root/'main.json').read_text());main['paired_input_qc']={'RBC_1':{'r1':'remote/r1','r2':'remote/r2'}}
            main['input_provenance']={'remote/r1':'0'*64}
            (root/'bad_main.json').write_text(json.dumps(main))
            bad=copy.deepcopy(a);bad.analysis_report=str(root/'bad_main.json');bad.output_dir=str(root/'bad_out')
            with self.assertRaises(ValueError):d.run(bad)
            self.assertFalse((root/'bad_out/rbc_diagnostics.json').exists())
            a.min_aligned_fraction=.9
            with self.assertRaises(ValueError):d.run(a)


if __name__=='__main__': unittest.main(verbosity=2)
