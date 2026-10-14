"""Persistent counting regression tests; set CSRNA_TEST_TOOLS=1 for real alignment."""
import csv
import gzip
import json
import os
from pathlib import Path
import random
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
import csrna_pipeline as p


class CountCacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        rng = random.Random(72)
        self.reference = ''.join(rng.choice('ACGT') for _ in range(500))
        (self.root/'ref.fa').write_text('>ribo|test rRNA\n'+self.reference+'\n')
        self.sequences = [self.reference[20:80],self.reference[90:150],self.reference[20:80],self.reference[180:240],self.reference[90:150]]
        for mate in (1,2):
            with gzip.open(self.root/f'r{mate}.fq.gz','wt') as f:
                for i,s in enumerate(self.sequences):
                    seq = s if mate==1 else p.reverse_complement(s)
                    f.write(f'@read{i}/{mate}\n{seq}\n+\n'+chr(63+i)*len(seq)+'\n')
        (self.root/'samples.csv').write_text('sample_id,r1,r2,condition,replicate\nRBC_1,r1.fq.gz,r2.fq.gz,RBC,1\n')
        self.args = p.parser().parse_args(['--manifest',str(self.root/'samples.csv'),'--pair-handling','r1-only',
            '--analysis-mode','exploratory','--reference',str(self.root/'ref.fa'),'--count-cache',str(self.root/'cache'),
            '--count-checkpoint-records','2','--top-fragments','1','--max-fold-variants','1','--threads','1',
            '--output',str(self.root/'report.json')])
        self.rows = p.read_manifest(self.args.manifest)
        self.hashes = {Path(r[m]):p.sha256(r[m]) for r in self.rows for m in ('r1','r2')}

    def test_counts_quality_and_completed_resume(self):
        lib,qual,audit,stats = p.disk_count_inputs(self.rows,self.args,self.hashes)
        expected_q = {}; expected = p.load_counts(self.root/'r1.fq.gz',qualities=expected_q)
        self.assertEqual(lib['RBC_1'],expected);self.assertEqual(qual['RBC_1'],expected_q)
        self.assertEqual(stats['union_distinct_sequences'],3)
        self.assertEqual(audit['RBC_1']['validated_pairs'],5)
        with patch.object(p,'paired_r1_records',side_effect=AssertionError('revalidated completed sample')):
            self.assertEqual(p.disk_count_inputs(self.rows,self.args,self.hashes)[0],lib)

    def test_interrupted_batch_resumes_without_double_counts(self):
        real = p.paired_r1_records
        def interrupted(*a):
            for i,x in enumerate(real(*a)):
                if i==3: raise RuntimeError('interrupted after committed pair 2')
                yield x
        with patch.object(p,'paired_r1_records',interrupted):
            with self.assertRaises(RuntimeError):p.disk_count_inputs(self.rows,self.args,self.hashes)
        with sqlite3.connect(self.root/'cache/counts.sqlite') as db:
            self.assertEqual(db.execute('SELECT records FROM progress').fetchone()[0],2)
        observed = []
        def track(*a):
            observed.append(a[2]['records']);yield from real(*a)
        with patch.object(p,'paired_r1_records',track):lib,quality,_,_=p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.assertEqual(observed,[2]);self.assertEqual(sum(lib['RBC_1'].values()),5)
        q={};p.load_counts(self.root/'r1.fq.gz',qualities=q);self.assertEqual(quality['RBC_1'],q)

    def test_count_guard_does_not_abort_counting(self):
        self.args.max_distinct_sequences=1
        with self.assertRaisesRegex(ValueError,'Counts safely retained'):p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.args.count_only=True
        with patch.object(p,'paired_r1_records',side_effect=AssertionError('recount')):
            result=p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.assertEqual(result[3]['samples']['RBC_1']['reads'],5)
        self.args.count_only=False;self.args.max_distinct_sequences=3
        self.assertEqual(sum(p.disk_count_inputs(self.rows,self.args,self.hashes)[0]['RBC_1'].values()),5)

    def test_invalid_mate_never_completes(self):
        with gzip.open(self.root/'r2.fq.gz','at') as f:f.write('@extra/2\nAAAA\n+\nIIII\n')
        with self.assertRaises(ValueError):p.disk_count_inputs(self.rows,self.args,self.hashes)
        with sqlite3.connect(self.root/'cache/counts.sqlite') as db:
            self.assertEqual(db.execute('SELECT complete FROM progress').fetchone()[0],0)

    def test_invalidation_and_lock(self):
        p.disk_count_inputs(self.rows,self.args,self.hashes)
        changed=dict(self.hashes);changed[Path(self.rows[0]['r1'])]='0'*64
        with self.assertRaisesRegex(ValueError,'changed'):p.disk_count_inputs(self.rows,self.args,changed)
        self.args.quality_detail='counts-only'
        with self.assertRaisesRegex(ValueError,'changed'):p.disk_count_inputs(self.rows,self.args,self.hashes)
        (self.root/'cache/.count.lock').write_text('test')
        with self.assertRaisesRegex(ValueError,'locked'):p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.assertTrue((self.root/'cache/.count.lock').exists())

    def test_capacity_sampling_and_mode_validation(self):
        self.args.alignment_capacity_queries=2
        a=p.disk_count_inputs(self.rows,self.args,self.hashes)
        b=p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.assertEqual(a[0],b[0]);self.assertEqual(len(a[0]),2)
        self.assertEqual(a[3]['eligible_distinct_sequences'],3)
        self.args.count_only=True
        with self.assertRaises(ValueError):p.validate_args(self.args)

    def test_count_only_analyze_and_metadata_change_reuse(self):
        self.args.count_only=True
        r=p.analyze(self.args,{},self.root)
        self.assertEqual(r['schema_version'],'csrna-count-cache/1')
        (self.root/'samples.csv').write_text('sample_id,r1,r2,condition,replicate,pcr_cycles\nRBC_1,r1.fq.gz,r2.fq.gz,RBC,1,24\n')
        with patch.object(p,'paired_r1_records',side_effect=AssertionError('recount')):
            self.assertEqual(p.analyze(self.args,{},self.root)['counting']['samples'],r['counting']['samples'])

    def test_two_samples_and_global_guard(self):
        self.rows.append({**self.rows[0], 'sample_id':'RBC_2'})
        self.args.max_total_distinct_sequences=2
        with self.assertRaisesRegex(ValueError,'Counts safely retained'):
            p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.args.count_only=True
        report=p.disk_count_inputs(self.rows,self.args,self.hashes)[3]
        self.assertEqual(report['union_distinct_sequences'],3)
        self.assertEqual(sum(x['reads'] for x in report['samples'].values()),10)

    def test_counts_only_and_unrelated_database(self):
        self.args.quality_detail='counts-only'
        libraries, qualities, _, _=p.disk_count_inputs(self.rows,self.args,self.hashes)
        self.assertEqual(sum(libraries['RBC_1'].values()),5)
        self.assertEqual(qualities,{'RBC_1':{}})
        other=self.root/'unrelated';other.mkdir()
        with sqlite3.connect(other/'counts.sqlite') as db:db.execute('CREATE TABLE original (value TEXT)')
        self.args.count_cache=str(other)
        with self.assertRaisesRegex(ValueError,'refusing to modify'):
            p.disk_count_inputs(self.rows,self.args,self.hashes)
        with sqlite3.connect(other/'counts.sqlite') as db:
            self.assertEqual(db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall(),[('original',)])

    @unittest.skipUnless(os.environ.get('CSRNA_TEST_TOOLS')=='1','real tools opt-in')
    def test_real_pipeline_equivalence_and_capacity(self):
        tools=p.require_tools()
        for mode in ('disk','legacy','capacity'):(self.root/mode).mkdir()
        disk=p.analyze(self.args,tools,self.root/'disk')
        self.args.count_cache=None
        legacy=p.analyze(self.args,tools,self.root/'legacy')
        for key in ('features','count_matrices','parent_landscape','endpoint_profiles','support_reporting'):
            self.assertEqual(disk[key],legacy[key],key)
        self.args.count_cache=str(self.root/'cache');self.args.alignment_capacity_queries=2
        report=p.analyze(self.args,tools,self.root/'capacity')
        self.assertEqual(report['schema_version'],'csrna-alignment-capacity/1')
        self.assertEqual(report['capacity']['sampled_queries'],2)
        self.assertEqual(report['capacity']['mapped_queries'],2)


if __name__=='__main__':unittest.main(verbosity=2)
