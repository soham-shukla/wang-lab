import copy
import os
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).parent))
from model import validate, composition, feature_frame, comparison, compact, csv_bytes
from streamlit.testing.v1 import AppTest
import pandas as pd


def fixture(scope='full_library'):
    features=[]
    for i,category,counts in [('f1','rRNA',{'A':(10,2),'B':(20,4)}),('f2','snoRNA',{'A':(0,0),'B':(0,8)})]:
        features.append({'id':i,'parent_RNA':'parent_'+i,'category':category,'strand':'+','start':0,'end':30,
                         'counts':{s:{'exact':a,'inexact':b,'inclusive':a+b} for s,(a,b) in counts.items()},
                         'observed_sequence':'A'*30})
    samples={s:{'total_input_reads':n,'rpm_denominator':n,'mapped_reads':n//2,'unmapped_reads':n//2,
                 'category_counts':{'rRNA':n//4,'snoRNA':n//4,'unmapped':n//2},'category_retention':{}}
             for s,n in [('A',100),('B',200)]}
    return {'schema_version':'2.0.0','parameters':{'input_scope':scope,'denominator':'total_input','pseudocount_rpm':1},
            'samples':samples,'features':features,'reference':{'records':{'irrelevant':{'header':'large'}}},
            'count_matrices':{t:{'feature_ids':['f1','f2'],'sample_ids':['A','B'],
                                'data':[[f['counts'][s][t] for s in samples] for f in features]}
                              for t in ('exact','inexact','inclusive')}}


class ModelTests(unittest.TestCase):
    def test_validation_and_axes(self):
        r=fixture();validate(r)
        r['count_matrices']['exact']['sample_ids'].reverse()
        with self.assertRaises(ValueError):validate(r)

    def test_malformed_structure_and_zero_mapped(self):
        for broken in ([], {'schema_version':'2.0.0','samples':{'A':None},'features':[]},
                       {'schema_version':'99.0','samples':{},'features':[]}):
            with self.assertRaises(ValueError):validate(broken)
        r=fixture();r['parameters']['pseudocount_rpm']=float('nan')
        with self.assertRaises(ValueError):validate(r)
        r=fixture();r['features']=[]
        for m in r['count_matrices'].values():m['feature_ids']=[];m['data']=[]
        for row in r['samples'].values():
            row['mapped_reads']=0;row['unmapped_reads']=row['total_input_reads']
            row['category_counts']={'unmapped':row['total_input_reads']}
        validate(r);self.assertTrue(composition(r,'Mapped reads only').empty)

    def test_count_conservation(self):
        r=fixture();r['features'][0]['counts']['A']['inclusive']+=1
        with self.assertRaises(ValueError):validate(r)
        r=fixture();r['samples']['A']['category_counts']['unmapped']+=1
        with self.assertRaises(ValueError):validate(r)

    def test_denominators_and_hidden_rrna(self):
        r=fixture();a=composition(r);b=composition(r,omit_rrna=True)
        self.assertEqual(a[a.Sample=='A'].Percent.sum(),100)
        self.assertEqual(b[b.Sample=='A'].Percent.sum(),75)
        mapped=composition(r,'Mapped reads only')
        self.assertNotIn('unmapped',set(mapped.Category))
        self.assertEqual(mapped[mapped.Sample=='A'].Percent.sum(),100)

    def test_comparison_control_only_and_reciprocity(self):
        r=fixture();a=comparison(r,'A','B','inclusive').set_index('Feature')
        b=comparison(r,'B','A','inclusive').set_index('Feature')
        self.assertEqual(a.loc['f1','log2 fold change'],0)
        self.assertEqual(a.loc['f2','Detection'],'Denominator only')
        self.assertAlmostEqual(a.loc['f2','log2 fold change'],-b.loc['f2','log2 fold change'])
        with self.assertRaises(ValueError):comparison(fixture('selected_subset'),'A','B','inclusive')

    def test_identical_biological_profiles_have_no_ordination(self):
        from model import replicate_summary
        for size in [3,6,7,17]:
            r=fixture()
            names=[f'R{i}' for i in range(size)]
            r['samples']={s:copy.deepcopy(r['samples']['B']) for s in names}
            for f in r['features']:
                f['counts']={s:copy.deepcopy(f['counts']['B']) for s in names}
            r['sample_metadata']=[{'sample_id':s,'condition':'group','replicate_type':'biological'} for s in names]
            meta,corr,scores,variance,_=replicate_summary(r)
            self.assertTrue(scores.empty)
            self.assertEqual(variance,[])
            self.assertTrue((meta.replicate_type=='biological').all())
            self.assertAlmostEqual(corr.iloc[0,-1],1)

    def test_compaction_and_no_mutation(self):
        r=fixture();before=copy.deepcopy(r);small=compact(r)
        self.assertFalse(small['reference']['records']);feature_frame(r);composition(r);comparison(r,'A','B','exact')
        self.assertEqual(r,before)

    def test_empty_features_and_csv(self):
        r=fixture();r['features']=[]
        for m in r['count_matrices'].values():m['feature_ids']=[];m['data']=[]
        validate(r);self.assertTrue(feature_frame(r).empty);self.assertTrue(comparison(r,'A','B','exact').empty)
        self.assertIn("'=1+1",csv_bytes(pd.DataFrame({'x':['=1+1']})).decode('utf-8-sig'))


    def test_replicate_qc_numeric_and_degenerate_cases(self):
        import numpy as np
        from model import replicate_summary
        r=fixture()
        r['sample_metadata']=[{'sample_id': 'B', 'condition': 'group', 'replicate_type': 'unknown', 'pcr_cycles': '24'},
                              {'sample_id': 'A', 'condition': 'group', 'replicate_type': 'unknown', 'pcr_cycles': '15'}]
        original=copy.deepcopy(r)
        meta,corr,scores,variance,recurrence=replicate_summary(r)
        self.assertEqual(meta.sample_id.tolist(),['A','B'])
        self.assertAlmostEqual(corr.loc['A','B'],1)
        self.assertAlmostEqual(sum(variance),1)
        self.assertAlmostEqual(scores.PC1.sum(),0)
        self.assertEqual(recurrence.set_index('Feature').loc['f2','Samples detected'],1)
        self.assertEqual(r,original)
        # A single sample has no ordination; empty category/sample selections are safe.
        self.assertTrue(replicate_summary(r,selected=['A'])[2].empty)
        self.assertTrue(replicate_summary(r,categories=[])[1].isna().all().all())
        self.assertTrue(replicate_summary(r,selected=[])[0].empty)
        # Independent distance check of the PCA coordinates, using the documented transform.
        expected=np.log2(1+np.array([[120000,0],[120000,40000]]))
        self.assertAlmostEqual(np.linalg.norm(scores[['PC1','PC2']].iloc[0]-scores[['PC1','PC2']].iloc[1]),np.linalg.norm(expected[0]-expected[1]))


class AppTests(unittest.TestCase):
    def setUp(self):
        import tempfile,json
        self.temp=tempfile.TemporaryDirectory();self.path=Path(self.temp.name)/'report.json'
        self.path.write_text(json.dumps(fixture()))
        self.previous=os.environ.get('CSRNA_REPORT');os.environ['CSRNA_REPORT']=str(self.path)
    def tearDown(self):
        if self.previous is None:os.environ.pop('CSRNA_REPORT',None)
        else:os.environ['CSRNA_REPORT']=self.previous
        self.temp.cleanup()
    def test_navigation_and_empty_filters(self):
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        self.assertFalse(a.exception)
        for page in ['RNA composition','Replicate QC','Feature comparisons','Parent explorer','Feature library','Overview']:
            a.sidebar.radio[0].set_value(page).run();self.assertFalse(a.exception, page)
        a.sidebar.radio[0].set_value('Feature comparisons').run()
        a.multiselect[0].set_value([]).run();self.assertFalse(a.exception)
        a.sidebar.radio[0].set_value('Feature library').run()
        next(x for x in a.text_input if x.label.startswith('Find')).set_value('absent-query').run();self.assertFalse(a.exception)
    def test_empty_report_views_and_chart_export(self):
        import json
        r=fixture();r['features']=[]
        for m in r['count_matrices'].values():m['feature_ids']=[];m['data']=[]
        self.path.write_text(json.dumps(r))
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        next(b for b in a.button if b.label=='Prepare interactive HTML export').click().run()
        self.assertFalse(a.exception)
        for page in ['Replicate QC','Feature comparisons','Parent explorer','Feature library']:
            a.sidebar.radio[0].set_value(page).run();self.assertFalse(a.exception,page)

    def test_stacked_chart_offsets(self):
        import json
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        a.sidebar.radio[0].set_value('RNA composition').run()
        self.assertFalse(a.exception)
        fig=json.loads(a.get('plotly_chart')[0].proto.spec)
        self.assertEqual(fig['layout']['barmode'],'relative')
        self.assertTrue(all(t.get('offsetgroup')=='stack' for t in fig['data']))

    def test_invalid_and_subset_reports(self):
        import json
        self.path.write_text(json.dumps(fixture('selected_subset')))
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        a.sidebar.radio[0].set_value('Feature comparisons').run()
        self.assertFalse(a.exception);self.assertFalse(any(x.label=='Numerator sample' for x in a.selectbox))
        self.path.write_text('{invalid')
        a.run();self.assertFalse(a.exception);self.assertTrue(a.error)

    def test_replicate_empty_selection_and_exploratory_banner(self):
        import json
        r=fixture()
        r['parameters']['analysis_mode']='exploratory'
        r['sample_metadata']=[{'sample_id':s,'condition':'group','replicate_type':'unknown'} for s in r['samples']]
        r['paired_input_qc']={'A':{'validated_pairs':100,'analyzed_r1_reads':100}}
        self.path.write_text(json.dumps(r))
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        a.sidebar.radio[0].set_value('Replicate QC').run()
        self.assertFalse(a.exception)
        self.assertTrue(any('R1 analyzed' in x.value for x in a.caption))
        a.multiselect[1].set_value([]).run()
        self.assertFalse(a.exception)
        a.multiselect[0].set_value([]).run()
        self.assertFalse(a.exception)

    def test_fifteen_library_manifest_layout(self):
        import json
        r=fixture()
        names=[f'{group}_{i}' for group in ['BMDC','BMDC_CpG','BMDM','M1','RBC'] for i in [1,2,3]]
        r['samples']={s:copy.deepcopy(r['samples']['A']) for s in names}
        for f in r['features']:
            f['counts']={s:copy.deepcopy(f['counts']['A']) for s in names}
        r['sample_metadata']=[{'sample_id':s,'condition':s.rsplit('_',1)[0],'replicate_type':'unknown'} for s in names]
        for track,m in r['count_matrices'].items():
            m['sample_ids']=names
            m['data']=[[f['counts'][s][track] for s in names] for f in r['features']]
        self.path.write_text(json.dumps(r))
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        a.sidebar.radio[0].set_value('Replicate QC').run()
        self.assertFalse(a.exception)
        self.assertTrue(any('Ordination unavailable' in x.value for x in a.info))
        fig=json.loads(a.get('plotly_chart')[0].proto.spec)
        self.assertEqual(len(fig['layout']['annotations'])-1,5)

    def test_parent_chart_long_labels_and_ranking(self):
        import json
        r=fixture()
        long_id='GENCODE|' + 'long_reference_identifier|' * 8
        r['features'][0]['parent_RNA']=long_id
        self.path.write_text(json.dumps(r))
        a=AppTest.from_file(str(Path(__file__).with_name('app.py')),default_timeout=30).run()
        a.sidebar.radio[0].set_value('RNA composition').run()
        self.assertFalse(a.exception)
        fig=json.loads(a.get('plotly_chart')[-1].proto.spec)
        axis=fig['layout']['yaxis']
        self.assertEqual(axis['categoryarray'][0],long_id)
        self.assertEqual(axis['autorange'],'reversed')
        self.assertTrue(all(len(label)<=32 for label in axis['ticktext']))
        self.assertGreaterEqual(fig['layout']['height'],550)
        self.assertEqual([t['name'] for t in fig['data']],['A','B'])
        self.assertIn(long_id,fig['data'][0]['y'])
        self.assertEqual(list(fig['data'][0]['x']),[120000,0])

if __name__=='__main__':unittest.main(verbosity=2)
