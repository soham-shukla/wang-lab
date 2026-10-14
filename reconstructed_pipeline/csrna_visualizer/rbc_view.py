"""Read-only standalone diagnostic-report view. No analysis-program imports."""
import json
import math
from pathlib import Path
import pandas as pd
import plotly.express as px
import streamlit as st
from model import csv_bytes


def validate_diagnostics(report):
    if not isinstance(report,dict) or report.get('schema_version')!='csrna-rbc-diagnostics/1' or report.get('status')!='complete':
        raise ValueError('Choose a completed RBC diagnostic JSON, not the main pipeline JSON.')
    try:
        if report['input_scope'] not in ('full_library','selected_subset') or not report['samples']:
            raise ValueError('Invalid diagnostic scope/samples')
        metadata={m['sample_id']:m for m in report['sample_metadata']}
        observed=set()
        for key,s in report['samples'].items():
            if s['mate'] not in (1,2) or key!=f"{s['sample']}/R{s['mate']}" or metadata[s['sample']]['condition']!='RBC':
                raise ValueError('Diagnostic samples must be RBC and have valid mate identities')
            observed.add((s['sample'],s['mate']))
            n=s['reads']
            if type(n) is not int or n<=0 or type(s['distinct_sequences']) is not int or not 0<s['distinct_sequences']<=n:
                raise ValueError('Invalid diagnostic read totals')
            if any(type(v) is not int or v<0 for v in s['states'].values()) or sum(s['states'].values())!=n:
                raise ValueError('Alignment-mode categories do not conserve reads')
            for mode in ('end_to_end','local'):
                accepted=s['accepted_match_reads'].get(mode,0)
                reported=s['reported_match_reads'].get(mode,0)
                if not 0<=accepted<=reported<=n or sum(s['categories_by_mode'][mode].values())!=accepted:
                    raise ValueError('Invalid mapped-read accounting')
            for metric in s['complexity_by_state'].values():
                if not all(math.isfinite(metric[k]) for k in ('length_sum','gc_sum','entropy_sum','mean_phred_sum')):
                    raise ValueError('Nonfinite complexity summary')
            if sum(v['reads'] for v in s['complexity_by_state'].values())!=n:
                raise ValueError('Complexity counts do not conserve reads')
        for name,a in report['input_audit'].items():
            if (name,1) not in observed or (name,2) not in observed:
                raise ValueError('Both mates must be represented')
            for mode,states in report['pair_agreement'][name].items():
                if mode not in ('end_to_end','local') or any(type(n) is not int or n<0 for n in states.values()) or sum(states.values())!=a['analyzed_pairs']:
                    raise ValueError('Pair summary does not conserve pairs')
            for mate in (1,2):
                if report['samples'][f'{name}/R{mate}']['reads']!=a['analyzed_pairs']:
                    raise ValueError('Mate totals differ from analyzed pair count')
            if not 0<a['analyzed_pairs']<=a['validated_pairs']:
                raise ValueError('Invalid analyzed/validated pair counts')
        if {name for name,mate in observed}!=set(report['input_audit']):
            raise ValueError('Audit and summary samples differ')
        return report
    except (KeyError,TypeError,AttributeError) as e:
        raise ValueError('Incomplete RBC diagnostic report') from e


def render():
    st.title('RBC post-analysis diagnostics')
    st.caption('Independent-mate alignment and sequence diagnostics. Main pipeline fragment counts are unchanged.')
    path=st.text_input('RBC diagnostic JSON path',help='Select rbc_diagnostics.json from the post-analysis output directory.')
    upload=st.file_uploader('Or upload RBC diagnostic JSON',type=['json'],key='rbc_upload')
    try:
        raw=upload.getvalue() if upload is not None else Path(path).expanduser().read_bytes() if path.strip() else None
        if raw is None:
            st.info('Run csrna_rbc_diagnostics.py after the main pipeline, then open its rbc_diagnostics.json here.')
            return
        report=validate_diagnostics(json.loads(raw))
    except (OSError,ValueError) as e:
        st.error(f'Cannot load RBC diagnostics: {e}');return
    st.warning('RBC-associated preparation: no surface-localization or contamination verdict. Local alignment ends and library insert boundaries do not establish native RNA ends; preparation fragmentation remains a separate question.')
    with st.expander('Preparation metadata and interpretation'):
        st.json(report.get('sample_metadata', []))
        st.caption(report.get('metadata_policy', 'Metadata is the snapshot stored in this report; later lab documentation is not automatically applied.'))
        st.caption('Strand-specific preparation does not verify read orientation. UDI indexes do not establish UMIs. Review conflicting PCR sources; no numerical PCR correction or automatic deduplication is applied.')
        if report.get('main_analysis_sample_metadata'):
            st.write('Original main-analysis metadata (preserved)')
            st.json(report['main_analysis_sample_metadata'])
    st.caption(f"Scope: {report['input_scope']} · R1 and R2 aligned independently · No trimming or abundance deduplication")
    options=list(report['samples'])
    selected=st.multiselect('RBC libraries and mates',options,default=options)
    if not selected:st.info('Select a library/mate.');return
    def table(df,key):
        st.dataframe(df,hide_index=True,use_container_width=True)
        st.download_button('Download table · CSV',csv_bytes(df.assign(Input_scope=report['input_scope'],Unit='sequenced reads unless pair-labeled')),key+'.csv','text/csv',key=key)
    def chart(fig,key):
        fig.update_layout(template='plotly_white',legend_title_text='')
        fig.add_annotation(text=f"{report['input_scope']} · Independent mate reads · Diagnostic evidence only",xref='paper',yref='paper',x=0,y=-.25,showarrow=False)
        fig.update_layout(margin={'b':120})
        st.plotly_chart(fig,use_container_width=True,key=key,config={'toImageButtonOptions':{'format':'svg','filename':key}})
        if st.button('Prepare interactive HTML export',key='export_'+key):
            st.download_button('Download chart · HTML',fig.to_html(include_plotlyjs=True),key+'.html','text/html',key='html_'+key)
    mode=[];totals=[];complexity=[];labels=[];clipping=[];origins=[];top=[];composition=[];parents=[]
    for key in selected:
        s=report['samples'][key];n=s['reads']
        totals.append({'Library/mate':key,'Analyzed reads':n,'Distinct sequences':s['distinct_sequences'],
            'End-to-end accepted':s['accepted_match_reads'].get('end_to_end',0),'Local accepted':s['accepted_match_reads'].get('local',0),
            'Local reported before diagnostic gates':s['reported_match_reads'].get('local',0)})
        for alignment_mode in ('end_to_end','local'):
            for category,count in s['categories_by_mode'][alignment_mode].items():
                composition.append({'Library/mate':key,'Mode':alignment_mode,'Category':category,'Percent':100*count/n})
            composition.append({'Library/mate':key,'Mode':alignment_mode,'Category':'not_accepted','Percent':100*(n-s['accepted_match_reads'].get(alignment_mode,0))/n})
            for parent in s.get('top_parents_by_mode',{}).get(alignment_mode,[]):
                parents.append({'Library/mate':key,'Mode':alignment_mode,**parent})
        for state,v in s['states'].items():mode.append({'Library/mate':key,'State':state,'Reads':v,'Percent':100*v/n})
        for state,v in s['complexity_by_state'].items():
            count=v['reads'];complexity.append({'Library/mate':key,'State':state,'Reads':count,'Distinct sequences':v['distinct_sequences'],
                'Mean length':v['length_sum']/count,'Mean GC fraction':v['gc_sum']/count,'Mean entropy (bits)':v['entropy_sum']/count,
                'Mean read Phred':v['mean_phred_sum']/count})
        for label,v in s['labels'].items():labels.append({'Library/mate':key,'Label':label,'Read percent':100*v['reads']/n,
            'Distinct-sequence percent':100*v['distinct_sequences']/s['distinct_sequences']})
        for clip,v in s['local_clipped_bases'].items():clipping.append({'Library/mate':key,'Clipped bases':int(clip),'Reads':v})
        for label,status in report['screen_status'].items():origins.append({'Library/mate':key,'Screen':label,'Status':status,
            'Matching reads':s['origin_screen_matches'].get(label,0) if status=='tested' else None,
            'Matching reads in neither-accepted group':s['origin_screen_matches_neither'].get(label,0) if status=='tested' else None})
        for item in s['top_sequences']:top.append({'Library/mate':key,**item})
    table(pd.DataFrame(totals),'rbc_totals')
    chart(px.bar(pd.DataFrame(mode),x='Library/mate',y='Percent',color='State',title='Alignment-mode comparison'),'rbc_modes')
    st.caption('States use accepted best alignments after the displayed coverage/edit gates. “Neither accepted” can include weak reported matches; it does not mean no sequence similarity. Scores are not compared across modes.')
    with st.expander('Alignment and complexity thresholds'):
        st.json({k:report['parameters'][k] for k in ['min_aligned_length','min_aligned_fraction','max_edit_rate','homopolymer_length','entropy_threshold','kmer_threshold']})
    st.subheader('RNA assignments by alignment mode')
    chart(px.bar(pd.DataFrame(composition),x='Library/mate',y='Percent',color='Category',facet_row='Mode',height=650,title='Accepted category support as percent of analyzed mate reads'),'rbc_composition')
    st.caption('Header-derived categories; ambiguity is retained. Lack of accepted support is not biological absence.')
    if parents:
        table(pd.DataFrame(parents),'rbc_parent_examples')
        st.caption('Top 20 unambiguous raw-reference parents per mate and mode. Complete parent counts are in parent_counts.tsv. Do not pool mates or mix these values with the main endpoint counts.')
    st.subheader('Complexity and quality by alignment outcome')
    table(pd.DataFrame(complexity),'rbc_complexity')
    if labels:
        chart(px.bar(pd.DataFrame(labels),x='Library/mate',y='Read percent',color='Label',barmode='group',title='Overlapping complexity labels'),'rbc_labels')
        table(pd.DataFrame(labels),'rbc_label_counts')
    else:st.info('No sequences meet the configured complexity-label thresholds.')
    st.caption('Labels overlap and do not imply contamination. 3-mer complexity is length-dependent; mean quality is diagnostic and is not used for sequence-level alignment scoring.')
    if clipping:chart(px.bar(pd.DataFrame(clipping).sort_values('Clipped bases'),x='Clipped bases',y='Reads',color='Library/mate',title='Clipping in accepted local alignments'),'rbc_clipping')
    st.caption('Clipping histogram uses one representative accepted best alignment per read; tied placements can have different clipping.')
    st.subheader('Origin-screen evidence')
    table(pd.DataFrame(origins),'rbc_origins')
    st.caption('Not tested is not zero. Matches are nonexclusive; conserved or low-complexity matches require review. Genome matching alone cannot distinguish DNA contamination from RNA.')
    pairrows=[]
    for sample,modes in report['pair_agreement'].items():
        if not any(key.startswith(sample+'/') for key in selected):continue
        for mode,states in modes.items():
            for state,n in states.items():pairrows.append({'Sample':sample,'Mode':mode,'Pair outcome':state,'Pairs':n})
    st.subheader('Mate agreement — pairs counted once')
    table(pd.DataFrame(pairrows),'rbc_pairs')
    st.caption('Same parent means independently aligned mates share a single raw reference ID; it does not establish orientation, insert size, or concordant paired alignment.')
    st.subheader('Abundant sequence examples')
    table(pd.DataFrame(top),'rbc_top_sequences')
    st.caption('Top 20 sequences per mate, ties broken by sequence. Full sequences and all best-hit evidence are retained in the SQLite database and compressed evidence table; these examples are not a complete unmapped catalogue.')
    with st.expander('Detailed output and provenance'):
        st.write('Sequence evidence:',report['sequence_evidence_path'])
        st.write('Parent counts:',report['parent_counts_path'])
        st.json(report['input_audit'])
        st.json(report['stages'])
        for note in report['interpretation']:st.write(note)
