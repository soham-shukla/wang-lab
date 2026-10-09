"""Launch: streamlit run reconstructed_pipeline/csrna_visualizer/app.py -- --report /path/to/report.json"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import TRACKS, compact, comparison, composition, csv_bytes, feature_frame, load_path, provenance, validate, replicate_summary

st.set_page_config(page_title='csRNA pipeline data', page_icon='🧬', layout='wide')
COLORS = {'rRNA': '#487C97', 'lncRNA': '#C07842', 'mRNA': '#749C6A', 'snoRNA': '#9B78B0',
          'snRNA': '#C45F70', 'tRNA': '#558D85', 'mitochondrial_RNA': '#A38C45', 'pseudogene': '#8D87AD',
          'unmapped': '#B8BDC5', 'unclassified_reference': '#8A919A', 'multi_category_ambiguous': '#D3B47C',
          'below_min_read_length': '#D9DCE1'}


@st.cache_data(max_entries=2, show_spinner='Reading and validating report…')
def cached_path(path, mtime, size):
    return load_path(path)


@st.cache_data(max_entries=2, show_spinner='Reading and validating report…')
def cached_upload(data):
    return compact(validate(json.loads(data)))


def table(df, name, evidence=None):
    st.dataframe(df, use_container_width=True, hide_index=True)
    exported = df.assign(Report_scope=report['parameters']['input_scope'], Evidence_track=evidence or track,
                         RPM_denominator=report['parameters'].get('denominator', 'unspecified'))
    st.download_button('Download table · CSV', csv_bytes(exported), name + '.csv', 'text/csv', key='csv_' + name)


def chart(fig, title, key, subtitle='', evidence=None):
    # Explicit shared offsets keep stacking stable across Plotly frontend versions.
    if fig.layout.barmode in ('relative', 'stack'):
        fig.update_traces(offsetgroup='stack', selector={'type': 'bar'})
    fig.update_layout(title={'text': title, 'font': {'size': 18}}, template='plotly_white',
                      font={'family': 'Arial', 'color': '#243444'}, paper_bgcolor='rgba(0,0,0,0)',
                      margin={'l': 35, 'r': 25, 't': 65, 'b': 100}, legend_title_text='')
    # Exported figures carry scope and normalization, rather than relying on page context.
    fig.add_annotation(text=provenance(report, evidence or track).replace(' · Samples:', '<br>Samples:'),
                       xref='paper', yref='paper', x=0, y=-0.25, showarrow=False, align='left', font={'size': 9})
    if key == 'parents':
        fig.update_layout(margin={'l': 20, 'r': 20, 't': 100, 'b': 145},
                          legend={'orientation': 'h', 'y': 1.08, 'x': 0},
                          title_text='Top 15 parents by combined eligible RPM')
        fig.layout.annotations[-1].update(
            text=provenance(report, evidence or track).replace(' · ', '<br>'), y=-0.12)
    st.plotly_chart(fig, use_container_width=True, key=key,
                    config={'displaylogo': False, 'toImageButtonOptions': {'format': 'svg', 'filename': key}})
    if subtitle:
        st.caption(subtitle)
    if st.button('Prepare interactive HTML export', key='export_' + key):
        st.download_button('Download chart · HTML', fig.to_html(include_plotlyjs=True), key + '.html', 'text/html', key='html_' + key)


args = argparse.ArgumentParser(add_help=False)
args.add_argument('--report', default=os.environ.get('CSRNA_REPORT', ''))
opts, _ = args.parse_known_args()
st.sidebar.markdown('## csRNA pipeline data')
page = st.sidebar.radio('Explore', ['Overview', 'RNA composition', 'Replicate QC', 'Feature comparisons', 'Parent explorer', 'Feature library'])
track = st.sidebar.selectbox('Evidence track', TRACKS, index=2,
                            help='Inclusive = exact + accepted-inexact support. No track proves modification.')
with st.sidebar.expander('Open a pipeline report', expanded=not bool(opts.report)):
    path = st.text_input('Local JSON path', value=opts.report, help='Recommended for reports larger than 200 MB. Paths are read locally.')
    uploaded = st.file_uploader('Or upload JSON', type=['json'])
    st.caption('Read-only connection. The pipeline is never executed or modified.')
report = None
try:
    if uploaded is not None:
        report = cached_upload(uploaded.getvalue())
    elif path.strip():
        p = Path(path.strip()).expanduser().resolve()
        stat = p.stat()
        report = cached_path(str(p), stat.st_mtime_ns, stat.st_size)
except (OSError, ValueError, TypeError, KeyError, OverflowError) as error:
    st.error(f'Cannot open report: {error}')
    st.stop()
if report is None:
    st.title('Your RNA landscape, in context')
    st.write('Open a pipeline JSON report in the sidebar to explore RNA classes, observed fragments, and sample comparisons.')
    st.info('Use a local path for large reports. No sequencing files or pipeline installation are required.')
    st.stop()

samples = list(report['samples'])
subset = report['parameters']['input_scope'] == 'selected_subset'
st.title(page)
st.caption(provenance(report, track))
if subset:
    st.warning('Selected subset · Charts describe supplied reads only. Biological fold-change rankings are disabled.')
else:
    st.info('Descriptive comparison · No replicate statistics. Preparation differences may affect comparisons; surface specificity is not established.')
frame = feature_frame(report, track)
if report.get('paired_input_qc'):
    st.caption('Paired libraries · R1 analyzed; R2 validated only. Counts represent R1 reads, not independently established molecules or reconstructed inserts.')
if report['parameters'].get('analysis_mode') == 'exploratory':
    st.caption('Exploratory analysis · No controls designated; control-dependent candidates are disabled.')
categories = sorted(frame['Category'].unique())
with st.sidebar.expander('Report details'):
    st.write('Schema:', report['schema_version'])
    st.write('Created:', report.get('created_utc', 'Not recorded'))
    st.caption('Unknown or absent reference annotations are not evidence of biological absence.')
    st.json(report.get('sample_metadata', []))
    if report.get('performance'):
        st.write('Measured pipeline stage durations (seconds)')
        st.json(report['performance']['stage_seconds'])
        st.caption(report['performance']['scope'])

if page == 'Overview':
    a, b, c = st.columns(3)
    a.metric('Supplied reads', f"{sum(s['total_input_reads'] for s in report['samples'].values()):,}")
    b.metric('Shared observed features', f"{len(frame):,}")
    c.metric('RNA classes with features', len(categories))
    st.write('Start with composition, compare supported features, then inspect their parent-RNA boundaries.')
    rows = []
    for s, v in report['samples'].items():
        rows.append({'Sample': s, 'Input reads': v['total_input_reads'], 'Unmapped reads': v['unmapped_reads'],
                     'Unmapped (%)': v['unmapped_reads'] * 100 / v['total_input_reads'],
                     'Eligible feature reads': sum(f['counts'][s][track] for f in report['features'])})
    df = pd.DataFrame(rows)
    chart(px.bar(df, x='Sample', y='Unmapped (%)', text_auto='.1f', color_discrete_sequence=['#487C97']),
          'How much of the input remains unmapped?', 'unmapped', 'Percent of all supplied reads; not a taxonomy or contamination diagnosis.', evidence='all input reads')
    table(df, 'library_qc')
    st.subheader('Evidence retained for features')
    support = pd.DataFrame([{'Sample': s, 'Track': t, 'Reads': sum(f['counts'][s][t] for f in report['features'])}
                            for s in samples for t in ('exact', 'inexact')])
    chart(px.bar(support, x='Sample', y='Reads', color='Track', color_discrete_map={'exact': '#487C97', 'inexact': '#C07842'}),
          'Exact and accepted-inexact support', 'support', 'Disjoint counts; inclusive support is their sum. Read counts are not independent molecule counts.', evidence='exact + inexact (disjoint)')
    if report.get('support_reporting', {}).get('samples'):
        resolved_rows = [{'Sample': s, 'Evidence': key.replace('_', ' '), 'Reads': values['counts'][key]}
                         for s, values in report['support_reporting']['samples'].items()
                         for key in ('raw_unique_exact', 'raw_unique_inexact', 'equivalence_resolved_exact', 'equivalence_resolved_inexact')]
        with st.expander('Raw-unique versus equivalence-resolved support'):
            support_resolution = pd.DataFrame(resolved_rows)
            chart(px.bar(support_resolution, x='Sample', y='Reads', color='Evidence'),
                  'How feature support was resolved', 'support_resolution',
                  'These four components are disjoint. Equivalence resolves reference redundancy, not the genomic copy of origin.', evidence='four disjoint support components')
            table(support_resolution, 'support_resolution', evidence='four disjoint support components')
    if report.get('unmapped_qc'):
        screens = [{'Sample': s, 'Screen': label, 'Status': v['status'], 'Matching reads': v['matching_reads']}
                   for s, qc in report['unmapped_qc'].items() for label, v in qc['screens'].items()]
        with st.expander('Unmapped investigation status'):
            table(pd.DataFrame(screens), 'diagnostic_status')
            st.caption('Screen matches may overlap. Not tested does not mean no match.')
    with st.expander('Mapping ambiguity and exclusions'):
        retention = [{'Sample': s, 'Category': cat, **values} for s, v in report['samples'].items()
                     for cat, values in v.get('category_retention', {}).items()]
        if retention:
            table(pd.DataFrame(retention).fillna(0), 'retention')
        st.caption('Raw ambiguity and retained support overlap: do not add these columns. Reference equivalence can recover endpoint support.')
    with st.expander('Pipeline notes'):
        for warning in report.get('warnings', []):
            st.write(warning)

elif page == 'RNA composition':
    denominator = st.radio('Composition denominator', ['All supplied reads', 'Mapped reads only'], horizontal=True)
    omit = st.checkbox('Hide rRNA to inspect other classes', help='Keeps the selected denominator; remaining bars need not sum to 100%.')
    df = composition(report, denominator, omit)
    if df.empty or df['Percent'].notna().sum() == 0:
        st.info('No reads are available for this composition view. Try all supplied reads or show rRNA.')
    else:
        chart(px.bar(df, x='Sample', y='Percent', color='Category', color_discrete_map=COLORS, hover_data=['Reads']),
              'RNA classes represented in the reads', 'composition', denominator + (' · rRNA hidden; denominator unchanged.' if omit else '.'), evidence='all mapping categories')
    table(df, 'composition', evidence='all mapping categories')
    st.subheader('Abundance and feature diversity answer different questions')
    rows = []
    for cat in categories:
        group = [f for f in report['features'] if f['category'] == cat]
        for s in samples:
            rows.append({'Category': cat, 'Sample': s, 'Detected features': sum(f['counts'][s][track] > 0 for f in group),
                         'Feature reads': sum(f['counts'][s][track] for f in group)})
    if rows:
        diversity = pd.DataFrame(rows)
        chart(px.bar(diversity, x='Category', y='Detected features', color='Sample', barmode='group'),
              'Observed feature diversity', 'diversity', 'Counts depend on sequencing depth, reference annotation, and eligibility. These are not unbiased species-richness estimates.')
        table(diversity, 'category_support')
    if not frame.empty:
        st.subheader('Which parent RNAs dominate eligible support?')
        parents = frame.groupby(['Parent RNA', 'Category'], as_index=False)[[f'{s} RPM' for s in samples]].sum()
        parents['_total'] = parents[[f'{s} RPM' for s in samples]].sum(axis=1)
        parents = parents.sort_values(['_total', 'Parent RNA'], ascending=[False, True]).head(15)
        plot_data = parents.drop(columns='_total').melt(id_vars=['Parent RNA', 'Category'], var_name='Sample', value_name='Eligible RPM')
        plot_data['Sample'] = plot_data['Sample'].str.removesuffix(' RPM')
        # Keep full reference IDs in the data and hover; bounded tick labels leave room for bars.
        parent_order = parents['Parent RNA'].tolist()
        labels = [p if len(p) <= 32 else p[:13] + '…' + p[-18:] for p in parent_order]
        parent_fig = px.bar(plot_data, y='Parent RNA', x='Eligible RPM', color='Sample',
                            barmode='group', orientation='h', hover_data=['Category'],
                            category_orders={'Parent RNA': parent_order, 'Sample': samples},
                            height=max(550, 42 * len(parents) + 245))
        parent_fig.update_yaxes(tickmode='array', tickvals=parent_order, ticktext=labels,
                               title=None, automargin=True, categoryorder='array',
                               categoryarray=parent_order, autorange='reversed')
        parent_fig.update_xaxes(rangemode='tozero', tickangle=0)
        chart(parent_fig,
              'Top 15 parents by combined eligible RPM', 'parents',
              'Sums disjoint eligible features. This does not include all mapped or ambiguously assigned reads.')
        table(parents.drop(columns='_total'), 'parent_abundance')


elif page == 'Replicate QC':
    st.write('Compare individual libraries before interpreting group differences. No samples are pooled or automatically excluded.')
    chosen_samples = st.multiselect('Samples to inspect', samples, default=samples)
    chosen_categories = st.multiselect('RNA classes for replicate QC', categories, default=categories)
    metadata, correlations, scores, variance, recurrence = replicate_summary(report, track, chosen_samples, chosen_categories)
    table(metadata, 'sample_metadata')
    st.caption('Replicate type, PCR cycles, and preparation fields are reported as supplied. Unknown fields are not inferred. PCR cycles do not supply a numerical correction.')
    if not chosen_samples:
        st.info('Select at least one sample.')
        st.stop()
    comp = composition(report)
    comp = comp[comp['Sample'].isin(chosen_samples)].merge(metadata[['sample_id', 'condition']], left_on='Sample', right_on='sample_id')
    fig = px.bar(comp, x='Sample', y='Percent', color='Category', facet_col='condition', facet_col_wrap=2,
                 color_discrete_map=COLORS, height=max(450, 350 * ((metadata['condition'].nunique() + 1) // 2)))
    fig.update_xaxes(matches=None)
    chart(fig, 'RNA composition by group and library', 'replicate_composition',
          'All supplied reads; includes unmapped and ineligible reads. The RNA-class filter affects feature-based QC only.', evidence='all mapping categories')
    if correlations.notna().any().any():
        chart(go.Figure(go.Heatmap(z=correlations.values, x=correlations.columns, y=correlations.index,
                                  zmin=-1, zmax=1, colorscale='RdBu', colorbar_title='Pearson r')),
              'Similarity of feature abundance profiles', 'replicate_correlation',
              'Pearson correlation of log2(1 + RPM), using features detected in at least one selected sample. Constant profiles are blank; high correlation does not establish biological replication.')
        table(correlations.rename_axis('Sample').reset_index(), 'replicate_correlations')
    else:
        st.info('Correlation unavailable: too few features or constant sample profiles.')
    if not scores.empty:
        fig = px.scatter(scores, x='PC1', y='PC2', color='condition', hover_name='sample_id',
                         hover_data=['replicate_type', 'pcr_cycles'], labels={'PC1': f'PC1 ({variance[0]:.1%})', 'PC2': f'PC2 ({variance[1]:.1%})'})
        chart(fig, 'Library ordination', 'replicate_pca',
              'PCA of centered log2(1 + RPM); features are not variance-scaled. Separation is descriptive and may reflect composition or preparation. Axis signs are arbitrary.')
        table(scores, 'replicate_pca')
    else:
        st.info('Ordination unavailable: select multiple samples with differing feature profiles.')
    if not recurrence.empty:
        st.subheader('Features recurring within each group')
        table(recurrence.sort_values(['Condition', 'Samples detected', 'Mean sample RPM'], ascending=[True, False, False]), 'group_feature_recurrence')
        st.caption('Detection means at least one eligible read. Recurrence is across library columns and does not prove independent biological or molecular support.')
    st.subheader('Raw counts for downstream statistics')
    matrix = report['count_matrices'][track]
    counts = pd.DataFrame(matrix['data'], index=matrix['feature_ids'], columns=matrix['sample_ids'])
    st.download_button('Download raw count matrix · CSV', csv_bytes(counts[chosen_samples].rename_axis('feature_id').reset_index()), f'{track}_raw_counts.csv', 'text/csv')
    st.caption('All shared features are exported as raw integer counts in selected-sample order; no RPM or pseudocounts. Match metadata by sample_id. Statistical design and biological replication still need confirmation.')
    if report.get('paired_input_qc'):
        table(pd.DataFrame([{'Sample': s, **v} for s, v in report['paired_input_qc'].items()]), 'paired_input_qc')


elif page == 'Feature comparisons':
    if frame.empty:
        st.info('No eligible features were reported. Inspect mapping and retention in Overview.')
        st.stop()
    chosen = st.multiselect('RNA classes', categories, default=categories)
    limit = st.slider('Features displayed in heatmap', 5, 100, 25)
    selected = frame[frame['Category'].isin(chosen)].copy()
    rpm_cols = [f'{s} RPM' for s in samples]
    selected['_total'] = selected[rpm_cols].sum(axis=1)
    selected = selected.sort_values(['_total', 'Feature'], ascending=[False, True]).head(limit)
    if not selected.empty:
        import numpy as np
        values = selected[rpm_cols].to_numpy(dtype=float)
        fig = go.Figure(go.Heatmap(z=np.log10(1 + values), x=samples, y=selected['Feature'], customdata=values,
                                  colorscale='Blues', colorbar_title='log10(1 + RPM)',
                                  hovertemplate='%{y}<br>%{x}: %{customdata:.3f} RPM<extra></extra>'))
        fig.update_layout(height=max(400, len(selected) * 19 + 200))
        chart(fig, 'Shared features across libraries', 'heatmap', 'Same absolute color scale across rows. Zero means no observed support in this track.')
        table(selected.drop(columns='_total'), 'shared_features')
    else:
        st.info('No features match the selected classes.')
    if subset:
        st.info('Pairwise fold rankings require a full-library report. This subset heatmap is descriptive only.')
    elif len(samples) >= 2:
        a, b = st.columns(2)
        numerator = a.selectbox('Numerator sample', samples)
        denominator = b.selectbox('Denominator sample', [s for s in samples if s != numerator])
        df = comparison(report, numerator, denominator, track)
        df = df[df['Category'].isin(chosen)]
        if not df.empty:
            import numpy as np
            df = df.assign(**{'log10(1 + mean RPM)': np.log10(1 + df['Mean RPM'])})
            chart(px.scatter(df, x='log10(1 + mean RPM)', y='log2 fold change', color='Category',
                             color_discrete_map=COLORS, hover_data=['Feature', 'Numerator reads', 'Denominator reads', 'Detection']),
                  f'{numerator} relative to {denominator}', 'ma_plot',
                  f"Positive values favor {numerator}. RPM pseudocount: {report['parameters']['pseudocount_rpm']}. No significance test; preparation and cell identity may confound differences.")
            table(df, 'fold_comparison')
        else:
            st.info('No supported features for this pair, track, and class selection.')

elif page == 'Parent explorer':
    if frame.empty:
        st.info('No parent RNAs have eligible features.')
        st.stop()
    cat = st.selectbox('RNA class', categories)
    parents = sorted(frame.loc[frame['Category'] == cat, 'Parent RNA'].unique())
    parent = st.selectbox('Parent RNA', parents)
    strands = sorted(frame.loc[frame['Parent RNA'] == parent, 'Strand'].unique())
    strand = st.radio('Reference orientation', strands, horizontal=True)
    shown_samples = st.multiselect('Samples to show', samples, default=samples)
    st.caption('Coordinates are 0-based, half-open reference intervals. On the minus strand, reference end is the read’s 5′ boundary. Molecular insert completeness is not established.')
    group = [f for f in report['features'] if f['parent_RNA'] == parent and f['strand'] == strand]
    profiles = [p for p in report.get('endpoint_profiles', []) if p['parent_RNA'] == parent and p['strand'] == strand]
    if profiles:
        rows = []
        for sample in shown_samples:
            hist = profiles[0]['histograms_by_sample'][sample]
            for label, key in [('Reference start', 'reference_start'), ('Reference end', 'reference_end')]:
                for coordinate, count in hist[key].items():
                    rows.append({'Sample': sample, 'Boundary': label, 'Position': int(coordinate), 'Reads': count})
        if rows:
            chart(px.bar(pd.DataFrame(rows).sort_values('Position'), x='Position', y='Reads', color='Sample',
                          facet_row='Boundary', barmode='group'), 'Observed start and end positions · inclusive support', 'endpoints',
                  'Endpoint profiles always include all eligible exact and accepted-inexact reads, regardless of the sidebar track.', evidence='inclusive')
    rows = []
    for f in group:
        for sample in shown_samples:
            n = f['counts'][sample][track]
            if n:
                rows.append({'Feature': f['id'], 'Sample': sample, 'Start': f['start'], 'End': f['end'], 'Reads': n})
    if rows:
        df = pd.DataFrame(rows)
        fig = go.Figure()
        plotted = df.sort_values('Reads', ascending=False).head(100)
        for sample in shown_samples:
            for _, row in plotted[plotted['Sample'] == sample].iterrows():
                fig.add_trace(go.Scatter(x=[row['Start'], row['End']], y=[row['Feature']+' · '+sample]*2, mode='lines+markers',
                                        name=sample, showlegend=False, line={'width': 4, 'color': px.colors.qualitative.Safe[samples.index(sample) % len(px.colors.qualitative.Safe)]},
                                        hovertemplate=f"{row['Feature']}<br>{sample}: {row['Reads']} reads<extra></extra>"))
        fig.update_layout(height=min(1000, max(400, len(df) * 22 + 150)))
        chart(fig, 'Supported endpoint pairs', 'intervals', 'Showing up to 100 intervals ranked by read support. Each is a representative observed interval; the table includes all selected intervals.')
        table(df, 'parent_intervals')
    else:
        st.info('No support for the selected samples and evidence track.')
    patterns = [r for r in report.get('parent_landscape', []) if r['parent_RNA'] == parent and r['strand'] == strand and r['sample'] in shown_samples]
    if patterns:
        with st.expander('Descriptive parent patterns · inclusive support'):
            table(pd.DataFrame(patterns), 'parent_patterns', evidence='inclusive')
            st.caption('Diffuse support is compatible with multiple explanations; it is not proof of degradation.')

elif page == 'Feature library':
    if frame.empty:
        st.info('No eligible features are available.')
        st.stop()
    query = st.text_input('Find a feature, parent, or RNA class')
    chosen = st.multiselect('RNA classes', categories, default=categories)
    df = frame[frame['Category'].isin(chosen)]
    if query:
        mask = df[['Feature', 'Parent RNA', 'Category']].astype(str).apply(lambda col: col.str.contains(query, case=False, regex=False)).any(axis=1)
        df = df[mask]
    st.caption(f'{len(df):,} matching features · support is read abundance, not independent molecules.')
    table(df, 'feature_library')
    if not df.empty:
        selected = st.selectbox('Inspect a feature', df['Feature'].tolist())
        f = next(f for f in report['features'] if f['id'] == selected)
        st.subheader(f"{f['category']} · {f['parent_RNA']}")
        st.caption(f"{f['start']}–{f['end']} · strand {f['strand']} · sequenced-read boundaries")
        detail = pd.DataFrame([{'Sample': s, **f['counts'][s], 'Overlapping control reads': f.get('control_overlap_screen', {}).get(s, {}).get('count')}
                               for s in samples])
        table(detail, 'feature_support')
        st.caption('Overlapping control reads can belong to a different fragment species. Missing control-screen values mean not applicable, not zero.')
        st.code(f.get('observed_sequence', ''), language=None)
        metrics = f.get('sequence_metrics', {})
        if metrics:
            st.caption(f"Length: {metrics.get('length_nt', 'not reported')} nt · GC: {metrics.get('gc_fraction', 0):.1%} · Longest homopolymer: {metrics.get('longest_homopolymer', 'not reported')} nt")
        st.caption('Representative observed RNA sequence. Parent function is not automatically fragment function.')
        folds = [('Reference', f.get('reference_sequence', ''), f.get('reference_fold', {}))]
        folds += [(f"Observed variant {v['read_id']}", v['sequence'], v['fold']) for v in f.get('observed_variant_folds', [])]
        for label, sequence, fold in folds:
            with st.expander(label + ' · predicted structure'):
                if fold.get('status') != 'ok':
                    st.info('Not available: ' + fold.get('status', 'not reported'))
                else:
                    st.write(f"MFE {fold['mfe_kcal_mol']:.2f} kcal/mol · {len(sequence)} nt")
                    st.code(sequence + '\n' + fold['structure'], language=None)
                    structure = fold['structure']; stack=[]; pairs=[]
                    for i, ch in enumerate(structure):
                        if ch == '(': stack.append(i)
                        elif ch == ')' and stack: pairs.append((stack.pop(), i))
                    fig = go.Figure()
                    import numpy as np
                    for left, right in pairs:
                        theta=np.linspace(0, np.pi, 30); mid=(left+right)/2; radius=(right-left)/2
                        fig.add_trace(go.Scatter(x=mid+radius*np.cos(theta), y=radius*np.sin(theta), mode='lines',
                                                line={'color': '#487C97', 'width': 1}, showlegend=False, hoverinfo='skip'))
                    fig.add_trace(go.Scatter(x=list(range(len(sequence))), y=[0]*len(sequence), mode='markers', text=list(sequence),
                                            marker={'size': 4, 'color': '#C07842'}, hovertemplate='Position %{x}: %{text}<extra></extra>'))
                    fig.update_xaxes(title='Sequence position (0-based)'); fig.update_yaxes(visible=False)
                    chart(fig, label + ' · predicted base pairs', 'fold_' + hashlib.sha256(label.encode()).hexdigest()[:8],
                          'Canonical-base prediction; not measured structure, glycosylation, binding, or accessibility. Raw MFE is not a cross-length function score.', evidence='structural prediction')
