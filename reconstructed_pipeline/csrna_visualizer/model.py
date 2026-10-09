"""Read-only adapter for csRNA pipeline schema 2 reports. No pipeline imports."""
import json
import math
from pathlib import Path
import pandas as pd

TRACKS = ('exact', 'inexact', 'inclusive')


def _validate(report):
    if not isinstance(report, dict) or not str(report.get('schema_version', '')).startswith('2.'):
        raise ValueError('Choose a csRNA pipeline schema 2 JSON report.')
    samples = report.get('samples')
    features = report.get('features')
    if not isinstance(samples, dict) or not samples or not isinstance(features, list) or not all(isinstance(f, dict) for f in features):
        raise ValueError('Report must contain samples and a feature list.')
    metadata = report.get('sample_metadata', [])
    if metadata:
        names = [row['sample_id'] for row in metadata]
        if len(names) != len(set(names)) or set(names) != set(samples):
            raise ValueError('Sample metadata must match every sample exactly once.')
    scope = report.get('parameters', {}).get('input_scope')
    if scope not in ('full_library', 'selected_subset'):
        raise ValueError('Report has no recognized input scope.')
    pc = report['parameters'].get('pseudocount_rpm')
    if scope == 'full_library' and (not isinstance(pc, (int, float)) or not math.isfinite(pc) or pc <= 0):
        raise ValueError('Full-library report needs a finite positive RPM pseudocount.')
    for name, row in samples.items():
        total = row.get('total_input_reads')
        denom = row.get('rpm_denominator')
        if type(total) is not int or total <= 0 or not isinstance(denom, (int, float)) or not math.isfinite(denom) or denom <= 0:
            raise ValueError(f'Invalid read totals or normalization denominator for {name}.')
        categories = row.get('category_counts', {})
        if type(row.get('unmapped_reads')) is not int or row['unmapped_reads'] != categories.get('unmapped', 0):
            raise ValueError(f'Unmapped counts disagree for {name}.')
        if row.get('mapped_reads') != total - row['unmapped_reads'] - categories.get('below_min_read_length', 0):
            raise ValueError(f'Mapping totals disagree for {name}.')
        if any(type(n) is not int or n < 0 for n in categories.values()) or sum(categories.values()) != total:
            raise ValueError(f'Category counts do not conserve input reads for {name}.')
    ids = [f.get('id') for f in features]
    if any(not isinstance(i, str) for i in ids) or len(set(ids)) != len(ids):
        raise ValueError('Feature identifiers are invalid or duplicated.')
    for f in features:
        if not {'parent_RNA', 'category', 'strand', 'start', 'end', 'counts'} <= f.keys():
            raise ValueError('Feature annotations are incomplete.')
        for sample in samples:
            counts = f['counts'].get(sample, {})
            if any(type(counts.get(t)) is not int or counts[t] < 0 for t in TRACKS):
                raise ValueError('Feature counts must be nonnegative integers for every sample and track.')
            if counts['exact'] + counts['inexact'] != counts['inclusive']:
                raise ValueError('Exact + inexact support does not equal inclusive support.')
    for sample, summary in samples.items():
        if sum(f['counts'][sample]['inclusive'] for f in features) > summary['mapped_reads']:
            raise ValueError(f'Eligible feature support exceeds mapped reads for {sample}.')
    for track in TRACKS:
        m = report.get('count_matrices', {}).get(track, {})
        if m.get('feature_ids') != ids or m.get('sample_ids') != list(samples):
            raise ValueError('Count-matrix axes do not match features and samples.')
        if m.get('data') != [[f['counts'][s][track] for s in samples] for f in features]:
            raise ValueError('Count matrices disagree with feature counts.')
    return report


def validate(report):
    try:
        return _validate(report)
    except (AttributeError, KeyError, TypeError, IndexError) as error:
        raise ValueError('Report structure is incomplete or malformed; choose an unmodified pipeline JSON report.') from error


def compact(report):
    """Drop unused whole-reference metadata before caching; preserve source evidence."""
    keep = ('schema_version', 'created_utc', 'parameters', 'target', 'controls', 'samples', 'sample_metadata',
            'features', 'count_matrices', 'endpoint_profiles', 'parent_landscape', 'warnings', 'qc_alerts',
            'unmapped_qc', 'landscape_comparison', 'support_reporting', 'paired_input_qc', 'performance')
    out = {k: report[k] for k in keep if k in report}
    parents = {f['parent_RNA'] for f in report['features']}
    ref = report.get('reference', {})
    out['reference'] = {'sha256': ref.get('sha256'), 'records': {p: r for p, r in ref.get('records', {}).items() if p in parents}}
    return out


def load_path(path):
    p = Path(path).expanduser()
    if not p.is_file():
        raise ValueError('Report path does not point to a file.')
    with p.open(encoding='utf-8') as h:
        if not h.read(512).lstrip().startswith('{'):
            raise ValueError('Choose a pipeline JSON report, not a FASTQ, FASTA, or other input file.')
        h.seek(0)
        return compact(validate(json.load(h)))


def feature_frame(report, track='inclusive'):
    rows = []
    for f in report['features']:
        row = {'Feature': f['id'], 'Parent RNA': f['parent_RNA'], 'Category': f['category'], 'Strand': f['strand'],
               'Start': f['start'], 'End': f['end'], 'Length (nt)': len(f.get('observed_sequence', '')) or f['end'] - f['start']}
        for s, summary in report['samples'].items():
            row[f'{s} reads'] = f['counts'][s][track]
            row[f'{s} RPM'] = f['counts'][s][track] * 1e6 / summary['rpm_denominator']
        rows.append(row)
    columns = ['Feature', 'Parent RNA', 'Category', 'Strand', 'Start', 'End', 'Length (nt)'] + [v for s in report['samples'] for v in (f'{s} reads', f'{s} RPM')]
    return pd.DataFrame(rows, columns=columns)


def composition(report, denominator='All supplied reads', omit_rrna=False):
    rows = []
    for s, summary in report['samples'].items():
        counts = dict(summary['category_counts'])
        if denominator == 'Mapped reads only':
            counts = {k: v for k, v in counts.items() if k not in ('unmapped', 'below_min_read_length')}
        if omit_rrna:
            counts.pop('rRNA', None)
        # Excluding rRNA changes visibility, never the selected denominator.
        denom = summary['total_input_reads'] if denominator == 'All supplied reads' else summary['mapped_reads']
        for category, n in sorted(counts.items()):
            rows.append({'Sample': s, 'Category': category, 'Reads': n, 'Percent': n * 100 / denom if denom else None})
    return pd.DataFrame(rows, columns=['Sample', 'Category', 'Reads', 'Percent'])


def comparison(report, numerator, denominator, track):
    if report['parameters']['input_scope'] != 'full_library':
        raise ValueError('Biological fold comparisons are disabled for selected-subset reports.')
    if numerator == denominator or numerator not in report['samples'] or denominator not in report['samples']:
        raise ValueError('Choose two different samples.')
    pc = report['parameters'].get('pseudocount_rpm')
    if not isinstance(pc, (int, float)) or not math.isfinite(pc) or pc <= 0:
        raise ValueError('Report lacks a valid RPM pseudocount.')
    rows = []
    for f in report['features']:
        a, b = f['counts'][numerator][track], f['counts'][denominator][track]
        if not (a or b):
            continue
        ar = a * 1e6 / report['samples'][numerator]['rpm_denominator']
        br = b * 1e6 / report['samples'][denominator]['rpm_denominator']
        rows.append({'Feature': f['id'], 'Parent RNA': f['parent_RNA'], 'Category': f['category'],
                     'Numerator reads': a, 'Denominator reads': b, 'Numerator RPM': ar, 'Denominator RPM': br,
                     'Mean RPM': (ar + br) / 2, 'log2 fold change': math.log2((ar + pc) / (br + pc)),
                     'Detection': 'Both samples' if a and b else 'Numerator only' if a else 'Denominator only'})
    cols = ['Feature', 'Parent RNA', 'Category', 'Numerator reads', 'Denominator reads', 'Numerator RPM', 'Denominator RPM', 'Mean RPM', 'log2 fold change', 'Detection']
    return pd.DataFrame(rows, columns=cols).sort_values(['log2 fold change', 'Feature'], ascending=[False, True])


def provenance(report, track):
    return (f"Scope: {report['parameters']['input_scope']} · Track: {track} · "
            f"RPM denominator: {report['parameters'].get('denominator', 'unspecified')} · "
            f"Samples: {', '.join(report['samples'])} · Observed read endpoints; molecular ends unverified")


def csv_bytes(frame):
    safe = frame.copy()
    for col in safe.columns:
        safe[col] = safe[col].map(lambda x: "'" + x if isinstance(x, str) and x.lstrip().startswith(('=', '+', '-', '@')) else x)
    return safe.to_csv(index=False).encode('utf-8-sig')


def replicate_summary(report, track='inclusive', selected=None, categories=None):
    """Descriptive QC only: no pooling, tests, inferred replicates, or PCR correction."""
    import numpy as np
    samples = list(report['samples']) if selected is None else list(selected)
    if len(set(samples)) != len(samples) or not set(samples) <= set(report['samples']):
        raise ValueError('Select distinct report samples')
    metadata = {r['sample_id']: r for r in report.get('sample_metadata', [])}
    rows = []
    for sample in samples:
        row = dict(metadata.get(sample, {}))
        row['sample_id'] = sample
        for key in ('condition', 'replicate_type', 'pcr_cycles', 'library_preparation'):
            row[key] = row.get(key) or 'unknown'
        rows.append(row)
    meta = pd.DataFrame(rows, columns=list(dict.fromkeys(['sample_id', 'condition', 'replicate_type', 'pcr_cycles', 'library_preparation'] + [k for r in rows for k in r])))
    features = [f for f in report['features'] if categories is None or f['category'] in categories]
    raw = np.array([[f['counts'][s][track] for f in features] for s in samples], dtype=float).reshape(len(samples), len(features))
    denominators = np.array([report['samples'][s]['rpm_denominator'] for s in samples])
    transformed = np.log1p(raw * 1e6 / denominators[:, None]) / np.log(2) if samples else raw
    # Drop globally absent features. Constant profiles have undefined correlation.
    transformed = transformed[:, np.any(raw > 0, axis=0)]
    correlation = np.full((len(samples), len(samples)), np.nan)
    if transformed.shape[1] >= 2:
        valid = np.flatnonzero(np.ptp(transformed, axis=1) > 0)
        if len(valid):
            correlation[np.ix_(valid, valid)] = np.atleast_2d(np.corrcoef(transformed[valid])) if len(valid) > 1 else [[1.0]]
    scores, variance = pd.DataFrame(), []
    if len(samples) >= 2 and transformed.shape[1]:
        centered = transformed - transformed.mean(axis=0)
        if np.any(transformed != transformed[0]):
            u, singular, _ = np.linalg.svd(centered, full_matrices=False)
            explained = singular ** 2 / (singular ** 2).sum()
            values = np.zeros((len(samples), 2))
            n = min(2, len(singular))
            values[:, :n] = u[:, :n] * singular[:n]
            scores = meta.copy()
            scores['PC1'], scores['PC2'] = values[:, 0], values[:, 1]
            variance = list(explained[:2]) + [0.0] * (2 - n)
    recurrence = []
    for condition in meta['condition'].unique():
        members = meta.loc[meta['condition'] == condition, 'sample_id'].tolist()
        for f in features:
            counts = [f['counts'][s][track] for s in members]
            recurrence.append({'Condition': condition, 'Feature': f['id'], 'Category': f['category'],
                               'Parent RNA': f['parent_RNA'], 'Samples detected': sum(c > 0 for c in counts),
                               'Samples examined': len(members), 'Total reads': sum(counts),
                               'Mean sample RPM': sum(c * 1e6 / report['samples'][s]['rpm_denominator'] for c, s in zip(counts, members)) / len(members)})
    return meta, pd.DataFrame(correlation, index=samples, columns=samples), scores, variance, pd.DataFrame(recurrence)
