# csRNA pipeline data

A local, read-only interactive companion to `reconstructed_pipeline/csrna_pipeline.py`. It reads the pipeline JSON output; it never imports, executes, or edits the analysis pipeline. No FASTQ files are required.

## Start the app

From the project root (`/Users/soham/Documents/wang-lab`), using Python 3.10 or newer:

```bash
python -m pip install -r reconstructed_pipeline/csrna_visualizer/requirements.txt
python reconstructed_pipeline/csrna_visualizer/run.py --report /path/to/csrna_results.json
```

On the current workstation, `/opt/anaconda3/bin/python` already has the dependencies installed:

```bash
/opt/anaconda3/bin/python reconstructed_pipeline/csrna_visualizer/run.py --report /path/to/csrna_results.json
```

If your terminal is already in `reconstructed_pipeline`, use `python csrna_visualizer/run.py` instead. The launcher locates the app relative to its own file; report paths are resolved relative to your terminal’s working directory. After moving the folder, stop any server started from the old location with Ctrl+C and relaunch it from the new location.

Open **http://127.0.0.1:8501**. The launcher selects the light scientific theme, binds only to localhost, and disables Streamlit usage telemetry. Stop it with Ctrl+C. Use `--port 8502` if needed. Omit `--report` to choose a file in the sidebar.

Large pipeline reports should be opened by local path: the normal upload limit is 200 MB. The app has been exercised with the existing approximately 311 MB full-reference test report. It parses the original JSON once, validates its count matrices, then caches a smaller presentation model without unused whole-reference metadata or individual read alignments. Initial parsing still needs memory proportional to the report; this is not a streaming reader for arbitrarily large reports. A changed file is reloaded using its modification time and size. Uploaded files take precedence over the path; clear the upload to return to the path.

## Explore

| View | Use it to answer |
|---|---|
| Overview | How much input remains unmapped? How much exact/inexact evidence supports features? Which diagnostic screens were actually run? |
| RNA composition | Which RNA classes dominate reads? Does high abundance reflect many features or a few dominant parents? |
| Replicate QC | Do library profiles agree within conditions, which features recur, and do preparation differences accompany outliers? |
| Feature comparisons | Which shared features have support in each sample? For full-library reports, which are relatively elevated or depleted? |
| Parent explorer | Where do reads begin and end on a selected parent? Which endpoint pairs have support? |
| Feature library | Find an individual feature, inspect its counts, overlapping control evidence, sequence, and available structure predictions. |

The sidebar evidence track changes feature-based charts and tables. Inclusive support is exact plus accepted-inexact support. Composition and unmapped charts describe all input mapping categories, not just feature-eligible reads. Endpoint histograms and parent-pattern summaries are explicitly inclusive. Four-component support charts distinguish raw-unique and equivalence-resolved evidence.

## Interpretation

- **Scope is always visible.** Selected-subset reports retain tables and descriptive heatmaps but cannot produce biological fold-change rankings. They do not estimate the entire library automatically.
- Composition defaults to **all supplied reads**, including unmapped material. Mapped-only percentages are conditional on the mapped portion. Hiding rRNA does not renormalize the denominator.
- Feature diversity is observed feature count, affected by depth and eligibility. Read abundance does not establish independent molecule counts.
- Heatmaps use `log10(1 + RPM)` with one scale across all rows. They show the top selected features by summed RPM, not row-standardized enrichment. Zero means no observed support in the selected track.
- Full-library pairwise comparisons use the report’s denominator and RPM pseudocount. Both-zero features are omitted from the comparison, and control-only features remain visible. Ratios are descriptive, not significance tests or proof of surface specificity.
- Parent intervals are **observed read boundaries**. A sharp peak does not establish complete molecular ends. The interval chart displays up to 100 supported intervals; the accompanying table includes all selected intervals.
- Overlapping control reads may be a different fragment species. Missing control-screen values are not zero counts.
- Predicted folds are shown separately for the reference and observed variants. Arc plots show predicted base pairs, not measured extracellular structure, glycosylation, receptor binding, or function. Unfolded entries remain explicitly unavailable.
- RNA classes depend on reference annotation. Missing classes are not proof of biological absence. Sample/preparation differences and unconfirmed replication, when applicable, limit inference.

## Export

- **Download table · CSV** exports the displayed table with scope, evidence-track, and normalization labels. Use this button rather than the data grid’s built-in CSV shortcut when you need those labels. Text cells that could be spreadsheet formulas are prefixed with an apostrophe.
- Hover over a chart and use its camera icon to save **SVG**. Charts include scope and normalization annotations.
- **Prepare interactive HTML export** creates a downloadable, self-contained chart with embedded Plotly JavaScript. It can be opened offline. Exported artifacts contain the displayed experimental data; share them deliberately.
- Plotly supports zooming, hover details, and toggling individual legend entries. Tables support sorting and search.

## Compatibility and tests

The adapter accepts schema 2 pipeline reports, including reports preceding the Goal 2 addition. It checks sample and feature identities, matrix axes, integer count conservation, and denominators. Optional evidence sections are shown only when supplied. Counts and normalized values are calculated directly from the existing report; the pipeline remains the source of truth.

```bash
python reconstructed_pipeline/csrna_visualizer/test_visualizer.py
```

Tests cover normalization, rRNA visibility, count conservation, reciprocal comparisons, control-only support, subset restrictions, invalid reports, empty results and filters, CSV safety, chart export generation, and navigation through all six views. Real-report checks also exercise the existing 83-feature rRNA report and the 574-feature, six-category raw-read test report. These tests establish checked behavior; they do not prove that every possible input is bug-free.

Files: `run.py` is the local launcher; `app.py` is the interface; `model.py` is the read-only report adapter; `test_visualizer.py` contains isolated tests. The three analysis deliverables are not modified by this application.

Validation on 2026-10-04: 11 automated tests passed; all five views passed with the real 83-feature and 574-feature reports. Browser review caught and fixed a Plotly-version stacking issue; a regression test now checks shared offsets for stacked bars. A browser-downloaded composition CSV was independently verified to preserve all 3,000 diagnostic-subset reads, 100% per-sample composition totals, and scope labels. The preview uses a selected-subset report, not complete-library biological results.

## Replicate batch support (2026-10-07)

Schema 2.1 exploratory reports and paired-manifest R1 reports are supported alongside legacy reports. R1-only and no-control policies are displayed explicitly. See [pipeline run instructions](../README.md#replicate-libraries-and-paired-file-input--2026-10-07) and the [15-library manifest](../samples_recent_batch.csv).

Replicate QC preserves individual library columns. It provides condition-faceted composition, Pearson correlations and PCA on log2(1 + RPM) shared-feature profiles, recurrence tables, preparation metadata, paired-input audit, and integer raw-count export. Constant correlations are blank and identical profiles have no PCA; zero selections are handled explicitly. The RNA-class filter affects feature QC, not whole-read composition or the raw-count export. Unknown biological/technical status is never inferred. Neither PCR cycle correction nor inferential statistics is performed.

## Accuracy/performance audit (2026-10-08)

The current 15-library manifest records confirmed independent biological samples and supported PCR-cycle metadata; other preparation details remain unknown. Existing reports retain the metadata recorded at analysis time. New reports expose measured stage durations under **Report details**, excluding final assembly/serialization. Identical normalized profiles correctly have no PCA ordination. Sixteen frontend tests and navigation through all six views with the fresh 574-feature diagnostic report passed. See the [audit report](../PERFORMANCE_ACCURACY_AUDIT.md) for benchmarks and remaining full-library/R1-only limitations.
