# csRNA analysis pipeline

## 1. Purpose and functionality

The pipeline maps RNA sequencing reads to an annotated RNA reference, summarizes RNA classes and parent RNAs, groups observed read endpoints into features, compares support across samples, and predicts structures for selected features. The optional RBC extension investigates mapping failures using local/end-to-end alignment and sequence complexity.

Results are descriptive: read endpoints do not establish native RNA fragment boundaries, surface localization, contamination, or biological function. The paired main workflow analyzes R1 and validates R2; the RBC extension analyzes both mates independently.

## 2. Run the pipeline

Clone the repository and enter it:

```bash
git clone https://github.com/soham-shukla/wang-lab.git
cd wang-lab
```

Use Python 3.10+ with `bowtie2`, `bowtie2-build`, and `RNAfold`. With Conda already installed, create an environment once:

```bash
conda create -n csrna --override-channels -c conda-forge -c bioconda \
  python=3.12 bowtie2=2.5.5 viennarna=2.7.2
conda activate csrna
```

The following commands use macOS/Linux shell syntax. Set paths for your machine; keep the same terminal open:

```bash
CSRNA_CODE="$PWD/reconstructed_pipeline"
CSRNA_DATA="$HOME/sequencing_data"
CSRNA_MANIFEST="$CSRNA_DATA/samples_recent_batch.csv"
CSRNA_REFERENCE="$CSRNA_CODE/references/mouse_combined_with_rRNA_clean_dedup.fa"
mkdir -p "$CSRNA_DATA/results"
```

Supply the combined reference FASTA separately. Copy [the manifest template](reconstructed_pipeline/samples_recent_batch.csv) beside your FASTQs if needed, preserving existing metadata. Required columns are `sample_id,r1,r2,condition,replicate`; relative read paths resolve beside the manifest. Use extracted `.fastq.gz` files, not the sequencing archive or MultiQC HTML. Keep each biological replicate separate and use condition `RBC` for the optional extension.

Verify tools, then run a small biological pilot:

```bash
python "$CSRNA_CODE/csrna_pipeline.py" --self-test
python "$CSRNA_CODE/csrna_pipeline.py" --self-test-tools

python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_MANIFEST" --pair-handling r1-only \
  --reference "$CSRNA_REFERENCE" --analysis-mode exploratory \
  --input-scope selected_subset --pilot-pairs 1000 --sampling-seed 20261007 \
  --report-detail summary --quality-detail counts-only \
  --top-fragments 1 --max-fold-variants 1 --threads 4 \
  --output "$CSRNA_DATA/results/pilot.json" \
  > "$CSRNA_DATA/results/pilot.log" 2>&1
echo "Exit code: $?"
```

Proceed only after successful tests and exit code `0`. Sampling selects pairs across each library but still validates all input pairs. Review mapping and replicate QC before increasing `--pilot-pairs`; use new output filenames. Folding budgets limit structure predictions, not feature counts.

For large full-library runs, count once on persistent storage, then test alignment capacity:

```bash
CSRNA_CACHE="$CSRNA_DATA/count_cache/full_batch"

python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_MANIFEST" --pair-handling r1-only \
  --reference "$CSRNA_REFERENCE" --analysis-mode exploratory \
  --input-scope full_library --pilot-pairs 0 \
  --count-cache "$CSRNA_CACHE" --quality-detail counts-only --count-only \
  --output "$CSRNA_DATA/results/counting_complete.json"

python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_MANIFEST" --pair-handling r1-only \
  --reference "$CSRNA_REFERENCE" --analysis-mode exploratory \
  --input-scope full_library --pilot-pairs 0 \
  --count-cache "$CSRNA_CACHE" --quality-detail counts-only \
  --alignment-capacity-queries 10000 --threads 4 \
  --output "$CSRNA_DATA/results/alignment_capacity.json"
```

Repeat an interrupted counting command to resume. Keep input files and quality settings unchanged; never run two writers against the same cache. Retain the cache after downstream failures.

After measuring adequate RAM, disk, and runtime, run full analysis:

```bash
python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_MANIFEST" --pair-handling r1-only \
  --reference "$CSRNA_REFERENCE" --analysis-mode exploratory \
  --input-scope full_library --pilot-pairs 0 \
  --count-cache "$CSRNA_CACHE" --quality-detail counts-only \
  --report-detail summary --top-fragments 1 --max-fold-variants 1 --threads 4 \
  --output "$CSRNA_DATA/results/report.json" \
  > "$CSRNA_DATA/results/full_analysis.log" 2>&1
echo "Exit code: $?"
```

Counting is disk-backed; downstream analysis is not. Default memory guards may stop full analysis. Raise `--max-distinct-sequences` and `--max-total-distinct-sequences` only after resource assessment. Counting/capacity JSON files are not biological reports. Use `--help` for additional options.

## 3. Run the visualization app

The dashboard is read-only and does not require alignment tools or raw reads. On the viewing computer, use Python 3.11/3.12 and run from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r reconstructed_pipeline/csrna_visualizer/requirements.txt
python reconstructed_pipeline/csrna_visualizer/run.py
```

Open **http://127.0.0.1:8501/**. Leave the terminal running; press **Ctrl+C** to stop. To restart later, activate `.venv` and repeat the last command. Add `--port 8502` for a separate instance.

## 4. Load, navigate, and share results

In **Open a pipeline report**, enter the absolute path to a completed main-analysis JSON in **Local JSON path**, then press Enter. Remove any selected upload first. Alternatively, upload a small report. For RBC output, select **RBC diagnostics** and use its separate loader.

| Page | What to inspect |
|---|---|
| Overview | Scope, read totals, and analysis warnings |
| RNA composition | RNA-class abundance and unmapped reads |
| Replicate QC | Sample flags, profiles, and recurrence |
| Feature comparisons | Descriptive support differences where enabled |
| Parent explorer | Parent-RNA support and endpoint patterns |
| Feature library | Individual features and available structure evidence |
| RBC diagnostics | Mate-specific alignment modes, clipping, complexity, and origin-screen evidence |

Use page filters and exact/inexact/inclusive support tracks; download available CSV tables or chart exports. Subset reports describe only sampled reads. Untested origin screens are not negative findings, and predicted folds are not functional validation.

| Location/artifact | Contents and use |
|---|---|
| GitHub | Code, dependencies, and documentation |
| Analysis server | Raw FASTQs, reference, logs, count cache, and full outputs |
| Main output JSON | Completed main analysis; load in the main dashboard |
| RBC output directory | `rbc_diagnostics.json` for display; evidence TSVs, SQLite database, and logs for detailed review/resume |
| Box/shared storage | Downloadable copies of reports, labeled with run scope and code version |
| Slack | Links to code/results and instructions; not required to run the app |

Download shared reports completely before opening them. A local path refers to the machine **running Streamlit**: a deployed website cannot read `/Users/...` on your Mac. Large JSON files may exceed cloud or local RAM even when upload limits permit them. Prefer local path loading for large reports. A compact cloud-export format is not currently implemented.

## 5. Run optional RBC analysis

Run after a successful main pilot or full analysis, using the same reference and original paired FASTQs. Activate the analysis environment and set the paths from section 2. Keep the RBC script beside the main script.

```bash
python "$CSRNA_CODE/csrna_rbc_diagnostics.py" \
  --analysis-report "$CSRNA_DATA/results/pilot.json" \
  --manifest "$CSRNA_MANIFEST" --reference "$CSRNA_REFERENCE" \
  --output-dir "$CSRNA_DATA/results/rbc_pilot" \
  --index-cache "$CSRNA_CODE/reference_indexes" \
  --pilot-pairs 10000 --sampling-seed 20261007 --threads 4 \
  > "$CSRNA_DATA/results/rbc_pilot.log" 2>&1
echo "Exit code: $?"
```

Replace `pilot.json` with your completed main report as needed. Load `results/rbc_pilot/rbc_diagnostics.json` on the dashboard's **RBC diagnostics** page. Local matches do not alter main-pipeline counts; do not count the two mates as independent molecules.

Repeat the identical command to reuse completed stages. Changed inputs, code, or settings require a new output directory. After resource assessment, use `--pilot-pairs 0` with a new directory for full RBC libraries. See [RBC diagnostic details](reconstructed_pipeline/RBC_DIAGNOSTICS.md) for optional origin-reference screens and recovery instructions.

## 6. Report errors or bugs

Open a [GitHub issue](https://github.com/soham-shukla/wang-lab/issues) with:

- The command or dashboard steps, expected result, and actual result.
- Complete error/traceback and relevant log lines, not only the final health-check failure.
- Code commit (`git rev-parse HEAD`), OS, Python/tool versions, and local versus cloud deployment.
- Report schema/scope, file size, sample identifiers, and a small reproducible example where possible.

Do not post private sequencing data, credentials, or large reports publicly. Share approved evidence through restricted storage. For diagnostics, run from the repository root in the appropriate environment:

```bash
python reconstructed_pipeline/csrna_pipeline.py --self-test
python reconstructed_pipeline/csrna_pipeline.py --self-test-tools
CSRNA_TEST_TOOLS=1 python reconstructed_pipeline/test_count_cache.py
CSRNA_TEST_TOOLS=1 python reconstructed_pipeline/test_rbc_diagnostics.py
python reconstructed_pipeline/csrna_visualizer/test_visualizer.py
```

Passing tests cover the tested cases; they do not guarantee full-library capacity or biological validity.
