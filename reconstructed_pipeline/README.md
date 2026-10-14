# csRNA fragment discovery and landscape comparison — version 2

This version updates the existing pipeline in response to the supervisor's 09/22 feedback. The primary biological unit is an observed fragment population with supported paired endpoints, across RNA categories. Exact-match rRNA assembly is no longer the gatekeeper for discovery. A merged reference interval remains a coverage summary, not evidence of one long RNA molecule.

The three maintained deliverables are `csrna_pipeline.py`, `csrna_pipeline.ipynb`, and `README.md`. The notebook embeds the complete standalone script in ordered cells tagged `pipeline-source`. This document contains exactly the notebook's Markdown annotations. The core pipeline remains a standalone script. The optional RBC post-analysis companion described below adds its own diagnostic outputs.

## Resumable counting and alignment-capacity workflow — 2026-10-13

For the large paired batch, use **persistent counting before full alignment**. Version 2.2.0 adds `--count-cache DIRECTORY` to the main standalone script. It validates R1/R2 together and counts R1 into SQLite, committing counts, positional quality sums (when enabled), and input offsets atomically. Counting has no distinct-sequence cap. Existing memory guards apply only when loading completed counts into the downstream in-memory feature engine. Completed samples are reused; an interrupted sample resumes from the last committed batch. At most the uncommitted batch is recounted.

This currently supports **full-library paired manifests only** (`--pilot-pairs 0 --input-scope full_library`). Legacy single-end, FASTA/CSV, and uniform-pilot workflows remain unchanged. No mate merging, sampling, abundance filtering, or change to feature definitions is introduced by counting. The core script remains standalone and needs no Python package beyond the standard library for SQLite.

### A. Count and validate once

On the server, activate the environment and set paths as in the walkthrough below. Choose an empty dedicated cache directory on persistent storage with sufficient free space (not temporary job storage). Use the same directory on subsequent invocations.

```bash
CSRNA_CODE="$HOME/reconstructed_pipeline"
CSRNA_DATA="$HOME/sequencing_data"
CSRNA_MANIFEST="$CSRNA_DATA/samples_recent_batch.updated.csv"
CSRNA_REFERENCE="$CSRNA_CODE/references/mouse_combined_with_rRNA_clean_dedup.fa"
CSRNA_CACHE="$CSRNA_DATA/count_cache/full_batch"
mkdir -p "$CSRNA_DATA/results"

/usr/bin/time -v python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_MANIFEST" --pair-handling r1-only \
  --analysis-mode exploratory --input-scope full_library --pilot-pairs 0 \
  --reference "$CSRNA_REFERENCE" \
  --count-cache "$CSRNA_CACHE" --count-only \
  --quality-detail counts-only \
  --output "$CSRNA_DATA/results/counting_complete.json" \
  > "$CSRNA_DATA/results/counting.stdout.log" \
  2> "$CSRNA_DATA/results/counting.log"
echo "Exit code: $?"
```

Adjust the manifest filename if your working file has another name. `counts-only` validates FASTQ quality syntax but omits positional quality arrays. Choose `--quality-detail positional` here if those arrays are required later; all cache-reusing commands must keep the same quality policy. Changing it requires a new cache and recount. Count-only mode does not require Bowtie2/RNAfold, although the reference path is still required and hashed for provenance.

Repeat the command after an interruption. A hard kill may leave `.count.lock`: verify the previous process has stopped before removing that one stale lock. Do not delete `counts.sqlite`, `counts.sqlite-wal`, or other cache data to resolve a downstream failure. SQLite recovers committed transactions; storage failure or corruption is not automatically repaired. Do not run concurrent writers against one cache.

Every invocation verifies source SHA256 hashes; retries still read the compressed input bytes. Completed libraries skip decompression/FASTQ validation/counting. Resuming partway through `.gz` may decompress the earlier portion to seek to the stored offset; ordinary gzip cannot guarantee constant-time random access. Final report publication rechecks input hashes. Inputs, sample identities, paths, counting format/Python version, and quality policy must match; otherwise use a new cache. Annotation-only manifest changes may reuse counts. The checkpoint batch defaults to 10,000 pairs; `--count-checkpoint-records` changes checkpoint frequency, not results.

The count-stage JSON lists complete per-sample read/distinct totals and the union distinct count. It is **not a biological analysis report** and cannot be used by the dashboard or RBC extension.

### B. Test alignment capacity using the saved counts

After counting succeeds, align a reproducible uniform sample of distinct R1 sequences against the requested reference, with the same all-best-hit mapping policy as the main analysis:

```bash
/usr/bin/time -v python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_MANIFEST" --pair-handling r1-only \
  --analysis-mode exploratory --input-scope full_library --pilot-pairs 0 \
  --reference "$CSRNA_REFERENCE" \
  --count-cache "$CSRNA_CACHE" --quality-detail counts-only \
  --alignment-capacity-queries 10000 --sampling-seed 20261007 --threads 4 \
  --output "$CSRNA_DATA/results/alignment_capacity_10000.json" \
  > "$CSRNA_DATA/results/alignment_capacity_10000.stdout.log" \
  2> "$CSRNA_DATA/results/alignment_capacity_10000.log"
echo "Exit code: $?"
```

Inspect index-build/alignment seconds, mapped queries, total retained best hits, and maximum/p95 best-hit multiplicity in the JSON. Inspect process resource measurements in the log and monitor allocated-job RAM/disk on the server. Alignment measurements include index construction; indexes are rebuilt for each main-script capacity run. This is a distinct-sequence sample, not an abundance-weighted read pilot. Rare repetitive queries may be missed, and nonlinear hit volume can dominate larger runs. No full-run runtime/RAM guarantee is inferred.

If feasible, repeat with 50,000 and then 100,000 queries and new output/log names, reusing the same cache and seed. Capacity testing bypasses the materialization guards by loading only the selected queries. It does not construct endpoint features, fold variants, or test final-report memory. A failed capacity run leaves committed counts reusable. Its JSON also has a distinct diagnostic schema and is not accepted as a main or RBC biological report.

### C. Continue the main pipeline from counts

Use the normal full-library main command below, adding:

```bash
--count-cache "$CSRNA_CACHE" --quality-detail counts-only
```

Omit `--count-only` and `--alignment-capacity-queries`. Use exactly one `--quality-detail` argument matching the cache. Completed counts and paired validation are reused; alignment, global feature discovery, folding and reporting run normally. The resulting schema-2 report opens in the existing dashboard and can precede RBC diagnostics as before. New reports record cache provenance under `counting_checkpoint`.

The default 250,000-per-sample and 500,000-union guards still protect downstream memory. If exceeded, the message identifies the retained database; no counts are discarded. Increase `--max-distinct-sequences` and `--max-total-distinct-sequences` only after assessing downstream resources. This change does **not** make alignment, feature construction, or JSON serialization disk-backed. If those stages cannot fit, further engineering is needed; repeatedly raising limits is not a capacity solution.

Use new output filenames for repeat measurements. `--overwrite` is only needed when intentionally replacing a previously published JSON, not for recovering an incomplete counting stage. Keep the cache until successful analysis and backup; it is an intentional intermediate and may be large.

### Verification

```bash
CSRNA_TEST_TOOLS=1 python "$CSRNA_CODE/test_count_cache.py"
python "$CSRNA_CODE/csrna_pipeline.py" --self-test-tools
CSRNA_TEST_TOOLS=1 python "$CSRNA_CODE/test_rbc_diagnostics.py"
```

The counting tests compare positional quality sums and counts with the legacy reader, inject an interruption between checkpoints, check completed-cache reuse, input-policy invalidation, malformed mates, locks, both memory guards, unrelated-file protection, uniform query sampling, and real-tool feature/count equivalence. These fixtures test behavior, not server-scale throughput.

## Preparation documentation update — 2026-10-12

The lab's `RNA_Sequencing_Documentation.docx` reports Watchmaker mRNA Library Prep Kit used without RNA selection, 1 ng RNA input, strand-specific construction, unique dual indexes, and Fragment Analyzer assessment. The batch manifest now records this provenance. Fragmentation is **unconfirmed**, read-to-RNA orientation remains **unverified**, and adapter trimming and UMI status remain unknown. UDI sample indexes do not establish molecular identifiers. The phrase “150 nt” does not resolve cycles per mate; verify FASTQ lengths and run settings.

The document reports 24 PCR cycles, conflicting with earlier sample-specific values for several groups. `pcr_cycles` retains the earlier recorded values (RBC remains blank pending sample-specific confirmation); `pcr_cycles_document_reported` records 24 separately. `pcr_cycles_review_status` marks conflicting, consistent, or unresolved sample-specific records. Neither record is silently preferred for count correction. Fragment Analyzer assessment is recorded as a procedure, not proof that each sample passed a DNA-contamination test.

Endpoint clusters remain descriptive observed-read features. Library fragmentation and cDNA synthesis can create boundaries; even a reconstructed library insert need not be a native RNA molecule. Endpoint peaks, discrete/diffuse patterns, and predicted fragment folds do not establish native csRNA species or surface localization. Prioritize RNA-category composition, parent support, coverage, and replicate consistency while preparation details are unresolved. No strand filtering, trimming, UMI inference, automatic deduplication, or numerical PCR correction was introduced.

**Applying the update:** copy the updated manifest to the server and merge its annotation columns into your working manifest while preserving correct FASTQ paths and any newer verified metadata. New main reports include the supplied metadata and preparation warnings. Existing JSON reports remain unchanged and retain their original metadata. For RBC analysis, the current manifest supplies diagnostic annotations, while `main_analysis_sample_metadata` preserves the original main-report snapshot. Sample identities and reference/input integrity checks remain enforced. Updated code or metadata requires a new RBC output directory; existing completed reports are not rewritten. The dashboard explains endpoint limitations even for older reports and exposes current versus original RBC metadata explicitly.

## Start here: main analysis, optional RBC extension, and dashboard

This walkthrough uses the current 15-library paired FASTQ batch. Run the main analysis first; RBC diagnostics is an **optional, separate post-analysis step**. Both reports open in the same read-only Streamlit application, using different pages. The browser does not run alignment or change counts.

1. [Prepare files and environment](#1-prepare-files-and-environment).
2. [Validate and run the main pilot](#2-validate-and-run-the-main-pilot).
3. [Decide whether to run the complete main libraries](#3-decide-whether-to-run-the-complete-main-libraries).
4. [Optionally run RBC diagnostics](#4-optionally-run-rbc-diagnostics-after-main-analysis).
5. [Open the reports in the web application](#5-open-the-reports-in-the-web-application).
6. [Explore, interpret, and export](#6-explore-interpret-and-export-results).
7. [Troubleshoot and rerun](#7-troubleshoot-and-rerun).

The dated audits later in this document describe earlier datasets and versions. Their historical sample counts, test counts, replication status, and report-size measurements are not findings from a newly submitted batch. For this batch, the 15 samples are confirmed independent biological samples; preparation, UMI, and native-end completeness remain unknown.

### 1. Prepare files and environment

On the Linux analysis server, the examples assume:

```text
~/reconstructed_pipeline/
    csrna_pipeline.py
    csrna_rbc_diagnostics.py
    test_rbc_diagnostics.py
    samples_recent_batch.csv
    csrna_visualizer/              # needed only where you run the dashboard
    references/
        mouse_combined_with_rRNA_clean_dedup.fa
~/sequencing_data/
    <the 30 original R1/R2 .fastq.gz files>
    samples_recent_batch.csv       # working copy of the manifest
    results/
```

Keep both analysis scripts beside each other: the RBC script imports input-reading helpers from the main script. Supply the actual combined reference FASTA; it is not downloaded or reconstructed automatically. The same reference must be used for main and RBC analysis. If your reference lives elsewhere, change `CSRNA_REFERENCE` below. Do not use the MultiQC HTML or the `.tar` archive as a read input; use the extracted FASTQs without decompressing the `.gz` files.

Activate your existing working environment. For the previously configured home-directory installation:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate "$HOME/envs/csrna"
```

Adjust the installation prefix if yours differs. If Conda is available but this environment does not yet exist, create it once, then activate it:

```bash
conda create --prefix "$HOME/envs/csrna" --override-channels \
  -c conda-forge -c bioconda \
  python=3.12 bowtie2=2.5.5 viennarna=2.7.2
conda activate "$HOME/envs/csrna"
```

Conda itself must already be available; an environment-creation command does not install Conda. If it is missing, use your institution's approved user-level installation or software modules. No system Python upgrade is needed. The main script needs Python >=3.10, `bowtie2`, `bowtie2-build`, and `RNAfold`; the RBC extension does not call RNAfold. Jupyter is optional.

In the same server terminal, define paths and prepare the working directory:

```bash
CSRNA_CODE="$HOME/reconstructed_pipeline"
CSRNA_DATA="$HOME/sequencing_data"
CSRNA_REFERENCE="$CSRNA_CODE/references/mouse_combined_with_rRNA_clean_dedup.fa"
cd "$CSRNA_DATA"
mkdir -p results
# Copy only if absent; preserve any metadata edits in an existing manifest.
if [ ! -e samples_recent_batch.csv ]; then
  cp "$CSRNA_CODE/samples_recent_batch.csv" samples_recent_batch.csv
fi
python --version
command -v bowtie2 bowtie2-build RNAfold
```

Review `samples_recent_batch.csv` before running. Required columns are `sample_id,r1,r2,condition,replicate`; relative mate paths resolve beside the manifest. The supplied file lists all 15 samples separately. Keep conditions `BMDC`, `BMDC_CpG`, `BMDM`, `M1`, and `RBC`; the extension selects **exactly `RBC`**. Retain confirmed biological replication and known PCR metadata, and leave unknown preparation details unknown. Sample indexes are not UMIs. Do not invent an adapter sequence or apply a guessed trimming procedure.

For an unrelated single-end dataset, use repeated `--sample NAME=/path/reads.fastq.gz` instead of `--manifest` and `--pair-handling`. See [Run the updated pipeline](#run-the-updated-pipeline) and [Full-library inputs and quality evidence](#full-library-inputs-and-quality-evidence) for FASTA/CSV rules. The RBC extension requires paired FASTQs and a paired manifest.

### 2. Validate and run the main pilot

Run the software checks before processing the batch:

```bash
python "$CSRNA_CODE/csrna_pipeline.py" --self-test
python "$CSRNA_CODE/csrna_pipeline.py" --self-test-tools
```

Proceed only when both commands succeed. These checks validate the implementation and installed tools; they do not validate the new biological libraries. Use the institution's compute-job allocation for long runs, with at least as many allocated CPUs as `--threads`.

Start with a reproducible 1,000-pair pilot per library:

```bash
/usr/bin/time -v python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_DATA/samples_recent_batch.csv" \
  --pair-handling r1-only \
  --analysis-mode exploratory \
  --input-scope selected_subset \
  --pilot-pairs 1000 --sampling-seed 20261007 \
  --reference "$CSRNA_REFERENCE" \
  --report-detail summary --quality-detail counts-only \
  --top-fragments 1 --max-fold-variants 1 \
  --threads 4 \
  --output "$CSRNA_DATA/results/recent_batch_pilot.json" \
  > "$CSRNA_DATA/results/recent_batch_pilot.stdout.log" \
  2> "$CSRNA_DATA/results/recent_batch_pilot.log"
echo "Exit code: $?"
```

`/usr/bin/time -v` is for Linux. On macOS, omit that prefix or use `/usr/bin/time -l`. The command validates **every R1/R2 pair**, then maps the sampled R1 observations; R2 is validated but not mapped by the main pipeline. Sampling reduces alignment work, not whole-file validation/checksum time. These are R1 read counts, not paired-insert or molecule counts. No controls were designated for this batch, so use `exploratory` mode rather than treating another cell type as a negative control.

The pilot produces up to 15,000 R1 observations. `summary` preserves feature counts and QC while omitting individual alignments and unused reference metadata. `counts-only` still checks FASTQ syntax but omits positional quality summaries. Small folding budgets limit structure predictions, not mapping or feature discovery. Use `--quality-detail positional`, `--report-detail complete`, or larger folding budgets only when those outputs are needed and memory permits.

Check the exit code immediately: **0 means success**. Inspect warnings and resource use:

```bash
tail -n 60 "$CSRNA_DATA/results/recent_batch_pilot.log"
python - "$CSRNA_DATA/results/recent_batch_pilot.json" <<'PY'
import json, sys
with open(sys.argv[1]) as handle:
    report = json.load(handle)
assert report.get('status') == 'complete', 'Report is not complete'
print('Complete report; scope:', report['parameters']['input_scope'])
print('Samples:', ', '.join(m['sample_id'] for m in report['sample_metadata']))
PY
```

Do not use an old output as evidence that a failed rerun succeeded. Review **Overview** and **Replicate QC** in the dashboard before scaling up. An existing successful main pilot is sufficient to start the optional RBC extension; completing the entire main dataset first is not required.

### 3. Decide whether to run the complete main libraries

For libraries that exceed counting limits, use the resumable counting and capacity workflow above first. The command below without `--count-cache` uses legacy in-memory counting.

First repeat the pilot with `--pilot-pairs 10000` and new output/log filenames. Compare runtime, maximum resident memory, mapping rates, and observed feature counts. Sampling variation is expected. Full-library capacity has not been established for this batch: the main pipeline stores distinct sequences and analysis structures in memory, and defaults to a guard of 250,000 distinct sequences **per sample**. Do not disable that guard simply to suppress a failure, or split one library into artificial biological samples.

When resources have been assessed, the complete R1 command is:

```bash
/usr/bin/time -v python "$CSRNA_CODE/csrna_pipeline.py" \
  --manifest "$CSRNA_DATA/samples_recent_batch.csv" \
  --pair-handling r1-only \
  --analysis-mode exploratory \
  --input-scope full_library --pilot-pairs 0 \
  --reference "$CSRNA_REFERENCE" \
  --report-detail summary --quality-detail counts-only \
  --top-fragments 1 --max-fold-variants 1 \
  --threads 4 \
  --output "$CSRNA_DATA/results/recent_batch_full.json" \
  > "$CSRNA_DATA/results/recent_batch_full.stdout.log" \
  2> "$CSRNA_DATA/results/recent_batch_full.log"
echo "Exit code: $?"
```

This retains the default memory guard and intentionally small folding budget. A guard failure requires resource/design review, not interpretation of partial counts. Main-pipeline runs restart analysis; they do not resume from RBC checkpoints. Keep pilot and full reports separately, and compare cell types using one jointly processed main report with consistent scope and settings.

### 4. Optionally run RBC diagnostics after main analysis

Skip this section if only the main landscape analysis is needed. Otherwise, keep the original paired FASTQs accessible and verify the companion:

```bash
CSRNA_TEST_TOOLS=1 python "$CSRNA_CODE/test_rbc_diagnostics.py"
```

Run the RBC pilot using the successful main report, original manifest, and identical reference:

```bash
/usr/bin/time -v python "$CSRNA_CODE/csrna_rbc_diagnostics.py" \
  --analysis-report "$CSRNA_DATA/results/recent_batch_pilot.json" \
  --manifest "$CSRNA_DATA/samples_recent_batch.csv" \
  --reference "$CSRNA_REFERENCE" \
  --output-dir "$CSRNA_DATA/results/rbc_diagnostic_pilot" \
  --index-cache "$CSRNA_CODE/reference_indexes" \
  --pilot-pairs 10000 --sampling-seed 20261007 --threads 4 \
  > "$CSRNA_DATA/results/rbc_diagnostic_pilot.stdout.log" \
  2> "$CSRNA_DATA/results/rbc_diagnostic_pilot.log"
echo "Exit code: $?"
```

Use `recent_batch_full.json` instead if that is the report you wish to associate with the diagnostics. Sample identities and reference checksums must agree. The extension validates all RBC pairs and samples both mates together, then aligns R1 and R2 **independently** in end-to-end and local modes. It does not reconstruct inserts or add local matches to the main feature counts.

After inspecting this pilot and available disk/RAM/time, run all RBC pairs:

```bash
/usr/bin/time -v python "$CSRNA_CODE/csrna_rbc_diagnostics.py" \
  --analysis-report "$CSRNA_DATA/results/recent_batch_pilot.json" \
  --manifest "$CSRNA_DATA/samples_recent_batch.csv" \
  --reference "$CSRNA_REFERENCE" \
  --output-dir "$CSRNA_DATA/results/rbc_diagnostic_full" \
  --index-cache "$CSRNA_CODE/reference_indexes" \
  --pilot-pairs 0 --sampling-seed 20261007 --threads 4 \
  > "$CSRNA_DATA/results/rbc_diagnostic_full.stdout.log" \
  2> "$CSRNA_DATA/results/rbc_diagnostic_full.log"
echo "Exit code: $?"
```

A full RBC run can follow a main **pilot** report: the extension reads the original FASTQs, not just the main pilot's selected sequences. Disk-backed SQLite reduces whole-library Python memory pressure, but storage and all-best-hit alignment volume can be substantial. Full-server throughput remains unverified.

Each RBC workspace contains `rbc_diagnostics.json` for the dashboard, `sequence_evidence.tsv.gz` for detailed sequence evidence, `parent_counts.tsv` for raw-parent counts/RPM, and `diagnostics.sqlite` plus logs/configuration/query data for audit and resume. Pilots also retain sampled paired FASTQs. Keep this workspace to resume. Repeat the exact command to reuse completed stages; changing inputs, code, or options requires a new directory. Never run two writers against one workspace.

Optional origin screens need appropriate supplied FASTAs. Add one or more of these flags to either RBC command, with real paths and a new output directory:

```bash
--screen-reference genome=/absolute/path/to/mouse_genome.fa \
--screen-reference microbial=/absolute/path/to/microbial_reference.fa \
--screen-reference technical=/absolute/path/to/technical_reference.fa
```

Only add `--adapter-sequence` when the actual adapter sequence is confirmed. Omitted screens remain `not_tested`, not negative. Default acceptance gates are `--min-aligned-length 30 --min-aligned-fraction 0.5 --max-edit-rate 0.1`; these are configurable diagnostic thresholds, not universal biological definitions. See [RBC_DIAGNOSTICS.md](RBC_DIAGNOSTICS.md) for all outcome labels, thresholds, provenance checks, and recovery details.

### 5. Open the reports in the web application

**Recommended: analyze on the server and visualize on your workstation.** Only JSON reports are needed by the dashboard; FASTQs, indexes, and SQLite do not need to be copied. Run the following on your **Mac**, not inside the server SSH session (adjust the hostname if your institution provides a different SSH address):

```bash
mkdir -p "$HOME/Documents/wang-lab/reconstructed_pipeline/results/server_batch"
scp sohamms3@huawang3-serv-01:~/sequencing_data/results/recent_batch_pilot.json \
  "$HOME/Documents/wang-lab/reconstructed_pipeline/results/server_batch/"
# Optional: only after RBC diagnostics completes.
scp sohamms3@huawang3-serv-01:~/sequencing_data/results/rbc_diagnostic_pilot/rbc_diagnostics.json \
  "$HOME/Documents/wang-lab/reconstructed_pipeline/results/server_batch/"
```

For full runs substitute the full report/workspace paths; keep pilot and full copies in separate local directories so the RBC filenames do not overwrite one another. The paths inside a copied report may still refer to the server; copying the JSON does not copy the detailed evidence tables.

Install the dashboard dependencies once in the local Python environment you will use, then launch:

```bash
cd "$HOME/Documents/wang-lab/reconstructed_pipeline"
python -m pip install -r csrna_visualizer/requirements.txt
python csrna_visualizer/run.py \
  --report "$PWD/results/server_batch/recent_batch_pilot.json"
```

The existing Mac installation can use `/opt/anaconda3/bin/python` instead of `python` for both commands. The dashboard requires Python >=3.10 and its requirements; Bowtie2 and RNAfold are not needed just to visualize reports.

1. Keep that terminal running and open **http://127.0.0.1:8501/** in your browser.
2. The main report is loaded by `--report`. To change it, expand **Open a pipeline report** in the sidebar and enter its local path or upload its JSON. Clear an uploaded file before switching back to a path; uploads take precedence. Prefer paths for reports above the default 200 MB upload limit.
3. To view the optional report, select **RBC diagnostics** under **Explore**. Enter the local path to `rbc_diagnostics.json` in **RBC diagnostic JSON path**, press Enter, or upload that file. This page uses its own loader; do not pass the RBC JSON to the launcher's main `--report` argument.
4. Return to another sidebar page to view the main report. With only an RBC report, launch `python csrna_visualizer/run.py` without `--report`, then select **RBC diagnostics**.
5. Stop with **Ctrl+C** in the launching terminal. To restart after a code update, run the same command again and refresh the browser. If port 8501 is already occupied, stop the known old launcher or use `--port 8502` and open http://127.0.0.1:8502/.

**Alternative: run the dashboard on the server through SSH.** Install the visualizer requirements in the server environment, launch `python "$CSRNA_CODE/csrna_visualizer/run.py" --report "$CSRNA_DATA/results/recent_batch_pilot.json"` there, and leave it running. In a separate **local** terminal create a tunnel:

```bash
ssh -N -L 8502:127.0.0.1:8501 sohamms3@huawang3-serv-01
```

Open http://127.0.0.1:8502/ locally. Use **server paths** in this dashboard, because its Python process runs on the server. Keep the SSH tunnel open; Ctrl+C closes it. The launcher binds to localhost, so no public Streamlit port is required. If your compute job runs on a different host, use your institution's approved tunnel route to that host.

### 6. Explore, interpret, and export results

For the main report, start with **inclusive** evidence (exact plus accepted-inexact), then compare exact-only support to see whether a pattern depends on mismatched reads.

| Dashboard page | What to inspect | What the evidence supports |
|---|---|---|
| **Overview** | Scope, unmapped fraction, exact/inexact retention, ambiguity, screens actually performed | Whether enough interpretable evidence exists to explore; a high unmapped fraction alone does not identify contamination |
| **RNA composition** | All-read versus mapped-only denominator, RNA categories, abundant parents, individual replicates | Descriptive makeup of the sequenced preparation; mapped-only percentages omit unresolved material |
| **Replicate QC** | Review flags, per-condition composition, metadata/PCR differences, feature- versus parent-level profiles, correlation/PCA, recurrence | Reproducibility and possible outliers; flags do not exclude samples, and PCA is not a significance test |
| **Feature comparisons** | Shared counts/RPM, selected samples and track, scope restrictions | Descriptive differences under consistent processing; subset fold rankings remain disabled |
| **Parent explorer** | Supported endpoint pairs, coverage, full-length-compatible/discrete/diffuse patterns | Locations of observed read evidence on a parent; peaks do not prove native molecular ends |
| **Feature library** | Support tiers, number of libraries detected, sequences, available reference and variant folds | Candidate features for follow-up; singleton evidence and predicted structures require particular caution |
| **RBC diagnostics** (optional) | Library/mate selection, mode outcomes, RNA composition, top parents, complexity, clipping, origin screens, pair agreement | Whether local matching recovers additional sequence evidence and which unresolved sequences merit investigation |

In RBC diagnostics, inspect each replicate and mate separately. `local_only` can reflect partial RNA matches, adapter sequence, or other clipped material; it is not automatically a valid fragment species. `neither_accepted` includes sequences with no match and sequences whose reported matches failed the gates. Short matches below 30 nt fail the default span gate, so absence from accepted counts does not establish biological absence. Homopolymer/low-complexity labels overlap and do not remove reads. Compare read-weighted versus distinct-sequence-weighted frequencies to distinguish repeated sequences from sequence diversity.

Genome/microbial/technical matches are non-exclusive evidence, not definitive source calls. Independent-mate agreement is not a native paired alignment or insert-length test. Do not sum R1 and R2 counts as independent molecules, merge RBC local-match counts into main endpoint matrices, or compare full-depth RBC diagnostics directly with other cells' small pilots. Use the same joint main report for cross-cell comparisons.

Neither report proves cell-surface localization, glycosylation, modification, degradation, or biological function. Replication supports consistency, but unknown preparation, PCR amplification, reference coverage, and lack of negative controls still constrain interpretation. Read boundaries remain observed endpoints until insert completeness is established.

To share results, use **Download table · CSV** for tables with scope/provenance labels, chart camera icons for SVG, and **Prepare interactive HTML export** for offline interactive charts. The data grid's own CSV shortcut may omit the added provenance columns. Export raw integer feature counts from **Replicate QC** for separately designed statistical analysis. Record report scope, evidence track, normalization, reference, and diagnostic thresholds alongside every figure; retain the original JSON and logs.

### 7. Troubleshoot and rerun

| Symptom | Action |
|---|---|
| `conda: command not found` or missing executable | Initialize the correct Conda installation, activate the environment, and rerun `command -v` and tool tests |
| Missing FASTQ, unequal mate counts, or mismatched names | Verify the manifest against actual files; preserve matching R1/R2 pairs and correct the input before retrying |
| Main output already exists | Choose a new filename; use `--overwrite` only when intentionally replacing that report |
| Distinct-sequence guard or memory failure | Assess RAM and increasing pilot sizes; a failed run is not a partial successful library analysis |
| RBC workspace configuration/checksum mismatch | Use a new output directory for changed settings, source code, or inputs; keep the shared index cache |
| RBC `.lock` remains after an interruption | Check that no process is using that workspace/cache before removing a stale lock; never remove an active run's lock |
| Report rejected by dashboard | Use main schema-2 JSON in the main loader and completed RBC diagnostic JSON on **RBC diagnostics**; do not load logs or TSV files as JSON |
| Empty/undefined correlation, PCA, or missing folds | Review eligible counts and filters; constant profiles have undefined correlations, identical profiles have no ordination, and folding is budget-limited |
| Dashboard is stale or port is occupied | Restart the correct launcher, refresh the browser, and confirm the report path; use an alternate port if necessary |

For tests after later code updates, run the two main checks above, `CSRNA_TEST_TOOLS=1 python "$CSRNA_CODE/test_rbc_diagnostics.py"`, and `python csrna_visualizer/test_visualizer.py` from the local `reconstructed_pipeline` directory. A passing test suite is evidence for the tested cases, not a guarantee about every input or server capacity.

## Technical reference and historical validation

The sections below preserve design details, annotations, and dated validation records. Use the walkthrough above for the current paired batch and optional RBC workflow.

## Complete-library validation and scale limits (2026-09-26)

The three supplied trimmed R1 FASTQ files total approximately 110 GiB. Complete FASTQ format/base/quality validation and representative end-to-end mapping are distinct checks: validating every record does not mean every read has been aligned. This implementation retains unique sequences, per-sequence positional quality summaries, mappings, and the complete JSON report in memory. It is not a bounded-memory production engine for hundreds of millions of reads.

A new `--max-distinct-sequences` limit (default 250,000 per library) aborts ingestion explicitly when exceeded. It never truncates data or publishes partial results. The limit also applies to FASTA and CSV; duplicate counts continue accumulating normally. Set `0` to disable or choose a higher limit only after assessing resources. This is an ingestion guard, not a RAM guarantee: reference indexing, all-best-hit alignment multiplicity, and report construction also consume memory. Complete libraries exceeding the limit require an adequately provisioned execution environment or a disk-backed redesign; dividing them into separately analyzed samples would change clustering and normalization and is not an equivalent workaround.

FASTQ quality aggregation now allocates quality arrays only for a new sequence, avoiding discarded allocations for every duplicate read. Counts and quality statistics are unchanged.

A parent-landscape reporting error was fixed: full-length-compatible support now uses each sample’s actual read coordinates. Previously, a pooled feature representative could make one sample’s near-boundary reads appear full-length because another sample supplied the representative. Shared feature identity, counts, exact/inexact tracks, and fold comparisons remain unchanged. A regression test covers two samples sharing endpoint clusters but differing in actual full-parent support.

Proceeding to biological fragment-species interpretation additionally requires confirmation that read endpoints correspond to insert endpoints (original read cycles, adapter-trimming evidence, paired-end availability, and library/UMI design). Enrichment rankings remain descriptive and extraction protocols differ. No genome or microbial reference or experiment-specific adapter sequence has been supplied for the unresolved-read investigation; those screens remain untested. Current software tests cannot establish an unmapped read’s identity or prove a discrete biological RNA molecule.

### Complete FASTQ audit results

All **345,055,492 records** were inspected in a streaming four-line FASTQ audit: complete records, matching sequence/quality lengths, valid A/C/G/T/N bases, printable Phred+33 quality characters, and valid header/plus lines. No format failures were found. The following length and quality totals cover every supplied record, not a subset.

| Library | Reads | Distinct sequences in first 1,000,000 reads | Reads 148–150 nt | Bases below Q20 |
|---|---:|---:|---:|---:|
| M1 | 97,962,614 | 616,073 | 59.21% | 0.58% |
| BMDM | 139,086,742 | 210,026 | 51.25% | 0.72% |
| BMDC | 108,006,136 | 180,079 | 41.29% | 0.58% |

Maximum observed read length is 150 nt in every library. Concentration near that maximum is compatible with sequencing-length censoring; it does not prove censoring or identify the trimming protocol. Without insert/adapter provenance, these boundaries remain sequenced-read endpoints rather than validated molecular ends. Distinct counts in the first million reads are descriptive and must not be extrapolated as whole-library totals.

The M1 full-file ingestion guard was exercised directly and stopped explicitly at 250,000 distinct sequences without producing a partial report. The older 229-sequence/83-feature fixture retained identical features, counts, folds, support reports, and parent-landscape results after the fixes. Regression fixtures specifically exercise the corrected cross-sample endpoint case.

### Full-reference tests and readiness decision

A pre-fix baseline run on the first 200 reads of each library completed against all 576,495 combined-reference records. After fixes, a deterministic uniform reservoir of 2,000 reads from across **each entire file** (6,000 reads total, seed 20260926) completed against the same full reference. No rRNA-only reference filter was used. These runs explicitly used `--input-scope selected_subset`; candidates and biological fold rankings were disabled. Folding budgets were `--top-fragments 1 --max-fold-variants 1` to keep the validation run focused; mapping and count eligibility were not restricted by those budgets.

| Across-file validation subset | Primary-unmapped reads | Fraction of 2,000 sampled reads |
|---|---:|---:|
| M1 | 1,071 | 53.55% |
| BMDM | 180 | 9.00% |
| BMDC | 553 | 27.65% |

These are subset estimates, not exact full-library mapping totals. The reservoir contained 4,314 distinct sequences and produced 1,133 fragment features. Independent checks re-read the FASTQs, reproduced per-sequence counts and positional quality summaries, verified feature membership and paired endpoints, reproduced all count matrices, and walked all 6,308 reported raw alignments against the reference to verify CIGAR spans, edit distances, and exactness. Reference sequences and successful fold lengths/energies were also checked. All checks passed. The three maintained deliverables are synchronized; 46 regression tests, four synthetic Bowtie2/RNAfold workflows, and fresh notebook execution passed.

**Readiness: do not yet treat Goal 1/Goal 2 biological conclusions as validated.** Whole-file format validation is complete, but full-library alignment/feature analysis has not been completed. The single JSON report for only this 6,000-read validation run is approximately 328 MiB, largely due to whole-reference metadata; full-scale memory and report storage need explicit engineering attention. The ingestion guard is a controlled failure, not a scalability solution.

The next implementation should use disk-backed sequence/quality aggregation, reusable verified reference indexes, bounded alignment batches preserving every tied best-score hit, one global endpoint-clustering pass, and streamed final reporting. Batches must not independently define features or sample denominators. Validate equivalence against the existing fixtures before a complete run. Alternatively, first profile progressively larger subsets to provision a sufficiently large-memory execution environment; simply disabling the guard is not evidence of capacity. Temporary on-disk state will be needed for a bounded-memory implementation and can be cleaned after successful completion.

Before interpreting fragment species, obtain original read-cycle and adapter/insert-completeness information and any paired-end/UMI metadata. Complete the unmapped-read investigation with appropriate genome/microbial references and actual library adapters. The high sampled M1 unmapped fraction remains unresolved; no origin or contamination label has been inferred. The supplied libraries have no biological replicates established in metadata, so later comparisons remain descriptive.

## What changed and what stayed reliable

| Design issue | Version 2 behavior |
|---|---|
| Detailed analysis restricted to rRNA | Paired-endpoint features and coverage summaries for every annotated RNA category, including explicitly unclassified references |
| Exact matching determines all target support | Disjoint exact and accepted-inexact tracks; inclusive support is their sum |
| Different eligibility rules in target and controls | Every feature is counted with identical assignment, endpoint, and edit-rate rules in all samples |
| Duplicate reference records treated as unusable biology | Identical full sequences grouped; exact precursor/mature projections validated; raw mapping ambiguity remains visible |
| RNA family known but coordinates unresolved | Family-level support retained once per read; no fabricated locus or fragment endpoint |
| Overlapping reads merged into assumed long molecules | Shared, paired endpoint features; merged coverage reported separately and not folded by default |
| Per-sequence RPM threshold loses sequence variants | All eligible reads enter feature discovery; abundance screening occurs after feature counts are aggregated |
| Unmapped data reduced to summary numbers | Every mapped-length sequence and its counts retained, with headline unmapped QC and optional secondary diagnostics |
| Folding a reference span conflates sequence/structure | Observed fragment variants and reference sequences folded separately |
| Single-sample outputs difficult to compare | Shared exact/inexact/inclusive matrices plus sample metadata embedded in the final JSON |
| Risk of comparing truncated Top100 datasets | Explicit selected-subset mode allows validation but disables biological fold rankings/candidates |

Checked tool execution, validated inputs, deterministic IDs, positive-overlap coverage merging, reference-bound checks, zero-based half-open coordinates, SHA-256 provenance, RNAfold validation, and atomic output publication are retained. The previous 17 regression checks remain, with additional tests for the redesign.

## Biological interpretation

The pipeline measures sequence evidence and descriptive abundance. It does not establish cell-surface localization, modification chemistry, RNA glycosylation, diagnostic exclusivity, probe accessibility, or the existence of one intact merged molecule. A substitution is retained as evidence, not labeled a modified base. Ordinary RNAfold models canonical-base sequences; it does not model glycan attachment or unknown modified nucleotides.

The supervisor reported a large unmapped M1 fraction and unequal mismatch burdens between M1 and BMDM. Those full-library observations motivated the redesign, but cannot be independently reproduced using the supplied Top100 tables. Different extraction protocols, PCR duplication, sequencing depth, reference coverage, and RNA composition still affect interpretation. Repeated sequence counts are not independent molecules without appropriate UMI/duplicate information.

The modified-base rationale is biologically plausible: glycoRNA studies established glycosylated cell-surface small RNAs and identified acp³U as an N-glycan attachment site; modification-sensitive sequencing studies show that modifications can affect reverse transcription. None of these findings establishes the cause of the mismatches in these particular samples. Sources: [Flynn et al.](https://pmc.ncbi.nlm.nih.gov/articles/9097497/), [acp³U study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11571744/), [ARM-seq](https://www.nature.com/articles/nmeth.3508).

## Requirements

- Python 3.10 or newer; the production script uses only the standard library.
- `bowtie2`, `bowtie2-build`, and `RNAfold` on the same terminal/kernel PATH.
- JupyterLab and an IPython kernel only if using the notebook.
- Complete single-end libraries, a suitable RNA reference, and optional sample metadata/diagnostic references.

The implementation was tested with Bowtie2 2.5.5 and RNAfold 2.7.2. Tool versions and executable paths are saved in each report.

```bash
conda create -n csrna --override-channels -c conda-forge -c bioconda python=3.12 bowtie2=2.5.5 viennarna=2.7.2 -y
conda activate csrna
python -m pip install jupyterlab ipykernel
```

The final command is unnecessary for standalone-script use. Existing working environments can be used without reinstalling packages.

## Full-library inputs and quality evidence

Supported extensions are `.fastq`/`.fq`, `.fasta`/`.fa`/`.fna`, and `.csv`, with optional `.gz` compression.

- FASTQ uses strict four-line records and Phred+33 quality encoding. Record structure, sequence/quality length, quality characters, and sequence alphabet are validated. Uppercase and U-to-T normalization precede counting.
- For each distinct FASTQ sequence in each sample, per-position Phred sums, counts below Q20, and number of observations are retained. Divide a position's sum by its record count for its mean quality. Mismatch events refer to the original read positions, including on reverse mappings.
- Mapping still deduplicates sequences and uses FASTA queries, applying the same quality-independent Bowtie2 scoring to all samples. Quality summaries are diagnostic evidence, not a claim that alignment is quality-weighted. Individual read names and individual quality strings are not retained.
- Multiline FASTA is supported; each record counts once. Abundance encoded in FASTA headers is not guessed. For collapsed FASTA, provide a correctly decoded full-library sequence-count CSV instead.
- CSV requires `sequence` and exactly one integer `count` or `counts` column. Duplicate normalized sequences are summed. Full-library CSV inputs require `--counts-complete`; zero rows are ignored, negative/fractional/missing counts rejected, and an empty library is an error.
- `master_counts.csv` and gene-level featureCounts tables are not sequence-count inputs. Legacy selected rRNA RPM tables are rejected in full-library mode, even with `--counts-complete`.

Adapter removal, quality trimming, UMI deduplication, paired-end merging, and demultiplexing are not silently performed. Use protocol-appropriate single-end inputs. The optional adapter diagnostic reports exact sequence evidence; it does not trim reads. Phred+64 FASTQ must be converted before use.

## Run the updated pipeline

From the repository root, substitute real full-library input paths:

```bash
python reconstructed_pipeline/csrna_pipeline.py \
  --sample M1=/absolute/path/M1.fastq.gz \
  --sample BMDM=/absolute/path/BMDM.fastq.gz \
  --sample BMDC=/absolute/path/BMDC.fastq.gz \
  --reference downstreamAnalysis/mouse_combined_with_rRNA_clean_dedup.fa \
  --output output/csrna_results.json
```

Add `--counts-complete` for complete raw sequence-count CSVs. The default target is M1 and the default controls are BMDM and BMDC. To use other sample names, supply `--target NAME` and repeat `--control NAME`; every declared target/control must exist. Additional samples may participate in the shared feature table. Replicates are stored as separate sample columns and are not pooled automatically.

The script builds the matching reference index in a temporary directory. `--reference-id ID` may be repeated to restrict the reference universe, but such a restriction limits interpretation. An rRNA-only reference cannot characterize other RNA classes or establish genome-wide specificity. Quote IDs containing shell metacharacters such as `|`.

Optional metadata:

```csv
sample_id,condition,replicate,library_preparation
M1,activated,1,csRNA_surface
BMDM,resting,1,Trizol
BMDC,dendritic,1,record_actual_protocol
```

Supply it with `--sample-metadata /path/samples.csv`. IDs must cover exactly the provided sample names, with no duplicate rows. Without it, condition defaults to the sample name, replicate is null, and preparation is unspecified. These defaults are not fabricated biological annotations. This version produces descriptive fold rankings only, not DESeq2 tests, p-values, significance, or inferred replicates.

## How reads become shared fragment features

1. Count all supplied reads before alignment or selection. Keep short reads in input totals and an explicit excluded-sequence table.
2. Map each distinct sequence once across all samples using Bowtie2 `--end-to-end --very-sensitive -a --seed 0`. Retain all reported equal-best AS hits. A primary alignment alone is not uniqueness evidence.
3. Resolve identical sequence records and validated coordinate projections. Keep family-level identity separate from coordinate identity and raw genomic-locus ambiguity.
4. Admit a read to endpoint features only when it resolves to one oriented parent interval, has a supported first and last aligned query base, and meets the edit-rate limit. Exactness is a separate property, not a universal eligibility gate.
5. Pool eligible endpoints across samples. Cluster reference starts and ends separately with a configurable tolerance, then create a feature only for a start/end cluster pair actually linked by a spanning read. No unobserved endpoint combinations are manufactured.
6. Count every sample against this same feature definition. Count each sequence's abundance once per feature; a sequence can belong to at most one feature. Produce exact, inexact, and inclusive matrices with identical row/column order.
7. Apply abundance/enrichment thresholds after aggregation, using symmetric eligibility for target and control feature counts. Run the broader coordinate-overlap screen separately.
8. Characterize and fold selected observed variants and their reference counterpart. Retain the complete evidence and publish one final JSON.

Endpoint clusters use a deterministic bounded-diameter rule: the sorted cluster span is at most `2 × tolerance`, and its midpoint is the cluster center. Every member is within the tolerance of that center. This avoids chaining gradually drifting endpoints into one large cluster. It is a reproducible heuristic, not a fitted model of biological cleavage. Changing tolerance or adding new endpoint observations can change feature boundaries/IDs; compare samples using one joint run, not independently discovered rows.

Endpoints are initially sequencing-read boundaries. If sequencing stops at a fixed cycle limit (for example, 150 nt), the read end need not be the original RNA molecule end. Check insert-length/adapter and library-preparation evidence before interpreting a paired-read feature as a discrete physical RNA species. The report explicitly labels this distinction; no automatic completeness claim is made.

By default the tolerance is 2 nt. A feature's representative start/end pair comes from an actual supporting read, chosen by pooled abundance with a deterministic tie-break. Minus-orientation features reverse the interpretation of reference start/end as molecular 5′/3′ boundaries. All coordinates remain zero-based, half-open on the named reference. Read orientation is not necessarily biological RNA strand without library-preparation information.

## Assignment hierarchy and reference equivalence

| Assignment state | Interpretation | Used in endpoint features? |
|---|---|---|
| `unique_locus` | One reported best raw locus | Yes, if boundary/edit criteria pass |
| `equivalent_parent_resolved` | Multiple raw hits resolve to one sequence-equivalent oriented parent interval | Yes, if boundary/edit criteria pass; genomic copy remains unresolved |
| `family_resolved_coordinates_ambiguous` | RNA family supported, but positions/orientations remain ambiguous | No; counted once in family support |
| `unresolved_parents` | More than one unresolved parent/family remains | No; full evidence retained |
| `unmapped` | No alignment in the primary reference under this search | No; retained for diagnostics |

Identical **full reference sequences within the same category** are grouped after direct sequence equality verification. Genomic loci remain listed in the raw alignments. Automatic rRNA projections are restricted to annotated precursor/repeating-unit sources and annotated mature 18S/28S/5.8S/5S targets whose entire sequence is an exact substring of that source, in the indicated orientation. A read must lie completely inside the verified target interval to be projected. No similar-looking coordinate or partial sequence similarity is treated as an equivalence proof. Palindromic orientation ambiguity is retained.

For the supplied reference, full sequence verification finds:

| Source | Target | Zero-based source offset | Target length |
|---|---|---:|---:|
| `RIBO|BK000964.3` | `RIBO|NR_003278.3` (18S) | 4007 | 1870 |
| `RIBO|BK000964.3` | `RIBO|NR_003279.1` (28S) | 8122 | 4730 |

The other supplied RIBO records do not satisfy this exact full-target containment rule and are not assigned an invented transform. Mature and precursor molecules are not thereby declared biologically identical: the read sequence supports a mature-RNA coordinate region but may not establish which precursor/mature molecule or genomic copy produced it.

Use `--no-auto-rrna-projections` to disable automatic containment projections. Other validated projections can be supplied explicitly, for example:

```bash
--reference-projection 'RIBO|BK000964.3=RIBO|NR_003279.1:8122:+'
```

The syntax is `SOURCE=TARGET:OFFSET:STRAND`, with zero-based offset on the source; the complete oriented target sequence must match the source slice. Invalid, cross-category, cyclic/chained projections fail. `--parent-family 'REFERENCE_ID=FAMILY'` declares a family grouping without asserting coordinate equivalence. Obvious rRNA subtype labels and GENCODE gene identifiers are used for family annotation when available; this is not a taxonomic assignment or a modification call.

## Step 3 — exact-only and inexact sensitivity reporting

The terminal now prints `SUPPORT TRACKS` for each sample. The JSON adds `support_reporting.samples[SAMPLE]` with counts, RPM, inexact fraction of eligible support, and inclusive-to-exact support ratio, both overall and `by_category`. Scope labels preserve the distinction between complete supplied libraries and selected subsets. Zero denominators for support fractions/ratios produce `null`, never infinity or an invented zero. The RPM denominator is exactly the existing sample denominator, shared by all tracks.

The existing `count_matrices.exact`, `.inexact`, and `.inclusive` remain unchanged. Exact-only plus accepted-inexact equals inclusive support. Reporting uses the same shared fragment features; it does not reconstruct different endpoint clusters for each track. Inexact means non-exact alignments accepted under the existing endpoint and edit-rate rules, including allowed substitutions or internal indels. Rejected or unresolved mappings contribute to neither feature-support track and remain in the existing exclusion/assignment reports. Counts describe reads, not independent molecules or demonstrated RNA modifications.

To make “unique-but-inexact” explicit without removing correct reference-equivalence handling, `support_reporting.assignment_count_matrices` contains four disjoint matrices:

| Matrix | Evidence |
|---|---|
| `raw_unique_exact` | One raw best-score alignment; exact and feature eligible |
| `raw_unique_inexact` | One raw best-score alignment; non-exact and feature eligible |
| `equivalence_resolved_exact` | Multiple raw alignments resolved to one oriented parent interval; exact and feature eligible |
| `equivalence_resolved_inexact` | Multiple raw alignments resolved to one oriented parent interval; non-exact and feature eligible |

These matrices share feature/sample order with the existing count matrices. Summing the two exact components reproduces the exact matrix; summing the two inexact components reproduces the inexact matrix. The four components together reproduce inclusive counts. The strict raw-unique sensitivity track is available separately; accepted-inexact totals also retain equivalence-resolved evidence to avoid silently discarding redundant-reference rRNA support. Family-only or unresolved-coordinate assignments are never promoted to unique support.

All sample and category totals are read-weighted, not distinct-sequence counts. Every feature/sample component sum is checked against existing counts before publication. A mapped category with no feature-eligible support has explicit zero totals; categories represented in the shared features may also have zero support in an individual sample. Ratios are descriptive measures of added read support, not treatment fold changes. Existing exact/inclusive candidate comparisons, normalization, folding, mapping, and all prior output fields retain their behavior. No additional permanent output files are created.

Three new regression tests verify disjoint component accounting, RPM and ratio definitions, zero support, selected-subset labeling, rejection of invalid members/counts, and non-mutation of existing features. Real Bowtie2/RNAfold fixtures verify exact/inexact totals and reference-equivalence components.

## Evidence tracks, comparison, and structural interpretation

The feature-count tracks are:

- `exact`: a resolved interval with all supporting best raw alignments verified as zero-edit, ungapped, exact reference sequence matches.
- `inexact`: the same assignment and endpoint requirements, but one or more accepted differences from the reference; default `NM / query_length <= 0.10`.
- `inclusive`: exact + inexact, without double counting.

Internal insertions/deletions can contribute endpoint evidence when both query ends are aligned and the edit limit passes. Terminal indels, clipping, and reference skips do not establish the required endpoint pair. Substitution and indel events remain in each raw alignment, with CIGAR-aware positions. Observed query sequences are never reconstructed by concatenating mismatched suffixes or silently corrected to the reference.

The default RPM denominator is every supplied positive-count input read, including unmapped and short reads. `--denominator mapped` explicitly uses reads with at least one primary-reference alignment, counted once even if multi-mapped. A zero mapped denominator fails. Category RPM always uses the input total, independently of the feature denominator setting. Neither denominator corrects extraction bias or explains unmapped biology.

For each feature and for both exact and inclusive tracks, the descriptive ratio is:

```text
(target feature RPM + pseudocount RPM) / (control feature RPM + pseudocount RPM)
```

This is calculated separately against every control. The default candidate flag uses inclusive feature support, minimum RPM 100, and minimum smoothed enrichment 4 against each control. `--require-control-absence` additionally requires no assigned control reads for that same feature. `--max-homopolymer` is an optional screen on the representative observed sequence; all sequence variants remain available and must be evaluated separately for actual probe design.

`control_overlap_screen` is independent: any best-mapped control read overlapping the feature's observed endpoint envelope contributes once, including inexact and ambiguously placed reads. Default screening includes both orientations; `--control-strand same` requires compatible stranded libraries. A feature may be descriptively enriched even while this screen detects substantial control overlap. Neither candidate status nor feature absence asserts probe specificity. Near-match off-target binding outside the supplied reference/best-hit evidence is not exhaustively searched.

`family_count_matrix` preserves reads that resolve to a family even if endpoint assignment or edit-rate eligibility fails. It answers a broader assignment question than the feature matrices; the two must not be added together.

Folding is performed separately for observed variants and representative reference slices, within a per-category budget. A reference slice can have a different length from an indel-containing read. Canonical-base RNAfold does not account for unspecified nucleotide modifications. MFE and unpaired fraction in one predicted structure are descriptive, not measured probe accessibility. Regions formed only by overlapping coverage are not folded unless `--fold-coverage-regions` is explicitly set, and remain labeled coverage summaries.

## Step 2 — headline unmapped QC and source investigation

The terminal prints `HEADLINE QC` per sample: primary-unmapped read count divided by **all supplied input reads**, irrespective of the RPM normalization option. Below-minimum-length reads are shown separately: they were not aligned and are not called unmapped. Selected inputs are labeled `SELECTED SUBSET ONLY`; their fraction cannot estimate the whole-library unmapped fraction.

The JSON `unmapped_qc[SAMPLE]` provides headline counts/fractions, a disjoint `partition_counts` table, and four diagnostic screens: `genome`, `microbial`, `adapter`, and `other`. Secondary references require explicit categories; names alone never establish taxonomy. Existing uncategorized diagnostic references remain supported as `other`.

```bash
# Append these options to your normal full-input pipeline command:
--diagnostic-reference host=/path/to/mouse_genome.fa \
--diagnostic-reference-category host=genome \
--diagnostic-reference microbes=/path/to/lab_approved_microbes.fa \
--diagnostic-reference-category microbes=microbial \
--adapter-sequence ACGTACGTACGT
```

The adapter above is an illustrative placeholder, not a recommended experimental adapter. Use the actual library adapter sequence (at least 8 bases). Adapter screening detects exact internal or terminal substrings in either orientation; it does not trim, detect all partial adapters, or exclude reads. Secondary alignment uses the existing unspliced end-to-end Bowtie2 search; a no-match result is conditional on supplied references and search settings. Large genome references can require substantial indexing memory and temporary disk space. References are local, hashed, and temporary alignment files are cleaned as before.

Each primary-unmapped fragment adds `unmapped_diagnostic_labels`: zero or more evidence labels; mapped fragments use `null`. Sample `screens` gives read-weighted matching counts and fractions of primary-unmapped reads, with `not_tested` (reference/adapter absent), `tested`, or `no_unmapped_reads` (configured but no queries). Unprovided screens have null counts and fractions, not zero negatives. Fractions with a zero unmapped denominator are null.

`partition_counts` counts each unmapped read exactly once: `genome_evidence`, `microbial_evidence`, `adapter_evidence`, `other_evidence`, `multiple_source_evidence`, `no_match_in_supplied_screens`, or `not_tested`. Screen counts can overlap and must not be summed. Reads matching both genome and microbial sources retain both labels; they are not assigned a definitive origin or called contamination. Multiple references in one category count once per read. A conservation checkpoint requires partition totals to equal the original unmapped count. Source categories are recorded in `diagnostic_reference_categories`; raw secondary alignments and source hashes remain available in existing diagnostic fields.

This adds reporting only. Primary assignments, feature eligibility, matrices, folds, and candidates are unchanged. Determining the identity of M1’s unmapped population requires full library input and appropriate lab-approved references/adapters. Without them, screens remain explicitly untested. Synthetic real-tool fixtures verify genome/adapter overlap; regression checks cover microbial overlap, missing screens, empty denominators, invalid labels, and read accounting.

## Unmapped and ambiguous evidence

The terminal and JSON report give read-weighted unmapped fractions per input library. Per-category retention tables show raw ambiguous reads, retained exact/inexact feature reads, and excluded reads. **Step 1 of the supervisor’s next steps: prominent category ambiguity QC.** Each sample/category now emits a `HIGH MAPPING AMBIGUITY` warning to stderr when the read-weighted raw ambiguity fraction reaches `--ambiguity-warning-fraction` (default 0.5, inclusive boundary). A conspicuous summary banner follows the warnings. The numerator counts reads with multiple raw best-score loci; the denominator is mapped reads assigned to that category, before equivalence handling. These are read counts, not counts of distinct sequences. Cross-category mappings stay in `multi_category_ambiguous`, counted once, without redistribution among categories.

The JSON adds `samples[SAMPLE].category_ambiguity_qc` for every mapped category and top-level `qc_alerts` for triggered warnings. Each alert has stable code `HIGH_CATEGORY_MAPPING_AMBIGUITY`, severity, sample, category, input scope, threshold, ambiguous count/fraction, retained feature count/fraction, and message. Messages also remain in the existing `warnings` list. For example, inspect `[a["message"] for a in report["qc_alerts"]` after loading the JSON.

A high raw ambiguity rate may coexist with high resolved-parent retention; both are reported. Retained support includes all feature-eligible reads in the category, not just the ambiguous subset. Exclusion can reflect ambiguity or other eligibility filters, so the difference is not labeled unresolved ambiguity. Review `category_retention`, `assignment_state_counts`, and per-read evidence before interpretation. Missing/zero mapped support has no defined fraction (`null`) and does not warn; a zero threshold still requires at least one ambiguous read. Unmapped and below-length reads are excluded from this denominator and reported separately. `SELECTED SUBSET ONLY` explicitly labels warnings from selected inputs. The feature is diagnostic only: it does not change alignments, eligibility, counts, normalization, candidates, folding, or exit status, and creates no additional output files.

All sequences of mapping length are stored once with counts in every sample. `unmapped_read_ids` identifies those without a primary-reference hit; `unresolved_read_ids` also includes mapped sequences excluded from endpoint features. Short excluded sequences have their own count table. Absence from a Top100 display is never a reason to discard evidence.

Optional local diagnostics for primary-unmapped reads:

```bash
--diagnostic-reference genome=/path/mouse_genome.fa \
--diagnostic-reference microbial=/path/curated_microbes.fa \
--adapter-sequence ACGTACGTACGT
```

These are user-supplied references/sequences; the script does not download databases, transmit reads, or invent identities. Secondary alignment uses the same Bowtie2 search and records hits from all supplied diagnostic sources independently. Diagnostic matches do not silently change primary assignments or denominators. Adapter checks are exact substring checks including reverse complements, with a minimum 8-base supplied sequence; they are neither trimming nor a comprehensive adapter detector. Genome alignment here is unspliced, so a failed diagnostic match cannot rule out genomic origin. A microbial-reference match alone does not prove contamination or organism identity.

## Parameters and resource controls

| Option | Default | Purpose |
|---|---|---|
| `--input-scope` | `full_library` | `selected_subset` permits selected-table validation but disables comparisons/rankings/candidates |
| `--denominator` | `total_input` | Feature RPM population; `mapped` is explicit and conditional |
| `--min-read-length` | 15 | Reads below this length are counted but not aligned |
| `--endpoint-tolerance` | 2 | Endpoint clustering tolerance in nt, applied separately to each reference/strand |
| `--max-edit-rate` | 0.10 | Maximum NM/query length for endpoint feature eligibility |
| `--min-rpm` | 100 | Feature-level target abundance screen, after aggregation |
| `--min-enrichment` | 4 | Minimum smoothed feature enrichment against each control |
| `--pseudocount-rpm` | 1 | Smoothing in RPM units |
| `--control-strand` | `both` | Policy for the separate control-overlap screen |
| `--ambiguity-warning-fraction` | 0.5 | Warn at or above this fraction of mapped category reads with multiple raw best-score loci; diagnostic only |
| `--temperature` | 37 | RNAfold temperature in Celsius |
| `--max-fold-length` | 2000 | Explicit skipped status and null energy beyond this length |
| `--fold-timeout` | 120 | Seconds per RNAfold call; failures/timeouts stop publication |
| `--top-fragments` | 100 | Per-category/per-sample feature and read folding/ranking budget |
| `--max-fold-variants` | 5 | Maximum abundance-ranked observed sequence variants folded per selected feature |
| `--min-landscape-reads` | 10 | Minimum eligible parent support for pattern description |
| `--min-species-reads` | 3 | Read-count support for a recurrent feature; not an independent-molecule threshold |
| `--landscape-dominance` | 0.5 | Fraction needed for full-length-compatible or recurrent-fragment descriptions |
| `--threads` | 4 | Bowtie2/index threads |
| `--temp-dir` | system temp | Optional scratch parent; run directory is unique and cleaned |
| `--overwrite` | off | Permit atomic replacement of a completed output |

No fixed homopolymer or HCR start/end constraint is invented. The optional sequence screen and all thresholds are explicit in the report. Vary endpoint/edit thresholds deliberately to assess sensitivity; record the settings rather than interpreting every threshold-dependent feature as a biological species.

Parent landscape labels are descriptive and per sample/reference/orientation. Low parent support gives `insufficient_support`. Predominantly full-length-compatible features give `full_length_compatible`; otherwise a sufficient fraction of read support in recurrent paired-endpoint features gives `discrete_fragment_pattern`. Remaining support is labeled `diffuse_or_multiple_fragments`: it is compatible with a degradation smear but cannot distinguish degradation from multiple real fragment populations or technical artifacts. Raw counts without UMI information do not establish independent recurrence.

## Output and checkpoint guarantees

A successful run publishes one strict JSON report. It embeds the feature matrices and sample-metadata table rather than writing extra intermediate CSVs. Each matrix contains `feature_ids`, `sample_ids`, and a two-dimensional integer `data` array with matching row/column order. These raw counts can later be exported for an appropriate replicate-aware model; no such statistical test is performed now.

| Report field | Meaning |
|---|---|
| `schema_version`, `created_utc`, `status` | Format version 2.0.0 and run completion |
| `parameters`, `tools`, `input_provenance`, `reference` | Configuration, versions, every input hash, annotations, verified projections and aliases |
| `samples`, `sample_metadata`, `warnings` | Read-weighted QC, assignment/retention categories, denominator, condition/protocol metadata |
| `features`, `count_matrices` | Shared paired-endpoint features and disjoint exact/inexact plus inclusive matrices |
| `endpoint_profiles`, `parent_landscape` | Per-sample endpoint histograms and descriptive parent patterns |
| `family_count_matrix` | Family-level counts including coordinate-unresolved support, counted once per read |
| `fragments` | Every mapping-length sequence, raw counts/RPM, qualities, best raw hits, projected assignment, mismatch/indel events and folding |
| `unmapped_read_ids`, `unresolved_read_ids`, `excluded_short_sequences` | Explicit access to unassigned/excluded evidence |
| `diagnostic_reference_provenance` | Diagnostic sources tested and their hashes |
| `top_feature_ids`, `top_fragment_ids` | Per-category display/folding selections, never count filters |
| `candidate_feature_ids`, `fold_change_ranking` | Inclusive descriptive feature screening and ranking; empty for selected subsets |
| `regions` | Secondary merged coverage summaries, explicitly not observed molecules |

The legacy `parent_rRNA` field remains an alias on raw fragment/coverage records for compatibility; use `parent_RNA` and category in new code. Fragment raw coordinates refer to raw alignments; canonical/projected coordinates are in `assignment.locations` and the shared feature records. `candidate_region_ids` remains an empty compatibility field: version 2 candidates refer to fragment features, not coverage regions. Consumers must check `schema_version` and migrate candidate interpretation accordingly. A unique raw locus and a resolved RNA-level interval are different fields by design.

Runtime checkpoints require input totals to equal category/state accounting, each sequence to belong to no more than one endpoint feature, and inclusive matrix totals to equal eligible read counts. Reference slices and CIGAR spans are checked, exactness is verified against sequence, and RNAfold output must have valid sequence/structure/energy. Hashes are checked again before publication. Invalid inputs, failed tools, malformed folds, and accounting errors stop publication rather than yielding fabricated zero energies or partial results.

Only the final report persists. SAM is streamed; temporary reference/query FASTA, index files, and diagnostic logs are removed on normal/error exits. RNAfold uses `--noPS`; no BAM, SAM, PNG, spreadsheet, or PostScript result is created. Atomic publication briefly needs one temporary JSON sibling. Without `--overwrite`, existing outputs cannot be clobbered; with it, an old report survives a failed replacement. Empty successful runs write explicit empty matrices/features. A forced OS kill or power failure may leave a scratch directory; an old report is not a successful new run.

Reference sequences, distinct reads, quality summaries, best-hit evidence, and the final report are held in memory. Retaining unresolved evidence and all-hit mapping can be expensive on very large/repetitive libraries. Disk reduction is not a claim of constant-memory processing. No undocumented hit cap or truncated control library is used to improve speed.

## RBC post-analysis and review improvements — 2026-10-11

Run the separate **`csrna_rbc_diagnostics.py` after the main analysis**. It reads the main JSON plus the original paired manifest and selects only samples whose condition is RBC in both sources. It does not alter the core counts. See [RBC_DIAGNOSTICS.md](RBC_DIAGNOSTICS.md) for complete server commands, reference screens, output definitions, resume behavior, limitations, and tests.

The companion compares independent R1/R2 end-to-end and local alignments, retains a reproducible paired pilot, measures sequence complexity and homopolymers, and optionally screens secondary reference FASTAs. Counts, pair membership, and best hits are disk-backed in SQLite. Reference indexes are cached; inputs, configuration, source code, and retained pilots are checked before resuming. It requires `csrna_pipeline.py` beside it for validated input helpers, but importing that script does not execute main analysis. Full-library disk/time capacity remains to be measured on the actual server.

Local alignment boundaries are not native fragment ends. The companion uses raw reference IDs without the main pipeline's equivalence projections, so its parent counts must not be merged with strict endpoint counts. The current scope is independent-mate alignment, not native paired alignment or insert reconstruction. Origin matches are nonexclusive evidence, not contamination verdicts. Genome screens do not provide splice-aware or exon/intron annotation. Complexity labels never automatically remove reads, and no adapter/UMI protocol is inferred.

The dashboard adds **RBC diagnostics** as a separate report view. Main-report views now distinguish catalogue membership from recurrence, show singleton/one-library/multiple-library support tiers, provide parent-level replicate profiles, and flag high unmapped fractions, low eligible counts, and large within-condition mapping differences. Review thresholds are displayed and do not automatically exclude samples. Biological-replicate status is taken from metadata.

Main script version 2.1.2 changes documentation and renames the misleading performance field `mapped_query_count` to `queried_sequence_count`. Mapping, assignment, feature construction, and count normalization are unchanged. Legacy reports remain supported.

## Accuracy and performance audit — 2026-10-08

Version 2.1.1 records confirmed biological replication while preserving unknown preparation details. The [audit report](PERFORMANCE_ACCURACY_AUDIT.md) documents primary-source research, methodological decisions, benchmarks, and remaining limits.

Compiled input checks preserve accepted sequence alphabets and Phred+33 validation. On 100,000 existing real reads, matched counting/validation benchmarks improved 5.54×/6.14×; these are parsing measurements, not end-to-end speedups. Category matrices now avoid repeated feature scans. Reports and the frontend expose measured processing-stage durations, excluding final report assembly/serialization.

`--max-total-distinct-sequences` defaults to 500,000 across libraries, supplementing the 250,000 per-library guard. Exceeding either guard fails explicitly rather than discarding reads. These cardinality guards are not RAM guarantees; setting either to 0 disables that guard. Native paired alignment, insert reconstruction, and bounded-memory full-library processing remain unimplemented. R2 is validated; mapping uses R1 only.

Validation includes 59 backend tests, six actual Bowtie2/RNAfold workflows, and 16 frontend tests. A fresh 3,000-read run against the full combined reference reproduced all checked feature counts, category matrices, endpoint profiles, parent landscapes, support summaries, and family counts (574 features). This checks the existing local diagnostic subset, not the new server-resident batch. Identical biological profiles now correctly produce no PCA ordination. Biological replication does not establish native molecular ends, UMI identity, or absolute surface abundance.

## Replicate libraries and paired-file input — 2026-10-07

This extension supports the recent 15-library batch without designating other cell types as negative controls. **It is an R1-only analysis of paired libraries, not native paired-end alignment or insert assembly.** Both mates are validated across the entire input, but only R1 is aligned and counted. This preserves the existing Goal 1 spanning-read endpoint definition and Goal 2 shared feature table. Native molecular ends, biological replication, fragmentation, and UMI handling remain unconfirmed. No statistical significance, surface localization, or PCR correction is inferred.

### Input and analysis choices

- `--analysis-mode exploratory` permits any positive number of samples without target/controls. Target/control candidates and enrichment rankings are disabled, even for full libraries. All original mapping, ambiguity, exact/inexact, endpoint, folding, and category QC still run. Existing `target_control` behavior is the default for legacy commands.
- `--manifest FILE.csv --pair-handling r1-only` accepts a manifest with `sample_id,r1,r2,condition,replicate`. Sample IDs must be unique; each FASTQ path can occur only once. Relative paths resolve beside the manifest, independent of the terminal's directory. Gzipped FASTQ is read directly. `--manifest` cannot be mixed with `--sample` or `--sample-metadata`.
- Every R1/R2 record is checked for valid sequence and Phred+33 syntax, equal mate counts, and matching read IDs. `/1` and `/2` suffixes and Illumina mate fields are checked when present. Mismatches, truncation, invalid bases, and empty libraries fail the run before any new final report is published. Equal IDs without mate suffixes are accepted. Matching IDs do not establish the biological correctness of the sample sheet.
- `--pilot-pairs N --input-scope selected_subset` selects up to N R1 observations per library with uniform reservoir sampling across all validated pairs. It never takes just the first N reads. Sampling is deterministic for a given seed, sample ID, and file order. `--sampling-seed` defaults to 20261007. All pairs are still read and validated, so a small pilot reduces alignment work, not initial file scanning time. Scratch pilot FASTQs are removed with other temporary data. Without sampling, the original R1 is read directly after pair validation; no whole-library copy is made.
- `paired_input_qc` records source paths, validated pair totals, analyzed R1 counts, sampling policy, and seed. The denominator is the analyzed R1 population, never the sum of both mates and never the full-library total for a subset. `input_provenance` hashes both FASTQs and the manifest; inputs are checked again before publication.
- Optional manifest columns are retained as metadata. Use `replicate_type=biological`, `technical`, or `unknown`; do not infer type from numeric names. `pcr_cycles` must be an integer or blank. Additional useful columns are `library_preparation`, `fragmentation`, `trim_status`, `umi_status`, `strandedness`, and `batch`. PCR cycles are annotations, not a numerical abundance correction. Do not put sample indexes in a UMI field.
- `--report-detail summary` omits individual-read alignments/quality records, short-sequence records, and unused reference annotations/whole-reference alias and family maps. Feature counts, endpoint profiles, diagnostic summaries, and the shared matrices remain available. `report_omissions` explicitly records missing detail. Use `complete` for detailed read-level mismatch auditing. Both formats load in the frontend.
- `--quality-detail counts-only` validates FASTQ qualities but avoids per-position aggregate quality arrays; `positional` remains the default. It changes no alignment score or feature count. Missing quality summaries are explicitly unavailable, not evidence of good quality.
- Exploratory reports defer exhaustive pairwise rankings to the frontend rather than materializing every sample pair in JSON. Raw matrices and normalized tables remain available. Full-library descriptive sample comparisons are calculated on demand by the frontend. No biological fold ranking is enabled for a subset.

**Resource limit:** these options reduce avoidable memory/report overhead but do not convert mapping and feature construction into a bounded-memory engine. The per-library 250,000-distinct-sequence guard remains active, and the combined reference, union of sequences, and all-best-hit alignments also need RAM. Do not disable the guard simply to force a full-library run. A cluster with sufficient resources or a further disk-backed mapping/feature redesign is needed when the complete data exceed capacity. Running samples independently and joining their feature IDs is not equivalent to defining the shared features jointly.

### Run the recent batch on the Linux server

Copy the updated `reconstructed_pipeline` folder to the server and make the same combined reference available there. The new `samples_recent_batch.csv` names all 30 FASTQs supplied in the directory listing, with 15 distinct sample IDs and five conditions. All 15 libraries are confirmed independent biological samples. Preparation details remain unknown. The supplied PCR table supports 15 cycles for BMDC, BMDC_CpG, and M1, and 24 cycles for BMDM. RBC cycle values remain blank because the index-to-FASTQ mapping is unconfirmed; evidence is recorded in the manifest.

From the server's `sequencing_data` directory, set the actual project location and copy the manifest beside the FASTQs:

```bash
CSRNA_PROJECT=/absolute/path/to/wang-lab
cp "$CSRNA_PROJECT/reconstructed_pipeline/samples_recent_batch.csv" ./samples_recent_batch.csv
```

Edit metadata in that CSV when confirmed. Do not change mate paths unless the actual files differ. Confirm adapter/UMI preparation with the sequencing core: the pipeline does not trim adapters, remove UMIs, merge pairs, or reconstruct native ends. If preprocessing is needed, use protocol-appropriate prepared mates and update the paths.

Use Python 3.10+ with `bowtie2`, `bowtie2-build`, and `RNAfold` on PATH. An optional environment setup is:

```bash
conda create -n csrna -c conda-forge -c bioconda python=3.12 bowtie2=2.5.5 viennarna=2.7.2
conda activate csrna
python "$CSRNA_PROJECT/reconstructed_pipeline/csrna_pipeline.py" --self-test-tools
```

Start with 1,000 pairs sampled from each library (up to 15,000 R1 observations):

```bash
/usr/bin/time -v python "$CSRNA_PROJECT/reconstructed_pipeline/csrna_pipeline.py" \
  --manifest "$PWD/samples_recent_batch.csv" \
  --pair-handling r1-only \
  --analysis-mode exploratory \
  --input-scope selected_subset \
  --pilot-pairs 1000 \
  --sampling-seed 20261007 \
  --reference "$CSRNA_PROJECT/downstreamAnalysis/mouse_combined_with_rRNA_clean_dedup.fa" \
  --report-detail summary \
  --quality-detail counts-only \
  --top-fragments 1 --max-fold-variants 1 \
  --threads 4 \
  --output "$PWD/results/recent_batch_pilot.json" \
  2> "$PWD/recent_batch_pilot.log"
```

The small folding budgets limit structural predictions only; they do not filter mapping or feature counts. Inspect the log's failure/success status and maximum resident memory. A second pilot can use `--pilot-pairs 10000` and a new output name, subject to server capacity. Retain the same sampling seed for reproducibility. Optional diagnostic-reference and adapter flags described below still apply; absent screens remain `not_tested`. A high unmapped fraction cannot be identified as microbial, genomic, or adapter-derived from the fraction alone.

Only after resource and preparation review, a complete R1 analysis uses the same command with **both** `--pilot-pairs 1000` and `--input-scope selected_subset` replaced by `--pilot-pairs 0 --input-scope full_library`, with a new output filename. Increase folding budgets if more structural predictions are needed. An explicit memory-guard failure is not a successful full-library result. No new report is published on failure; existing outputs are protected unless `--overwrite` is explicitly given.

### Visualize and interpret

Copy the finished JSON to the Mac; the frontend does not need the FASTQs. From the project root:

```bash
/opt/anaconda3/bin/python reconstructed_pipeline/csrna_visualizer/run.py --report /absolute/path/to/recent_batch_pilot.json
```

The frontend remains a standalone read-only companion. Begin with **Overview** for unmapped and ambiguity QC; then **RNA composition** for all-category makeup. **Replicate QC** shows:

- Composition per individual library, arranged by condition, using all supplied reads.
- Pearson correlation of `log2(1 + RPM)` feature profiles, with features absent from every selected sample omitted. Constant profiles have undefined correlations and are blank. The selected RNA classes affect feature-based QC only.
- PCA of sample-by-feature `log2(1 + RPM)` values, centered per feature and not variance-scaled. Axis labels show variance explained; axis signs are arbitrary. It is descriptive ordination, not a significance test or batch correction.
- Per-condition feature recurrence across library columns. Detection means at least one eligible read; recurrence does not prove independent biological or molecular support. Technical replicates are not automatically pooled and unknown replicate types stay unknown.
- Sample metadata and pair/sampling audit tables. Inspect PCR differences and preparation information alongside potential outliers; do not automatically exclude a library because it differs.
- Raw integer feature-count CSV export for the selected evidence track and selected samples. All shared features are retained, even if the display's RNA-class filter excludes some. Use `sample_id` to match the metadata table. No pseudocount or RPM is added to raw counts.

Then use **Parent explorer** to compare observed endpoints and **Feature library** to inspect observed variants and available reference/variant folds. Goal 1 remains a shared definition of observed spanning-read endpoint features, with full-length-compatible/discrete/diffuse patterns; Goal 2 remains all-category shared counts and metadata. Subset displays describe only the sampled reads. Neither group separation, repeated endpoints, nor inexact support demonstrates surface specificity, degradation, modification, or biological function. Confirm biological replicate identities before downstream inferential statistics.


### Verification of this extension

57 backend unit/regression checks and six real Bowtie2/RNAfold workflows passed, including R1 equivalence, no-control mode, paired-input corruption, sampling, and count conservation. Fifteen frontend tests cover all six views, fifteen-library grouping, PCA distance preservation, constant/empty profiles, filtering, and exports. All six views also loaded the existing real 3,000-read/574-feature report. A fresh full-reference analysis of those reads reproduced the previous count matrices, endpoint profiles, parent landscape classifications, support tracks, family counts, and sample QC counts exactly. Notebook source and standalone documentation are synchronized and all code cells were executed sequentially with self-tests enabled.

The new server-resident 30 FASTQs were not available locally and have not been run here. Synthetic paired FASTQs test ingestion; the existing real R1 subset tests biological count preservation. This does not establish full-library resource capacity or native paired-end correctness. The user subsequently confirmed independent biological samples; the current manifest records this confirmation. Preparation details remain unknown. Native paired-end mapping/insert reconstruction and disk-backed whole-library processing remain future work.

## Goal 2 — shared RNA landscape comparison (2026-10-04)

The pipeline now packages its existing all-category shared features into a `landscape_comparison` report. It describes **observed-read endpoint features**; the molecular-end and full-library scale limitations documented below still apply. No new output files, alignment filters, candidate rules, or inferential statistics are introduced.

| Report field | Meaning |
|---|---|
| `feature_table` | One annotation row per shared feature: ID, parent RNA, RNA category, orientation, representative coordinates, endpoint clusters, and endpoint-evidence label |
| `raw_count_matrices_path` | Points to existing top-level `count_matrices`: unchanged integer exact, inexact, and inclusive matrices |
| `sample_metadata` | Metadata explicitly ordered to match every matrix’s sample columns; experimental fields are retained without inventing replicate identities |
| `normalized_rpm_matrices` | Exact, inexact, and inclusive counts normalized using each sample’s existing declared denominator |
| `normalization` | Denominator method, per-sample denominators, and the RPM pseudocount used only for ratios |
| `category_feature_count_matrices` | Eligible feature support summed by RNA category and sample; excludes ambiguous/ineligible reads, so it is not total RNA abundance |
| `pairwise_contrasts` | Every unordered pair of supplied samples, including control/control comparisons, with separate exact, inexact, and inclusive rankings |

The feature table and matrices have matching feature order and sample order. Counts include every shared feature observed in any sample; a feature found only in a control remains in the table with zero target counts. No RNA category is selectively dropped because it is not rRNA. RNA classes still depend on the supplied reference annotations, and excluded evidence remains in existing QC fields.

For a pair of samples A and B, each feature’s descriptive fold change is `(RPM_A + pseudocount_rpm) / (RPM_B + pseudocount_rpm)`. Rankings descend by signed log2 fold change. The configured target is the numerator whenever present; other pairs follow input sample order. Large positive values favor the numerator; negative values favor the denominator. Ties use larger observed RPM followed by feature ID for deterministic ordering. Counts and RPMs are included so a high ratio based on few reads is visible. Pseudocounts are not added to raw count matrices or RPM matrices.

Each ranked row labels detection as `both_samples`, `numerator_only`, or `denominator_only`, and direction as `higher_in_numerator`, `higher_in_denominator`, or `equal`. Zero/zero features stay in the matrices but are omitted from that pair’s track ranking. Rankings are descriptive: there are no p-values, adjusted p-values, significance labels, or new candidate-selection decisions. The existing target/control comparisons and legacy `fold_change_ranking` remain unchanged; use the new pairwise report when inspecting both directions and control-only features.

`selected_subset` runs provide the shared tables, RPM matrices, and category counts, but have empty pairwise contrasts and status `disabled_selected_subset`. Subset RPM describes only supplied reads and cannot support full-library biological fold comparisons. Full-library runs enable descriptive rankings without asserting that endpoints have been validated as molecular ends. The existing ingestion limit remains in force.

For a future replicate analysis, export the existing raw integer matrix and its aligned sample metadata. Specify genuine biological replicate, condition, and preparation/batch information when available. RPM matrices and pseudocount-adjusted ratios are not raw-count inputs. Endpoint features must be defined jointly across the compared samples; adding samples can change clusters and feature IDs, so independently generated feature tables should not be joined as if their boundaries were identical.

Example inspection after a run:

```python
import json
with open('csrna_results.json') as handle:
    report = json.load(handle)
landscape = report['landscape_comparison']
annotations = landscape['feature_table']
raw_counts = report['count_matrices']['inclusive']
metadata = landscape['sample_metadata']
assert raw_counts['sample_ids'] == [row['sample_id'] for row in metadata]
assert raw_counts['feature_ids'] == [row['id'] for row in annotations]
for contrast in landscape['pairwise_contrasts']:
    print(contrast['numerator_sample'], contrast['denominator_sample'])
    print(contrast['rankings']['inclusive'][:5])
```

Four additional regression tests cover cross-category feature union, normalization with unequal library sizes, both directions of comparison, control-only features, all sample pairs, metadata ordering, matrix validation, empty features/tracks, subset restrictions, and non-mutation. Real Bowtie2/RNAfold fixtures also validate the new report against actual mapped feature counts.

## Segment 1 — shared data contracts

`Hit` records a raw reference interval, orientation, CIGAR, edit distance, and independently checked exactness. `Mapping` retains all reported equal-best AS hits. DNA/RNA normalization and hashing helpers are shared by every stage. These contracts preserve raw evidence separately from the later biological assignment model.

## Segment 2 — input, qualities, and reference annotation

FASTA, FASTQ, and sequence-count parsers validate inputs and aggregate full supplied counts. FASTQ quality evidence is accumulated per normalized sequence and original read position. Full-library and explicitly selected-subset CSV modes are separate. Reference IDs must be unique and requested IDs must exist; a reference need not contain rRNA. Header/biotype classification supports multiple RNA categories, treats pseudogenes explicitly, and leaves unknown references unclassified instead of guessing a biological type.


`read_manifest`, `paired_r1_records`, and `prepare_manifest_inputs` validate paired FASTQs and prepare uniform pilot samples while preserving one R1 observation per pair. This is explicitly not paired alignment or merging. Relative paths resolve beside the manifest; all source inputs are hashed.

## Segment 3 — checked alignment and raw ambiguity

The index is built from precisely the selected reference. Bowtie2 maps deduplicated sequences once across all samples; stdout SAM is parsed incrementally and stderr retained in temporary logs to avoid pipe deadlocks. Query IDs, reference bounds, CIGAR spans, NM/AS tags, reverse orientation, and exact sequence agreement are checked. Raw tied hits remain available even after equivalence resolution. Bowtie2 is heuristic: uniqueness means uniqueness among reported best hits, not a proof against all possible alignments.

`unique_exact_rrna` is retained as a legacy regression helper; it is no longer the production eligibility gate.

## Segment 4 — coordinate primitives and compatibility helpers

`IntervalIndex` finds interval overlaps without recounting a sequence at several possible loci. `merge_intervals` merges positive overlaps only, separately by parent and orientation. `sequence_for_region` extracts the reference slice with orientation conversion; `rpm` rejects invalid denominators. These working primitives remain in use for coverage and overlap diagnostics.

The old asymmetric `compare_region` helper remains for historical regression checks only. Production feature abundance uses `compare_counts` below; it never calls the old helper.

## Segment 5 — folding and sequence metrics

RNAfold runs with `--noPS`, explicit temperature, timeout, and checked return status. Parsing verifies sequence identity, balanced dot-bracket structure, length, and finite energy, including space-padded zero values. Deliberate length skips use null energy and a status; malformed output fails. Duplicate sequences share cached folds within a run. GC fraction and longest homopolymer are descriptive sequence properties, not assay validation.

## Segment 6 — biological assignment and mismatch evidence

`reference_model` verifies full-sequence aliases and precursor/mature coordinate transforms. Family membership is separate from coordinate equivalence. `projected_hits` applies only transforms containing the complete aligned interval, with correct reverse coordinates. `resolve_mapping` separates raw locus certainty, resolved RNA/family identity, endpoint confidence, edit-rate acceptance, and exactness.

`alignment_differences` walks the CIGAR to preserve substitutions and indels. Substitutions include the original read position so the per-sample Phred summaries can be examined. No observed base is changed to the reference, and no mismatch is called a chemical modification.

## Segment 7 — shared fragment species and sample comparisons

`endpoint_clusters` uses bounded spans to prevent chain merging. `build_features` uses actual paired observations to construct one feature set across all samples, records sequence members, and counts exact/inexact support without preselecting high-RPM sequences. It keeps per-sample endpoint histograms. Representative sequences and coordinates come from a real supporting observation.

`compare_counts` applies the same count semantics to target and controls and disables inference in selected-subset mode. `parent_landscape` reports evidence-limited full-length, recurrent-fragment, or diffuse patterns without declaring degradation as established fact.

## Segment 8 — unresolved evidence and sample metadata

`sample_metadata` validates a one-to-one mapping from sample IDs to metadata rows while preserving condition/protocol fields. `unresolved_diagnostics` performs optional local secondary-reference and exact adapter checks only on primary-unmapped sequences. Sources and hashes are recorded. Evidence can match several diagnostic sources and is never promoted automatically to confirmed organism identity, contamination, or a revised primary count.

## Segment 9 — complete run and publication

`analyze` connects input validation, counting, mapping, projection, QC, feature construction, symmetric comparison, separate overlap screening, folding, unresolved diagnostics, and final accounting. Each sample is measured against the same union of endpoint features. All mapping-length sequences and all short-read counts remain in the report. Top-N limits only folding and presentation, within categories.

The matrices, metadata, family counts, endpoint patterns, observed/reference sequences, and coverage summaries are assembled into one schema-versioned report. All input files—including optional metadata/diagnostic references—are protected from being overwritten and hashed for provenance. `atomic_report` publishes only after a successful run. The parser and numeric validation make every threshold/policy explicit.


`landscape_comparison` validates shared matrix axes, integer count conservation, and metadata membership; builds annotation and RPM tables; aggregates eligible support by RNA category; and ranks descriptive pairwise sample contrasts. It consumes existing features without changing mapping, feature membership, thresholds, or candidate selection.

## Segment 10 — preserved regression checks

The original 17 deterministic checks still exercise input rejection, normalization, exactness, balanced-indel handling, strand extraction, interval overlap, control accounting, folding error handling, and safe publication. These protect working low-level behavior. Historical helper tests are not an assertion that the old asymmetric comparison remains the production analysis.

## Segment 11 — redesign checks

Additional tests cover non-rRNA references, category annotation, sequence-equivalence resolution, verified forward/reverse projections, palindrome ambiguity, family versus coordinate resolution, accepted mismatches, terminal-indel rejection, bounded endpoint clusters, paired features, count conservation, post-aggregation thresholds, symmetric ratios, selected-subset protection, Phred summaries, CIGAR mismatch positions, metadata validation, and cautious parent-pattern labels. Three additional tests cover read-weighted category ambiguity, inclusive threshold boundaries, empty/unambiguous categories, subset labeling, cross-category accounting, and non-mutating QC. Real-tool fixtures also verify alerts after reference equivalence resolution.

## Segment 12 — real-tool and paired-input checkpoints

`test_with_tools` runs six workflows with actual Bowtie2/RNAfold, including no-control paired-manifest R1 counting and uniform pilot sampling. It verifies that summary reports, omitted positional quality summaries, and R1-only manifest ingestion preserve count matrices, endpoint profiles, and parent patterns. `ReplicateInputTests` tests manifest identity, mate corruption, unequal lengths/counts, sampling reproducibility, full-file validation, metadata, and explicit mode requirements.

## Segment 13 — entry point and notebook safety

`main` validates options, checks tools, runs the selected tests or analysis in unique scratch space, and publishes the report. The script guard avoids treating Jupyter kernel arguments as pipeline arguments. The notebook defines the identical implementation before its final execution-control cell; production input paths are placeholders until explicitly supplied by the user.

## Reproducible verification

From the repository root:

```bash
python reconstructed_pipeline/csrna_pipeline.py --self-test
python reconstructed_pipeline/csrna_pipeline.py --self-test-tools
```

The first runs 50 unit/regression tests without external tools. The second repeats those and executes four real-tool synthetic workflows. These tests establish the checked behavior; they are not a mathematical proof that no possible input can reveal a bug.

The available repository sample data are selected Top100 rRNA tables, not full libraries. A real end-to-end validation run uses:

```bash
python reconstructed_pipeline/csrna_pipeline.py \
  --input-scope selected_subset \
  --sample M1=M1_longRNA.Top100_rRNA_sequences_RPM.csv \
  --sample BMDM=BMDM_longRNA.Top100_rRNA_sequences_RPM.csv \
  --sample BMDC=BMDC_trimmed.Top100_rRNA_sequences_RPM.csv \
  --reference downstreamAnalysis/mouse_combined_with_rRNA_clean_dedup.fa \
  --reference-id 'RIBO|BK000964.3' \
  --reference-id 'RIBO|NR_003278.3' \
  --reference-id 'RIBO|NR_003279.1' \
  --reference-id 'RIBO|NR_003280.1' \
  --reference-id 'RIBO|X00686.1' \
  --output /tmp/csrna_sample_validation.json
```

Use a fresh output path or explicitly add `--overwrite`. The reference restriction here tests the historical rRNA material; use an appropriate broader reference for the lab's multi-category production analysis.

With the documented defaults, the supplied sample run processes **229 distinct sequences and produces 83 shared fragment features**. Its read-weighted feature totals are:

| Supplied subset | Exact feature support | Inexact feature support | Inclusive feature support | Total supplied counts |
|---|---:|---:|---:|---:|
| M1 | 4,435,911 | 5,583,281 | 10,019,192 | 10,085,005 |
| BMDM | 19,223,516 | 14,749,891 | 33,973,407 | 34,117,094 |
| BMDC | 15,951,730 | 11,941,210 | 27,892,940 | 27,892,940 |

These are counts represented in selected tables, not full-library totals or independent molecule counts. Zero unmapped reads within these preselected rRNA tables cannot reproduce or contradict the supervisor's full-library unmapped fractions. No biological candidate or fold ranking is emitted in this mode.

Independent validation rebuilt counts directly from the CSVs, checked every sample CIGAR edit distance against the reference, verified feature membership and endpoint bounds, checked reference extraction and representative observed sequences, reconciled all three matrices, and checked fold lengths/finite energies. The notebook is also executed from a fresh kernel, and source/documentation identity is checked automatically during delivery preparation. Full raw-library and wet-lab validation remain separate requirements.

## Run the annotated notebook

Run all cells from top to bottom in a kernel with the required tools on PATH. The final cell defaults to `RUN_SELF_TESTS = True`, `RUN_TOOL_TESTS = False`, and `RUN_PIPELINE = False`. This defines the implementation and runs the deterministic checks without reading production libraries.

Enable `RUN_TOOL_TESTS` for the synthetic real-tool fixtures. For actual analysis, edit `PIPELINE_ARGUMENTS` with real paths and options and set `RUN_PIPELINE = True`. CSVs need the appropriate full-library assertion or explicit selected-subset mode. Nothing depends on an imported older copy of the script; the notebook contains all implementation segments. If editing the script later, keep the tagged source cells synchronized. Notebook metadata records the embedded source's SHA-256.

## Reasoning check

Two observations spanning `[100,150)` and `[140,190)` share coverage but have different endpoint pairs. They remain different features at the default tolerance, even though the secondary coverage summary can merge them. A control observation `[120,170)` may support neither feature but still appear in each feature's conservative overlap screen. This supports a distinction between fragment-species enrichment and the risk of a probe binding a shared region.

If a read maps equally to a mature 28S record and its sequence-identical region in an annotated precursor, a validated projection can resolve its RNA-level coordinates. It does not establish which precursor/mature molecule or genomic copy was sequenced. If the equivalent interval is palindromic and orientation remains unresolved, it is not forced into one strand's endpoint feature.

## Limitations and interpretation safeguards

- Header-based categories and families depend on supplied annotation. Unknowns remain labeled unknown; biological reference curation still matters.
- Shared endpoints and sequence variants define operational features under explicit tolerances. Protocol-specific clipping/trimming biases, PCR duplicates, or different library strandedness can affect those boundaries.
- Gapped/inexact reads retain observed sequence evidence; neither their edit burden nor standard RNAfold identifies modification chemistry.
- Family-level and feature-level counts answer different assignment questions and must not be added together.
- Primary-reference unmapped reads remain unexplained if no suitable diagnostic references/context are supplied. Optional diagnostics provide clues, not definitive contaminant labels.
- Candidate flags and ranked ratios are descriptive, especially with different extraction protocols and no replicate model. No universal normalization makes such designs directly comparable.
- Coverage summaries and parent-pattern labels do not prove intact long molecules or degradation mechanisms.
- Existing outputs are preserved on failure. Check exit status, timestamps, scope, and schema version before interpreting a result.

## References used for this reconstruction

The supervisor's `csrna_pipeline_0922.docx`, the current reconstructed implementation, `ProjectDocumentation.ipynb`, the original downstream mapping/structural scripts and visualization assumptions, supplied Top100/structural/contig tables, the BMDM SAM, and the combined reference FASTA informed this update. Statements in documents were assessed as reference material; executable changes follow the user's explicit authorization.

Tool behavior follows the official [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) and [RNAfold manual](https://viennarna.readthedocs.io/en/latest/man/RNAfold.html). The previously supplied project URL returned HTTP 404 during the earlier review; it was not used to invent additional biological requirements. The older featureCounts/edgeR/gene-translation branch remains a different gene-level analysis and is not silently mixed with sequence-level fragment counts.

## Goal 2 validation results (2026-10-04)

The extension passed 50 regression tests, four real Bowtie2/RNAfold synthetic workflows, and fresh execution of all 14 notebook code cells. A before/after run of the existing 229-sequence/83-feature rRNA fixture showed that every prior biological output was unchanged; only the new landscape report was added (timestamps and output paths differ).

A fresh real-data test selected 250 consecutive reads at each of four file positions (beginning, approximately 25%, 50%, and 75% of file bytes) from each supplied trimmed FASTQ: 1,000 reads per library, 3,000 total. This is a deterministic diagnostic subset, not a uniform random sample or a complete-library analysis. It ran against all 576,495 records of the combined reference, with `selected_subset` scope and folding budgets of one feature per category and one variant per feature. Folding budgets did not filter counts.

That run contained 2,314 distinct input sequences and produced 574 eligible observed-read features across six RNA categories. Independent checks reproduced input read counts, verified retained quality-record totals and paired endpoint membership, checked exact/inexact/inclusive aggregation, and confirmed feature/metadata axes, normalized RPM values, category-count conservation, and subset-ranking restrictions. All checks passed. Pairwise biological rankings were correctly disabled for this selected subset; full-library ranking behavior was tested in the synthetic complete-input workflows. This test does not establish whole-library feature counts, molecular fragment ends, or inferential significance.
