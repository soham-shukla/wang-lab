# RBC diagnostics after the main pipeline

`csrna_rbc_diagnostics.py` is a separate post-analysis command. Keep it beside `csrna_pipeline.py`, whose validated FASTQ/FASTA reading and header-classification helpers it imports. Importing the main script does not run the pipeline. The original analysis report and all raw inputs are read-only.

This program addresses the discrepancy between RBC pilot libraries by aligning **both mates independently**, in **end-to-end and local modes**, to the same combined reference. It adds sequence-complexity and optional origin-screen evidence. It does not trim, infer UMIs, reconstruct inserts, perform native paired alignment, establish contamination, or alter the main fragment counts.

## Server prerequisites and first run

Use the existing Conda environment with Python >=3.10, Bowtie2 and bowtie2-build. RNAfold is not needed by this post-analysis script. Copy the updated main script and new diagnostic script to `~/reconstructed_pipeline/`. Keep the original main report and the manifest beside the FASTQs as before.

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate "$HOME/envs/csrna"
cd ~/sequencing_data
mkdir -p results

/usr/bin/time -v python ~/reconstructed_pipeline/csrna_rbc_diagnostics.py \
  --analysis-report "$PWD/results/recent_batch_pilot.json" \
  --manifest "$PWD/samples_recent_batch.csv" \
  --reference "$HOME/reconstructed_pipeline/references/mouse_combined_with_rRNA_clean_dedup.fa" \
  --output-dir "$PWD/results/rbc_diagnostic_pilot" \
  --index-cache "$HOME/reconstructed_pipeline/reference_indexes" \
  --pilot-pairs 10000 \
  --sampling-seed 20261007 \
  --threads 4 \
  > "$PWD/results/rbc_diagnostic_pilot.stdout.log" \
  2> "$PWD/results/rbc_diagnostic_pilot.log"
echo "Exit code: $?"
```

Only manifest rows with **condition exactly `RBC`** are selected, and they must have matching RBC identities in the completed main report. Other cell FASTQs need not be accessible. Supplied mate checksums are verified against the main report when its provenance records them. The primary reference must have the same checksum. Do not relabel other cells as RBC to reuse this command.

The default pilot validates every pair and retains a uniform, reproducible paired subset. This scans the complete RBC inputs, not just the first reads. Both mates from a selected pair are retained with original quality strings and headers. The main report may be a pilot: its role is to establish sample identity and provenance, not to limit the post-analysis read population.

## Outputs and resume

The output directory contains:

- `rbc_diagnostics.json`: compact dashboard report, with scope, thresholds, provenance, read/mate counts, origin-screen status, and stage timings.
- `sequence_evidence.tsv.gz`: full per-distinct-sequence evidence for each library/mate, including multiplicities, mean quality, complexity, representative accepted alignments, clipping, and mode outcomes.
- `parent_counts.tsv`: per-library/mate/mode counts and RPM for assignments to a single raw reference parent. Ambiguous parent assignments are excluded. These are **not** the main pipeline's equivalence-resolved endpoint counts.
- `diagnostics.sqlite`: disk-backed pair membership, sequence counts and all best reported alignments. Inspect with SQLite; do not edit while a run is active.
- Paired pilot FASTQs (pilot mode only), query FASTA, configuration and aligner logs. These are intentional reusable intermediates, not hidden excess copies of the full raw libraries.

Repeat the exact command to resume. Input checksums and configuration must match. Completed ingestion and alignment stages are reused; summaries are rebuilt and the JSON is atomically published after integrity checks. Interrupted ingestion restarts that sample, interrupted alignment restarts that reference/mode. Changing settings or reference inputs requires a **new output directory**. An existing report is not silently deleted on a failed invocation.

A `.lock` prevents concurrent writers. After a hard process kill, verify no process is using that workspace/cache before manually removing a stale lock. Retained pilot FASTQs are checksum-checked on resume; corruption or deletion requires a new workspace. Resume still reads original inputs to verify checksums.

## Full RBC analysis

After checking the pilot and measured server disk/RAM/time, use the same command with:

```bash
--pilot-pairs 0 \
--output-dir "$PWD/results/rbc_diagnostic_full"
```

Use new log names too. A complete run examines every RBC pair even when the main report was a pilot. Counts and best-hit records are held in SQLite rather than a Python dictionary containing the whole library. SQLite disk usage and repetitive-reference alignment volume can still be large. Full-scale throughput and capacity have **not** been established locally for the server's approximately 86 million RBC pairs. Do not launch multiple processes against one output directory or remove ambiguity handling to force completion.

## Alignment evidence and thresholds

Both modes use Bowtie2 `-a` and retain all **reported best-score** loci on disk. The search is still heuristic. No cross-mode raw-score comparison is made. Diagnostic acceptance requires aligned length >=30 nt, aligned fraction >=0.5, and edit/aligned-length fraction <=0.1 by default:

```bash
--min-aligned-length 30 --min-aligned-fraction 0.5 --max-edit-rate 0.1
```

These are configurable diagnostic gates, not validated universal biological cutoffs. Short-RNA matches below these gates remain reported in SQLite but are not accepted. Inspect gate sensitivity and rejected matches before claiming an RNA class is absent. Best-score selection precedes these gates; the program does not promote a weaker hit just because the stronger hit fails a gate.

Mode outcomes are mutually exclusive: `both_same_parent`, `both_different_parent`, `both_parent_ambiguous`, `end_to_end_only`, `local_only`, and `neither_accepted`. The latter can contain weak reported alignments, not just entirely unmatched sequences. Counts must sum to each analyzed mate's read total. Local clipping is converted back to the original sequenced-read orientation. Detailed examples use one accepted best alignment; the SQLite table retains all ties.

For consistent sequence-count deduplication, alignment is performed on FASTA distinct sequences with multiplicities preserved. Base qualities do not affect alignment scoring; mean read Phred is a separate diagnostic statistic. A future quality-aware alignment policy must be validated separately. Mean clipped-base quality is not currently calculated.

Raw parent IDs are used conservatively. Reference aliases and precursor/mature projections from the main feature algorithm are not applied here. Independent mate agreement is **not** a concordance or insert-length test. R1 and R2 abundance must not be added and described as molecule counts.

## Complexity and origin screens

Metrics include read length, GC and N fraction, nucleotide entropy (ACGT probabilities over total length), longest homopolymer for each base, terminal runs on both ends, and distinct-3-mer fraction. Default labels use runs >=10 nt, entropy <1.2 bits, or distinct-3-mer fraction <0.1. Thresholds are descriptive and length-sensitive; no reads are removed. Read-weighted and distinct-sequence-weighted label percentages are separate and labels can overlap.

Optional secondary references are local FASTAs with unique labels:

```bash
--screen-reference genome=/absolute/path/to/mouse_genome.fa \
--screen-reference microbial=/absolute/path/to/microbial_reference.fa \
--screen-reference technical=/absolute/path/to/technical_reference.fa \
--adapter-sequence ACTUAL_CONFIRMED_ADAPTER_SEQUENCE
```

Replace the adapter placeholder with a confirmed A/C/G/T/N sequence of at least eight bases, or omit it. Do not guess preparation sequences. `other=FASTA` is also supported. Secondary references are screened in local mode against all analyzed sequences, allowing overlapping evidence and comparison with primary outcomes. Summary counts separately identify matches among `neither_accepted` reads. No automatic reference downloads occur. Adding screens changes configuration: choose a new output directory; existing reference indexes can be reused.

`not_tested` never means zero. Genome matches are unspliced sequence evidence, not exon/intron annotation or proof of DNA contamination. Microbial or technical matches need specificity, coverage, low-complexity, alternative-host-match and laboratory-context review. This program does not assign definitive taxonomic identity or contamination status. It does not yet perform splice-aware genome mapping or BED/GTF overlap annotation.

## Visualize and compare

Copy `rbc_diagnostics.json` to the Mac, select **RBC diagnostics** in the dashboard sidebar, and enter its local path or upload it. This view uses a separate report loader and does not need the main report to be loaded. It shows alignment-mode outcomes, complexity by outcome, overlapping label percentages, clipping, origin screens, paired agreement, and abundant sequence examples. The complete tables remain available in the server output directory.

The main dashboard also adds review flags (>=80% unmapped, <50 eligible reads, >=50 percentage-point within-condition mapping spread), parent-level replicate profiles, and observed feature support tiers. These are inspection aids, not automatic sample exclusions or inferential tests. Parent profiles in the **main** dashboard sum eligible endpoint support; diagnostic parent counts use a different raw-reference policy and must not be merged blindly.

Compare RBC with other cells using the main report's consistently processed endpoint/parent profiles. Do not compare full-depth diagnostic RBC results directly against other cells' 1,000-read pilots or mix local-rescued counts into strict endpoint matrices. A matched diagnostic analysis of other cell types would require a separately scoped extension; this script deliberately remains RBC-only.

## Validation commands

```bash
python ~/reconstructed_pipeline/csrna_pipeline.py --self-test-tools
CSRNA_TEST_TOOLS=1 python ~/reconstructed_pipeline/test_rbc_diagnostics.py
# On the Mac with the visualizer dependencies:
python reconstructed_pipeline/csrna_visualizer/test_visualizer.py
```

Tests cover mode partitions, real local-only recovery, raw-reference ambiguity, abundance conservation, malformed mates, clipping orientation, diagnostic gates, partial-run restart, pilot reproducibility and integrity, input provenance, and atomic publication. Existing local sequences and report-derived fixtures support software checks; artificial mates/qualities are never described as a new analysis of the original server RBC libraries.

## Implementation validation — 2026-10-11

- Main pipeline: 59 unit tests and six real-tool workflows passed. Mapping, assignment, and counting logic are unchanged; the version and documentation were updated, and the performance field `mapped_query_count` was corrected to `queried_sequence_count` because it includes unmapped queries.
- RBC companion: 10 tests passed, including real Bowtie2 end-to-end/local alignment, local-only recovery, an optional technical-reference screen, ambiguous reference parents, and restart/integrity checks.
- Dashboard: 20 tests passed. All seven pages additionally opened with the existing 15-library pilot report; the RBC page loaded the diagnostic fixture, rendered its plots, and handled an empty sample selection.
- Notebook: all code cells executed with production analysis disabled; 59 tests passed, and embedded source, source checksum, and standalone documentation matched.
- Existing-data check: 112 test pairs were constructed from RBC sequences in the existing report with artificial reverse-complement mates and artificial qualities, using 14 corresponding reference entries. Both alignment modes completed and conserved read/pair counts. This is a software fixture, not an analysis of original paired reads or evidence about library preparation.

The complete original RBC FASTQs were not available locally. Full-library throughput, storage requirements, and biological conclusions remain unverified until the server run. Genome/microbial screens also require the relevant reference FASTAs. Run configurations record code hashes as well as input hashes; changed code or settings require a new output directory rather than silently reusing old alignments.

## Preparation documentation update — 2026-10-12

The lab's `RNA_Sequencing_Documentation.docx` reports Watchmaker mRNA Library Prep Kit used without RNA selection, 1 ng RNA input, strand-specific construction, unique dual indexes, and Fragment Analyzer assessment. The batch manifest now records this provenance. Fragmentation is **unconfirmed**, read-to-RNA orientation remains **unverified**, and adapter trimming and UMI status remain unknown. UDI sample indexes do not establish molecular identifiers. The phrase “150 nt” does not resolve cycles per mate; verify FASTQ lengths and run settings.

The document reports 24 PCR cycles, conflicting with earlier sample-specific values for several groups. `pcr_cycles` retains the earlier recorded values (RBC remains blank pending sample-specific confirmation); `pcr_cycles_document_reported` records 24 separately. `pcr_cycles_review_status` marks conflicting, consistent, or unresolved sample-specific records. Neither record is silently preferred for count correction. Fragment Analyzer assessment is recorded as a procedure, not proof that each sample passed a DNA-contamination test.

Endpoint clusters remain descriptive observed-read features. Library fragmentation and cDNA synthesis can create boundaries; even a reconstructed library insert need not be a native RNA molecule. Endpoint peaks, discrete/diffuse patterns, and predicted fragment folds do not establish native csRNA species or surface localization. Prioritize RNA-category composition, parent support, coverage, and replicate consistency while preparation details are unresolved. No strand filtering, trimming, UMI inference, automatic deduplication, or numerical PCR correction was introduced.

**Applying the update:** copy the updated manifest to the server and merge its annotation columns into your working manifest while preserving correct FASTQ paths and any newer verified metadata. New main reports include the supplied metadata and preparation warnings. Existing JSON reports remain unchanged and retain their original metadata. For RBC analysis, the current manifest supplies diagnostic annotations, while `main_analysis_sample_metadata` preserves the original main-report snapshot. Sample identities and reference/input integrity checks remain enforced. Updated code or metadata requires a new RBC output directory; existing completed reports are not rewritten. The dashboard explains endpoint limitations even for older reports and exposes current versus original RBC metadata explicitly.

## Main counting checkpoints

The main script now supports resumable SQLite counting. Only its completed schema-2 **analysis** JSON may be supplied to this extension. Count-only and alignment-capacity JSON files are diagnostic artifacts, not main reports. The RBC extension keeps its own workspace; main count caches do not contain the R2 sequence evidence needed by RBC diagnostics. Updated main-script source changes RBC workspace provenance, so use a new RBC output directory after installing this version. See [main workflow](README.md#resumable-counting-and-alignment-capacity-workflow--2026-10-13).
