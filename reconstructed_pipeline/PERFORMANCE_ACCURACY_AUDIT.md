# Performance and accuracy audit — 2026-10-08

The pipeline is suitable for a controlled exploratory pilot, with tested improvements described below. It is **not certified as globally optimal, a bounded-memory full-library engine, or a native paired-end analyzer**. New server-resident FASTQs are not available in this local workspace; no result here is a measurement of their mapping or biological composition.

## Experimental facts now recorded

The user confirmed that the libraries are independent biological samples. All 15 rows in `samples_recent_batch.csv` now have `replicate_type=biological`, with the confirmation recorded in `replicate_evidence`. Other preparation details remain unknown. Historical reports retain their original metadata and are not rewritten.

The supplied PCR table supports 15 cycles for all BMDC, BMDC_CpG, and M1 libraries, and 24 cycles for all BMDM libraries. These group-consistent values are now recorded with their source. RBC has one 15-cycle library and two 24-cycle libraries, but the index-to-FASTQ correspondence is unconfirmed; its per-file cycle fields remain blank with an explanation. We do not invent a sample mapping, assign a UMI protocol from a UDI index name, or correct abundance by powers of two.

## Search and methodological decisions

| Area reviewed | Evidence and decision |
|---|---|
| Mapping accuracy vs speed | [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml): `-a` can be expensive, especially for repetitive references; capped reporting does not guarantee the retained alignments are the best possible. Preserve `--very-sensitive -a` and existing exact/inexact eligibility. Do not substitute a primary alignment or a small `-k` as proof of uniqueness. |
| Version-specific options | The installed, tested Bowtie2 is 2.5.5. Current online documentation includes newer deterministic-seed options; these have not been substituted into the established workflow. Options that change search behavior require separate sensitivity validation. |
| Paired input and trimming | [Cutadapt guide](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads): paired processing checks synchronization and paired filtering affects both mates. Preserve strict paired identity/count checks. No guessed adapters or automatic fixed-length trimming have been introduced; those could alter the endpoint evidence central to Goal 1. |
| Unknown strand protocol | [RSeQC documentation](https://rseqc.sourceforge.net/#infer-experiment-py): strandedness can be estimated from appropriately annotated alignments. It cannot be recovered from this MultiQC summary alone. The current transcript-reference alignment orientation remains labeled as read-to-reference orientation, not a verified biological strand. No strandedness estimate has been fabricated. |
| Duplication/PCR | [Parekh et al., 2016](https://www.nature.com/articles/srep25533): repeated RNA-seq reads can originate from amplification or distinct molecules. Preserve their abundance; do not deduplicate by sequence as a substitute for molecular UMI correction. Counting identical sequences once computationally while retaining their multiplicities is unchanged. |
| Replicate statistics | [DESeq2 vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): integer count matrices and explicit experimental designs are the appropriate starting point for replicate analysis. Keep raw counts separate from RPM. Biological replication is now confirmed, but formal differential testing is not implemented by this audit, and a pilot subset is not a full-library result. |
| Structure calculations | [RNAfold manual](https://www.tbi.univie.ac.at/RNA/RNAfold): parallel jobs increase memory requirements. Preserve per-sequence caching, explicit length/time limits, and separate reference/observed-variant folds rather than increasing concurrency without measurements. |

An important interpretation limit is the PCR difference between BMDM and the groups prepared with 15 cycles. For a BMDM-versus-M1 contrast, cycle count and condition coincide; the two effects cannot be independently estimated from that contrast alone. Other groups do not automatically remove that confounding. BMDC_CpG versus BMDC has matched reported cycle counts, but unknown preparation effects remain possible. Biological replicates support assessment of biological variability; they do not remove these design limitations.

## Implemented changes

1. **Compiled sequence/quality validation.** Replace Python character-by-character predicates with compiled patterns that enforce the same query/reference alphabets and printable Phred+33 range. Uppercase and U-to-T normalization, malformed-record checks, and positional-quality aggregation are unchanged. Exhaustive byte-range and selected Unicode tests check both acceptance and rejection.
2. **Whole-run cardinality guard.** `--max-total-distinct-sequences` defaults to 500,000, checked after each library is counted. This supplements the existing per-library 250,000 guard, preventing fifteen individually acceptable libraries from silently multiplying the mapping population without a whole-run check. It fails explicitly without discarding reads or publishing a report. `0` disables it, but disabling either guard is not recommended without resource assessment. Neither guard bounds reference, quality-array, alignment-multiplicity, or report memory.
3. **Less repeated category work.** Partition shared features by category once before constructing category matrices. All integer counts, category order, sample order, and normalization remain unchanged.
4. **Measured stage durations.** Reports now include `performance.stage_seconds`; matching messages are printed to stderr. The frontend exposes them under Report details. Stages cover input hashing/validation/counting, reference model/query union, index construction/alignment, assignments/features/folding, unresolved diagnostics, and read reporting/integrity recheck. Final report assembly and serialization are explicitly excluded. These are elapsed times, not per-process memory or server-capacity estimates.
5. **PCA numerical edge case.** Identical normalized sample profiles do not receive an ordination merely because floating-point centering leaves tiny residuals. Tests cover repeated identical biological profiles at several sample counts. Genuine nonidentical profiles still use the documented centered log2(1 + RPM) PCA.
6. **Confirmed experimental metadata.** Biological replication and the defensible PCR-table facts described above are recorded without changing unknown preparation fields.

## Matched ingestion benchmark

The benchmark used the first 100,000 records from the existing local M1 trimmed FASTQ. This is a computational workload, not a representative biological sample. Five alternating before/after trials used the same uncompressed file and Python 3.13.5 on macOS arm64. Median timings:

| Operation | Before | After | Relative speed |
|---|---:|---:|---:|
| Count reads without positional quality arrays | 0.7282 s | 0.1314 s | 5.54× |
| Validate FASTQ with names and normalized sequences | 0.7192 s | 0.1171 s | 6.14× |

All count dictionaries and validated record/length totals matched. The sorted sequence-count checksum was `12bfb006272da6ce624fbea65940f0b8c6ba8feace909ce58908aa82ca2401a1`.

These measurements exclude gzip decompression, pair-name comparison, hashing, mapping, and folding, and do not measure the default positional-quality mode. They are **not a 5–6× end-to-end speedup claim**. The `performance` fields and Linux `/usr/bin/time -v` on the actual server should identify the dominant stages for the new batch.

## Remaining performance and accuracy limits

- Reference indexing and all-best-hit alignment remain expensive. No unsafe hit cap, threshold relaxation, local clipping, or reference-category restriction was introduced to make a benchmark faster.
- Native paired alignment and insert reconstruction remain absent: R2 is validated but not used to improve mapping. This is the largest unresolved use-of-information limitation for the new batch. R1 results must not be presented as paired-insert or intact-native-RNA results.
- Read endpoint recurrence can reflect biological fragments, library fragmentation, or amplification. Independent biological replication strengthens reproducibility evidence but cannot identify the preparation mechanism by itself.
- A disk-backed full-library mapping/count/report design remains needed when combined unique sequences and alignment multiplicity exceed available memory. The new guard prevents an unsafe launch; it does not solve that architectural limit.
- No transcriptome-only pipeline can establish the origin of all unmapped reads. Genome/microbial/adapter diagnostic screens remain optional, explicitly labeled, and nonexclusive. Unknown screens remain unknown.
- RPM is descriptive relative abundance. It does not estimate absolute surface molecules or remove composition, amplification, extraction, or batch effects. Raw count export remains available for a separately validated replicate-statistics workflow.

## Use on the server

Copy the updated script and manifest to the server, keep the confirmed biological metadata, and retain unknown preparation fields. The documented pilot command in `README.md` still applies. Use `/usr/bin/time -v` around it and inspect both that output and the new stage timing messages. Start with 1,000 sampled pairs per library, then increase pilot size only as resources allow. Do not remove ambiguity-aware alignment settings or disable memory guards to force completion.

## Final verification — 2026-10-08

- 59 backend tests passed, including strict validation equivalence and whole-run guard failure before alignment.
- Six actual Bowtie2/RNAfold integration workflows passed.
- 16 frontend tests passed, including identical-profile PCA and biological metadata retention.
- The synchronized notebook executed its ordered code cells and all 59 tests. Its 13 source segments exactly reconstruct the script, its SHA-256 matches, and README.md exactly matches its Markdown annotations.
- A fresh full-reference run on the existing 3,000-read local diagnostic subset produced 574 features. Feature count matrices, category feature matrices, endpoint profiles, parent landscapes, support reporting, and family counts were exactly equal to the saved baseline. All six frontend views loaded the new report without exceptions.
- That run recorded 2,314 distinct mapped queries and about 256 seconds for index construction plus alignment. This is a single diagnostic workload, not a server or full-library throughput benchmark. Summary JSON was approximately 2.8 MB.

The browser preview uses this existing local diagnostic subset. It is not a result from the 15 newly confirmed biological libraries. No claim of exhaustive bug absence or full-batch validation is made.
