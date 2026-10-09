#!/usr/bin/env python3
"""Multi-category csRNA fragment discovery; see README.md. Python >=3.10, standard library.

No statistical significance, surface localization, or probe specificity is inferred.
The only permanent analysis output is an atomically published JSON report.
"""

# %% Imports and data contracts
import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import csv
import random
import time
from itertools import zip_longest

VERSION = "2.1.1"
DNA = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
READ_BASES = re.compile(r'[ACGTN]+')
REFERENCE_BASES = re.compile(r'[ACGTRYSWKMBDHVN]+')
INVALID_PHRED33 = re.compile(r'[^!-~]')


def reverse_complement(sequence):
    return sequence.translate(DNA)[::-1]


def clean_sequence(value, reference=False):
    sequence = value.strip().upper().replace("U", "T")
    alphabet = REFERENCE_BASES if reference else READ_BASES
    if not alphabet.fullmatch(sequence):
        raise ValueError("Empty sequence or unsupported bases/embedded whitespace")
    return sequence


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).lower().endswith(".gz") else open(path, encoding="utf-8")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True, order=True)
class Hit:
    reference: str
    start: int
    end: int
    strand: str
    cigar: str
    edits: int
    exact: bool


@dataclass
class Mapping:
    score: int | None = None
    hits: set = field(default_factory=set)

    def add(self, score, hit):
        if self.score is None or score > self.score:
            self.score, self.hits = score, {hit}
        elif score == self.score:
            self.hits.add(hit)


# %% Input parsing and reference annotation
def fasta_records(handle):
    header, pieces = None, []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(pieces)
            header, pieces = line[1:].strip(), []
            if not header:
                raise ValueError("Empty FASTA header")
        elif header is None:
            raise ValueError("FASTA sequence before its header")
        else:
            pieces.append(line)
    if header is not None:
        yield header, "".join(pieces)


def fastq_records(handle, with_quality=False, with_name=False):
    """Strict four-line, single-end FASTQ; malformed records fail, never truncate."""
    record = 0
    while True:
        header = handle.readline()
        if not header:
            return
        record += 1
        sequence, plus, quality = [handle.readline() for _ in range(3)]
        if not header.startswith("@") or not header[1:].strip() or not plus.startswith("+") or not quality:
            raise ValueError(f"Malformed/truncated FASTQ record {record}")
        sequence, quality = sequence.rstrip("\r\n"), quality.rstrip("\r\n")
        if sequence != sequence.strip():
            raise ValueError(f'Whitespace in FASTQ sequence at record {record}')
        if not sequence or len(sequence) != len(quality) or INVALID_PHRED33.search(quality):
            raise ValueError(f"Invalid sequence/quality in FASTQ record {record}")
        if plus[1:].strip() and plus[1:].strip() != header[1:].strip():
            raise ValueError(f"FASTQ + header differs at record {record}")
        if with_name:
            yield header[1:].strip(), clean_sequence(sequence), quality
        else:
            yield (sequence, quality) if with_quality else sequence


def read_manifest(path):
    """Explicit R1/R2 manifest; relative paths are relative to the manifest."""
    path = Path(path).expanduser().resolve(strict=True)
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        if not {'sample_id', 'r1', 'r2', 'condition', 'replicate'} <= set(columns) or len(columns) != len(set(columns)):
            raise ValueError('Manifest needs unique sample_id,r1,r2,condition,replicate columns')
        rows = list(reader)
    if not rows:
        raise ValueError('Empty sample manifest')
    seen, files = set(), set()
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError('Malformed manifest row')
        row.update({k: v.strip() for k, v in row.items()})
        name = row['sample_id']
        if not re.fullmatch(r'[A-Za-z0-9_-]+', name) or name in seen or not row['condition']:
            raise ValueError('Manifest needs unique valid sample IDs and nonempty conditions')
        seen.add(name)
        for mate in ('r1', 'r2'):
            value = Path(row[mate]).expanduser()
            if not row[mate] or not str(value).lower().removesuffix('.gz').endswith(('.fq', '.fastq')):
                raise ValueError('Manifest mates must be FASTQ paths')
            value = (value if value.is_absolute() else path.parent / value).resolve(strict=True)
            if value in files:
                raise ValueError('A FASTQ file cannot represent multiple samples or both mates')
            files.add(value)
            row[mate] = str(value)
        kind = row.setdefault('replicate_type', 'unknown') or 'unknown'
        if kind not in ('biological', 'technical', 'unknown'):
            raise ValueError('replicate_type must be biological, technical, or unknown')
        row['replicate_type'] = kind
        if row.get('pcr_cycles') and not re.fullmatch(r'\d+', row['pcr_cycles']):
            raise ValueError('pcr_cycles must be a nonnegative integer or blank')
    return rows


def paired_r1_records(r1, r2):
    """Validate every mate and identity; emit one R1 observation per pair, never merge."""
    def identity(header, expected):
        parts = header.split()
        token = parts[0]
        if token.endswith(('/1', '/2')):
            if token[-1] != str(expected):
                raise ValueError('FASTQ mate suffix contradicts R1/R2 position')
            token = token[:-2]
        if len(parts) > 1 and re.match(r'^[12]:', parts[1]) and parts[1][0] != str(expected):
            raise ValueError('FASTQ Illumina mate field contradicts R1/R2 position')
        return token
    with open_text(r1) as one, open_text(r2) as two:
        for index, pair in enumerate(zip_longest(fastq_records(one, with_name=True), fastq_records(two, with_name=True)), 1):
            a, b = pair
            if a is None or b is None:
                raise ValueError(f'Mate record counts differ at pair {index}')
            if identity(a[0], 1) != identity(b[0], 2):
                raise ValueError(f'Mate identifiers differ at pair {index}')
            yield a[1], a[2]


def prepare_manifest_inputs(rows, args, workdir):
    """Stream validation and optional uniform reservoir; scratch FASTQs are removed on exit."""
    paths, audit = {}, {}
    for row in rows:
        sample = row['sample_id']
        # Stable per-sample seed: reordering the manifest cannot change sampling.
        seed = int(hashlib.sha256(f'{args.sampling_seed}:{sample}'.encode()).hexdigest(), 16)
        rng, reservoir, total = random.Random(seed), [], 0
        destination = Path(workdir) / f'{sample}.r1.fastq'
        def write_record(handle, index, record):
            sequence, quality = record
            handle.write(f'@{sample}_{index}\n{sequence}\n+\n{quality}\n')
        for total, record in enumerate(paired_r1_records(row['r1'], row['r2']), 1):
            if total % 1000000 == 0:
                print(f'PAIR QC — {sample}: checked {total:,} pairs…', file=sys.stderr)
            if not args.pilot_pairs:
                continue
            if total <= args.pilot_pairs:
                reservoir.append((total, record))
            else:
                slot = rng.randrange(total)
                if slot < args.pilot_pairs:
                    reservoir[slot] = (total, record)
        if not total:
            raise ValueError(f'Empty paired library: {sample}')
        if args.pilot_pairs:
            with open(destination, 'w', encoding='utf-8') as handle:
                for index, record in sorted(reservoir):
                    write_record(handle, index, record)
        used = min(total, args.pilot_pairs) if args.pilot_pairs else total
        paths[sample] = destination if args.pilot_pairs else Path(row['r1'])
        audit[sample] = {'r1': row['r1'], 'r2': row['r2'], 'validated_pairs': total,
                         'analyzed_r1_reads': used, 'sampling_seed': args.sampling_seed if args.pilot_pairs else None,
                         'sampling': 'uniform_reservoir_without_replacement' if args.pilot_pairs else 'all_pairs',
                         'pair_handling': 'r1_only_r2_validated', 'r2_used_for_alignment': False}
        print(f'PAIR QC — {sample}: validated {total:,} pairs; analyzing {used:,} R1 reads', file=sys.stderr)
    return paths, audit


def load_counts(path, counts_complete=False, allow_selected=False, qualities=None, max_distinct_sequences=0):
    """Count full libraries BEFORE mapping or abundance filtering."""
    counts = Counter()
    def add_count(sequence, count):
        if max_distinct_sequences and sequence not in counts and len(counts) >= max_distinct_sequences:
            raise ValueError(f'Distinct-sequence limit ({max_distinct_sequences:,}) reached in {path}. '
                             'No reads were silently discarded. This in-memory pipeline requires a larger-memory '
                             'execution environment or a disk-backed implementation for this library; '
                             'raise --max-distinct-sequences only after assessing memory requirements.')
        counts[sequence] += count
    name = str(path).lower().removesuffix(".gz")
    with open_text(path) as handle:
        if name.endswith(".csv"):
            if not counts_complete and not allow_selected:
                raise ValueError("CSV requires --counts-complete: all sequences in a full library, never top-N tables")
            reader = csv.DictReader(handle)
            columns = reader.fieldnames or []
            if len(columns) != len(set(columns)):
                raise ValueError("Duplicate CSV column names")
            if not allow_selected and ("RPM_total_mapped" in columns or "RPM_within_rRNA" in columns):
                raise ValueError("Legacy rRNA RPM tables are selected subsets, not complete raw-count inputs")
            count_columns = [c for c in ("count", "counts") if c in columns]
            if "sequence" not in columns or len(count_columns) != 1:
                raise ValueError("CSV needs sequence and exactly one of count/counts")
            for row in reader:
                if None in row or row["sequence"] is None:
                    raise ValueError("Malformed CSV row")
                value = row[count_columns[0]]
                if value is None or not re.fullmatch(r"\d+", value.strip()):
                    raise ValueError("Counts must be nonnegative integers")
                sequence = clean_sequence(row["sequence"])
                count = int(value)
                if count:
                    add_count(sequence, count)
        elif name.endswith((".fastq", ".fq")):
            for sequence, quality in fastq_records(handle, with_quality=True):
                sequence = clean_sequence(sequence)
                add_count(sequence, 1)
                if qualities is not None:
                    # Avoid allocating and discarding two arrays on every repeated read.
                    if sequence not in qualities:
                        qualities[sequence] = {"records": 0, "phred_sum": [0] * len(sequence),
                                               "below_q20": [0] * len(sequence)}
                    values = qualities[sequence]
                    values["records"] += 1
                    for index, char in enumerate(quality):
                        q = ord(char) - 33
                        values["phred_sum"][index] += q
                        values["below_q20"][index] += q < 20
        elif name.endswith((".fasta", ".fa", ".fna")):
            for _, sequence in fasta_records(handle):
                add_count(clean_sequence(sequence), 1)
        else:
            raise ValueError(f"Unsupported input extension: {path}")
    if not counts:
        raise ValueError(f"No positive read counts: {path}")
    return counts


def classify_reference(header):
    """Heuristic annotation of supplied reference headers, not a gene annotation service."""
    h = header.lower()
    if any(x in h for x in ("mitochond", "chrmt", "|mt-", "mtrnr", "mtrna")):
        return "mitochondrial_RNA"
    if "pseudogene" in h:
        return "pseudogene"
    if h.startswith("ribo|") or re.search(r'(^|[|\s_])rrna($|[|\s_])', h) or "ribosomal rna" in h:
        return "rRNA"
    if "snorna" in h:
        return "snoRNA"
    if "snrna" in h:
        return "snRNA"
    if "trna" in h:
        return "tRNA"
    if "mirna" in h or re.search(r"(^|[|\s])(?:mir-|let-7)", h):
        return "miRNA"
    if "pirna" in h:
        return "piRNA"
    if "protein_coding" in h or "mrna" in h:
        return "mRNA"
    if any(x in h for x in ("lncrna", "lincrna", "antisense", "processed_transcript", "sense_intronic", "sense_overlapping", "non_coding")):
        return "lncRNA"
    if any(x in h for x in ("y_rna", "yrna", "vault", "7sl", "7sk", "srp", "scarna", "misc_rna")):
        return "structured_small_ncRNA"
    return "unclassified_reference"


def load_reference(path, selected_ids=()):
    references, annotations = {}, {}
    selected = set(selected_ids)
    seen = set()
    with open_text(path) as handle:
        for header, sequence in fasta_records(handle):
            name = header.split()[0]
            if name in seen:
                raise ValueError(f"Duplicate reference ID: {name}")
            seen.add(name)
            if selected and name not in selected:
                continue
            references[name] = clean_sequence(sequence, reference=True)
            annotations[name] = {"header": header, "category": classify_reference(header)}
    if selected - references.keys():
        raise ValueError(f"Missing requested reference IDs: {sorted(selected - references.keys())}")
    if not references:
        raise ValueError("Empty reference selection")
    return references, annotations


# %% Checked alignment, streamed SAM, and explicit ambiguity
def require_tools():
    tools = {}
    for name in ("bowtie2", "bowtie2-build", "RNAfold"):
        executable = shutil.which(name)
        if executable is None:
            raise RuntimeError(f"Required executable not on PATH: {name}")
        result = subprocess.run([executable, "--version"], capture_output=True, text=True, check=True)
        tools[name] = {"path": executable, "version": (result.stdout or result.stderr).strip()}
    return tools


def checked_to_log(command, workdir):
    with tempfile.TemporaryFile(mode="w+t", dir=workdir) as log:
        result = subprocess.run(command, cwd=workdir, stdout=log, stderr=log)
        if result.returncode:
            log.seek(0)
            raise RuntimeError(f"Command failed ({result.returncode}): {command[0]}\n{log.read()[-12000:]}")


def parse_sam_record(line, sequences, references):
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11 or not re.fullmatch(r"q\d+", fields[0]):
        raise ValueError("Malformed SAM record/query name")
    query = int(fields[0][1:])
    if not 0 <= query < len(sequences):
        raise ValueError("Unexpected query in SAM")
    flag = int(fields[1])
    if flag & 4:
        return query, None, None
    if flag & 2048:
        raise ValueError("Unexpected supplementary alignment in end-to-end mapping")
    ref, start, cigar = fields[2], int(fields[3]) - 1, fields[5]
    if ref not in references:
        raise ValueError("SAM uses a reference outside the constructed index")
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if "".join(n + op for n, op in ops) != cigar or not ops or any(int(n) <= 0 for n, _ in ops):
        raise ValueError("Invalid CIGAR")
    span = sum(int(n) for n, op in ops if op in "MDN=X")
    query_span = sum(int(n) for n, op in ops if op in "MIS=X")
    end = start + span
    if start < 0 or end > len(references[ref]) or span <= 0 or query_span != len(sequences[query]):
        raise ValueError("SAM alignment outside reference or inconsistent query length")
    tags = dict(tag.split(":", 2)[::2] for tag in fields[11:] if tag.startswith(("AS:i:", "NM:i:")))
    if not {"AS", "NM"} <= tags.keys():
        raise ValueError("Mapped SAM record lacks AS/NM tags")
    edits, score = int(tags["NM"]), int(tags["AS"])
    if edits < 0:
        raise ValueError('Negative NM edit distance in SAM')
    strand = "-" if flag & 16 else "+"
    sequence = sequences[query] if strand == "+" else reverse_complement(sequences[query])
    ungapped = all(op in "M=" for _, op in ops)
    exact = ungapped and edits == 0 and sequence == references[ref][start:end] and set(sequence) <= set("ACGT")
    return query, score, Hit(ref, start, end, strand, cigar, edits, exact)


def align_sequences(sequences, references, tools, threads, workdir):
    """Map distinct sequences once across all samples, retaining all best-score hits."""
    workdir = str(Path(workdir).resolve())
    mappings = [Mapping() for _ in sequences]
    if not sequences:
        return mappings
    fasta = Path(workdir) / "reference.fa"
    queries = Path(workdir) / "queries.fa"
    prefix = str(Path(workdir) / "index")
    with fasta.open("w") as handle:
        for name, sequence in references.items():
            handle.write(f">{name}\n{sequence}\n")
    with queries.open("w") as handle:
        for query, sequence in enumerate(sequences):
            handle.write(f">q{query}\n{sequence}\n")
    checked_to_log([tools["bowtie2-build"]["path"], "--threads", str(threads), str(fasta), prefix], workdir)
    command = [tools["bowtie2"]["path"], "-x", prefix, "-f", "-U", str(queries),
               "--end-to-end", "--very-sensitive", "-a", "--seed", "0", "-p", str(threads)]
    # stderr goes to a temporary file so a full pipe cannot deadlock SAM streaming.
    with tempfile.TemporaryFile(mode="w+t", dir=workdir) as log:
        process = subprocess.Popen(command, cwd=workdir, stdout=subprocess.PIPE, stderr=log, text=True)
        seen = set()
        try:
            for line in process.stdout:
                if line.startswith("@"):
                    continue
                query, score, hit = parse_sam_record(line, sequences, references)
                seen.add(query)
                if hit is not None:
                    mappings[query].add(score, hit)
            status = process.wait()
            if status:
                log.seek(0)
                raise RuntimeError(f"Bowtie2 failed ({status}):\n{log.read()[-12000:]}")
            if len(seen) != len(sequences):
                raise RuntimeError("Bowtie2 output missing query records")
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.kill()
            process.wait()
    return mappings


def mapping_category(mapping, annotations):
    categories = {annotations[hit.reference]["category"] for hit in mapping.hits}
    if not categories:
        return "unmapped"
    return next(iter(categories)) if len(categories) == 1 else "multi_category_ambiguous"


def unique_exact_rrna(mapping, annotations):
    if len(mapping.hits) != 1:
        return None
    hit = next(iter(mapping.hits))
    return hit if hit.exact and annotations[hit.reference]["category"] == "rRNA" else None


# %% Strand-aware intervals, library normalization, and controls
class IntervalIndex:
    """Query overlap without rescanning all reads; each sequence ID returned once."""
    def __init__(self, entries):
        groups = defaultdict(list)
        for query, hit in entries:
            groups[(hit.reference, hit.strand)].append((hit.start, hit.end, query))
        self.groups = {}
        for key, values in groups.items():
            values.sort()
            starts, maxima, maximum = [], [], -1
            for start, end, _ in values:
                starts.append(start)
                maximum = max(maximum, end)
                maxima.append(maximum)
            self.groups[key] = (values, starts, maxima)

    def overlap(self, reference, strand, start, end):
        if end <= start:
            raise ValueError("Invalid interval")
        values, starts, maxima = self.groups.get((reference, strand), ([], [], []))
        position = bisect_left(starts, end) - 1
        found = set()
        while position >= 0 and maxima[position] > start:
            if values[position][1] > start:
                found.add(values[position][2])
            position -= 1
        return found


def rpm(count, denominator):
    if denominator <= 0:
        raise ValueError("RPM denominator must be positive")
    return count * 1_000_000.0 / denominator


def merge_intervals(entries):
    """Half-open coordinates; join only positive overlaps, never mere adjacency."""
    regions = []
    for query, hit in sorted(entries, key=lambda x: (x[1].reference, x[1].strand, x[1].start, x[1].end, x[0])):
        if (regions and regions[-1]["parent_rRNA"] == hit.reference and
                regions[-1]["strand"] == hit.strand and hit.start < regions[-1]["end"]):
            regions[-1]["end"] = max(regions[-1]["end"], hit.end)
            regions[-1]["member_ids"].append(query)
        else:
            regions.append({"parent_rRNA": hit.reference, "strand": hit.strand,
                            "start": hit.start, "end": hit.end, "member_ids": [query]})
    return regions


def sequence_for_region(region, references):
    sequence = references[region["parent_rRNA"]][region["start"]:region["end"]]
    if region["strand"] == "-":
        sequence = reverse_complement(sequence)
    if len(sequence) != region["end"] - region["start"]:
        raise ValueError("Reference extraction length mismatch")
    return sequence.replace("T", "U")


def compare_region(region, owner, controls, counts_by_id, denominators, exact_index, all_index, args):
    # Target support excludes ambiguous/inexact mappings. Control screening includes
    # every best-score locus, even ambiguous/inexact ones, counted once per region.
    support_ids = exact_index.overlap(region["parent_rRNA"], region["strand"], region["start"], region["end"])
    control_ids = all_index.overlap(region["parent_rRNA"], region["strand"], region["start"], region["end"])
    if args.control_strand == "both":
        opposite = "-" if region["strand"] == "+" else "+"
        control_ids |= all_index.overlap(region["parent_rRNA"], opposite, region["start"], region["end"])
    support = sum(counts_by_id[owner].get(query, 0) for query in support_ids)
    abundance = rpm(support, denominators[owner])
    comparisons = {}
    reasons = []
    if abundance < args.min_rpm:
        reasons.append("below_min_rpm")
    for control in controls:
        count = sum(counts_by_id[control].get(query, 0) for query in control_ids)
        control_rpm = rpm(count, denominators[control])
        enrichment = (abundance + args.pseudocount_rpm) / (control_rpm + args.pseudocount_rpm)
        comparisons[control] = {"overlap_count": count, "overlap_rpm": control_rpm,
                                "fold_enrichment": enrichment, "observed_absent": count == 0}
        if enrichment < args.min_enrichment:
            reasons.append(f"below_enrichment:{control}")
        if args.require_control_absence and count:
            reasons.append(f"observed_in_control:{control}")
    return {"support_count": support, "support_rpm": abundance,
            "controls": comparisons, "abundance_control_pass": not reasons,
            "rejection_reasons": reasons}


# %% RNAfold parsing and descriptive sequence characterization
FOLD_LINE = re.compile(r"^([.()]+)\s+\(\s*([-+]?\d+(?:\.\d+)?)\s*\)\s*$")


def parse_fold(stdout, sequence):
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if len(lines) != 2 or lines[0].upper().replace("T", "U") != sequence.upper().replace("T", "U"):
        raise ValueError("Unexpected RNAfold output/sequence")
    match = FOLD_LINE.fullmatch(lines[1])
    if match is None or len(match[1]) != len(sequence):
        raise ValueError("Invalid RNAfold structure/energy")
    balance = 0
    for symbol in match[1]:
        balance += (symbol == "(") - (symbol == ")")
        if balance < 0:
            raise ValueError("Unbalanced RNAfold structure")
    if balance:
        raise ValueError("Unbalanced RNAfold structure")
    energy = float(match[2])
    if not math.isfinite(energy):
        raise ValueError("Non-finite RNAfold energy")
    return {"status": "ok", "structure": match[1], "mfe_kcal_mol": energy,
            "mfe_per_nt": energy / len(sequence),
            "unpaired_fraction_in_mfe_structure": match[1].count(".") / len(sequence)}


def fold_sequence(sequence, tools, temperature, max_length, timeout, workdir, cache):
    if sequence in cache:
        return cache[sequence]
    if len(sequence) > max_length:
        result = {"status": "skipped_length_limit", "structure": None, "mfe_kcal_mol": None}
    else:
        process = subprocess.run([tools["RNAfold"]["path"], "--noPS", "-T", str(temperature)],
                                 input=sequence + "\n", capture_output=True, text=True,
                                 cwd=workdir, timeout=timeout)
        if process.returncode:
            raise RuntimeError(f"RNAfold failed ({process.returncode}): {process.stderr[-4000:]}")
        result = parse_fold(process.stdout, sequence)
    cache[sequence] = result
    return result


def sequence_metrics(sequence):
    longest = max((len(match[0]) for match in re.finditer(r"(.)\1*", sequence)), default=0)
    return {"length_nt": len(sequence), "gc_fraction": (sequence.count("G") + sequence.count("C")) / len(sequence),
            "longest_homopolymer": longest}


# %% Reference equivalence and per-read evidence

def rrna_kind(header):
    match = re.search(r'(?<![\w.])(5\.8s|18s|28s|5s)(?![\w.])', header.lower())
    return match[1] if match else None


def reference_model(references, annotations, args):
    """Collapse identical full sequences; project only validated identical intervals."""
    groups = defaultdict(list)
    for name, sequence in references.items():
        groups[(annotations[name]['category'], hashlib.sha256(sequence.encode()).hexdigest())].append(name)
    aliases = {}
    for names in groups.values():
        canonical = min(names)
        for name in names:
            # Hash equality alone is not a sequence-equivalence proof.
            if references[name] != references[canonical]:
                raise ValueError('Reference digest collision')
            aliases[name] = canonical
    projections = []
    families = {}
    for name in references:
        header = annotations[name]['header']
        kind = rrna_kind(header) if annotations[name]['category'] == 'rRNA' else None
        gene = re.search(r'\bENS[A-Z]*G\d+', header)
        families[name] = ('rRNA:' + kind if kind else 'gene:' + gene[0] if gene else 'parent:' + aliases[name])
    for specification in args.parent_family:
        name, sep, family = specification.partition('=')
        if not sep or name not in references or not family:
            raise ValueError('Use --parent-family REFERENCE_ID=FAMILY with an existing reference')
        families[name] = family
    for names in groups.values():
        if len({families[n] for n in names}) != 1:
            if any(spec.split('=', 1)[0] in names for spec in args.parent_family):
                raise ValueError('Identical references have conflicting explicit family annotations')
            for name in names:
                families[name] = 'sequence_equivalence:' + aliases[name]
    def add_projection(source, target, offset, strand, origin):
        if source not in references or target not in references or source == target:
            raise ValueError('Invalid source/target in reference projection')
        if annotations[source]['category'] != annotations[target]['category']:
            raise ValueError('Reference projection crosses RNA categories')
        expected = references[target] if strand == '+' else reverse_complement(references[target])
        if offset < 0 or references[source][offset:offset + len(expected)] != expected:
            raise ValueError('Reference projection is not an exact full-target sequence match')
        if aliases[source] == aliases[target] and offset == 0 and strand == '+':
            return
        row = {'source': source, 'target': target, 'offset': offset, 'strand': strand,
               'length': len(expected), 'origin': origin}
        if row not in projections:
            projections.append(row)
    if not args.no_auto_rrna_projections:
        precursors = [n for n in references if annotations[n]['category'] == 'rRNA' and
                      any(t in annotations[n]['header'].lower() for t in ('repeating unit', 'pre-ribosomal', 'preribosomal'))]
        mature = [n for n in references if annotations[n]['category'] == 'rRNA' and rrna_kind(annotations[n]['header']) and n not in precursors]
        for source in precursors:
            for target in mature:
                for strand, sequence in [('+', references[target]), ('-', reverse_complement(references[target]))]:
                    offset = references[source].find(sequence)
                    while offset >= 0:
                        add_projection(source, target, offset, strand, 'annotated_rRNA_exact_containment')
                        offset = references[source].find(sequence, offset + 1)
    for specification in args.reference_projection:
        try:
            source, target_spec = specification.split('=', 1)
            target, offset, strand = target_spec.rsplit(':', 2)
            if strand not in ('+', '-'):
                raise ValueError()
            offset = int(offset)
        except ValueError as error:
            raise ValueError('Use --reference-projection SOURCE=TARGET:ZERO_BASED_OFFSET:+/-') from error
        add_projection(source, target, offset, strand, 'explicit_verified_projection')
    source_ids = {aliases[p['source']] for p in projections}
    if any(aliases[p['target']] in source_ids for p in projections):
        raise ValueError('Chained/cyclic reference projections are unsupported; project directly to terminal parents')
    by_source = defaultdict(list)
    for row in projections:
        by_source[aliases[row['source']]].append(row)
    return {'aliases': aliases, 'families': families, 'projections': projections, 'by_source': by_source}


def projected_hits(hit, model):
    result = set()
    for projection in model['by_source'].get(model['aliases'][hit.reference], []):
        offset, length = projection['offset'], projection['length']
        if offset <= hit.start and hit.end <= offset + length:
            start, end, strand = hit.start - offset, hit.end - offset, hit.strand
            if projection['strand'] == '-':
                start, end = length - end, length - start
                strand = '-' if strand == '+' else '+'
            result.add((model['aliases'][projection['target']], start, end, strand))
    return result or {(model['aliases'][hit.reference], hit.start, hit.end, hit.strand)}


def resolve_mapping(mapping, model, annotations, length, max_edit_rate):
    locations = set().union(*(projected_hits(hit, model) for hit in mapping.hits)) if mapping.hits else set()
    parent_families = {model['families'][location[0]] for location in locations}
    category = mapping_category(mapping, annotations)
    if not locations:
        status = 'unmapped'
    elif len(locations) == 1:
        status = 'unique_locus' if len(mapping.hits) == 1 else 'equivalent_parent_resolved'
    elif len(parent_families) == 1:
        status = 'family_resolved_coordinates_ambiguous'
    else:
        status = 'unresolved_parents'
    # No clipping, terminal indels, or reference skips may define a measured endpoint.
    boundaries = bool(mapping.hits) and all(
        re.fullmatch(r'(?:\d+[MID=X])+', h.cigar) and
        re.match(r'^\d+[M=X]', h.cigar) and re.search(r'\d+[M=X]$', h.cigar)
        for h in mapping.hits)
    edit_ok = bool(mapping.hits) and all(h.edits >= 0 and h.edits / length <= max_edit_rate for h in mapping.hits)
    exact = bool(mapping.hits) and all(h.exact for h in mapping.hits)
    eligible = len(locations) == 1 and boundaries and edit_ok
    reasons = []
    if len(locations) != 1:
        reasons.append(status)
    if mapping.hits and not boundaries:
        reasons.append('uncertain_endpoints')
    if mapping.hits and not edit_ok:
        reasons.append('edit_rate_above_limit')
    return {'status': status, 'category': category, 'locations': sorted(locations),
            'family': next(iter(parent_families)) if len(parent_families) == 1 else None,
            'boundary_confident': boundaries, 'edit_rate_accepted': edit_ok,
            'exact': exact, 'feature_eligible': eligible, 'exclusion_reasons': reasons}


def alignment_differences(sequence, hit, references):
    """CIGAR-aware mismatches, with original-read positions for Phred diagnostics."""
    oriented = sequence if hit.strand == '+' else reverse_complement(sequence)
    query, reference = 0, hit.start
    events = []
    for count, op in re.findall(r'(\d+)([MIDNSHP=X])', hit.cigar):
        count = int(count)
        if op in 'M=X':
            for offset in range(count):
                if oriented[query + offset] != references[hit.reference][reference + offset]:
                    position = query + offset if hit.strand == '+' else len(sequence) - 1 - query - offset
                    events.append({'type': 'substitution', 'reference_position': reference + offset,
                                   'read_position': position, 'observed_oriented_base': oriented[query + offset],
                                   'reference_base': references[hit.reference][reference + offset]})
            query += count
            reference += count
        elif op == 'I':
            events.append({'type': 'insertion', 'reference_position': reference, 'length': count,
                           'oriented_query_offset': query, 'sequence': oriented[query:query + count]})
            query += count
        elif op in 'DN':
            events.append({'type': 'deletion' if op == 'D' else 'reference_skip',
                           'reference_position': reference, 'length': count})
            reference += count
        elif op == 'S':
            query += count
    return events


# %% Shared endpoint features and symmetric count matrices

def endpoint_clusters(coordinates, tolerance):
    """Bounded diameter prevents single-linkage chaining across a degradation smear."""
    assignment, clusters = {}, []
    for coordinate in sorted(set(coordinates)):
        if not clusters or coordinate - clusters[-1]['minimum'] > 2 * tolerance:
            clusters.append({'minimum': coordinate, 'maximum': coordinate})
        else:
            clusters[-1]['maximum'] = coordinate
        assignment[coordinate] = len(clusters) - 1
    for cluster in clusters:
        cluster['center'] = (cluster['minimum'] + cluster['maximum']) // 2
    return assignment, clusters


def build_features(resolved, sequences, counts_by_id, references, annotations, tolerance):
    groups = defaultdict(list)
    for query, evidence in enumerate(resolved):
        if evidence['feature_eligible']:
            parent, start, end, strand = evidence['locations'][0]
            groups[(parent, strand)].append((query, start, end))
    features, endpoints = [], []
    for (parent, strand), rows in sorted(groups.items()):
        start_map, starts = endpoint_clusters([s for _, s, _ in rows], tolerance)
        end_map, ends = endpoint_clusters([e for _, _, e in rows], tolerance)
        paired = defaultdict(list)
        for query, start, end in rows:
            paired[(start_map[start], end_map[end])].append(query)
        histograms = {}
        for sample, counts in counts_by_id.items():
            start_counts, end_counts = Counter(), Counter()
            for query, start, end in rows:
                if counts.get(query, 0):
                    start_counts[start] += counts[query]
                    end_counts[end] += counts[query]
            histograms[sample] = {'reference_start': dict(sorted(start_counts.items())),
                                  'reference_end': dict(sorted(end_counts.items()))}
        endpoints.append({'parent_RNA': parent, 'strand': strand,
                          'reference_start_clusters': starts, 'reference_end_clusters': ends,
                          'histograms_by_sample': histograms})
        for (start_cluster, end_cluster), queries in sorted(paired.items()):
            # Select a real observed pair, never a fabricated combination of marginal modes.
            representative = min(queries, key=lambda q: (-sum(c.get(q, 0) for c in counts_by_id.values()), sequences[q]))
            _, start, end, _ = resolved[representative]['locations'][0]
            identity = [parent, strand, starts[start_cluster], ends[end_cluster]]
            identifier = 'fragment_' + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
            row = {'id': identifier, 'parent_RNA': parent, 'category': annotations[parent]['category'], 'strand': strand,
                   'endpoint_evidence': 'sequenced_read_boundaries; molecular ends require library/insert-length evidence',
                   'start': start, 'end': end, 'start_cluster': starts[start_cluster], 'end_cluster': ends[end_cluster],
                   'five_prime_cluster': starts[start_cluster] if strand == '+' else ends[end_cluster],
                   'three_prime_cluster': ends[end_cluster] if strand == '+' else starts[start_cluster],
                   'member_ids': sorted(queries), 'representative_read_id': representative,
                   'observed_sequence': sequences[representative].replace('T', 'U'),
                   'reference_sequence': sequence_for_region({'parent_rRNA': parent, 'start': start, 'end': end, 'strand': strand}, references)}
            row['counts'] = {}
            for sample, counts in counts_by_id.items():
                exact = sum(counts.get(q, 0) for q in queries if resolved[q]['exact'])
                inexact = sum(counts.get(q, 0) for q in queries if not resolved[q]['exact'])
                row['counts'][sample] = {'exact': exact, 'inexact': inexact, 'inclusive': exact + inexact}
            features.append(row)
    return sorted(features, key=lambda row: row['id']), endpoints


def compare_counts(counts, denominators, target, controls, args):
    if getattr(args, 'analysis_mode', 'target_control') == 'exploratory':
        return {'status': 'disabled_no_controls', 'candidate': None, 'rejection_reasons': ['exploratory_analysis']}
    if args.input_scope == 'selected_subset':
        return {'status': 'disabled_selected_subset', 'candidate': None, 'rejection_reasons': ['incomplete_libraries']}
    abundance = rpm(counts[target], denominators[target])
    reasons = ['below_min_rpm'] if abundance < args.min_rpm else []
    comparisons = {}
    for control in controls:
        control_rpm = rpm(counts[control], denominators[control])
        enrichment = (abundance + args.pseudocount_rpm) / (control_rpm + args.pseudocount_rpm)
        comparisons[control] = {'count': counts[control], 'rpm': control_rpm, 'fold_enrichment': enrichment,
                                'log2_fold_enrichment': math.log2(enrichment), 'observed_absent': counts[control] == 0}
        if enrichment < args.min_enrichment:
            reasons.append('below_enrichment:' + control)
        if args.require_control_absence and counts[control]:
            reasons.append('observed_in_control:' + control)
    return {'status': 'descriptive_no_replicate_statistics', 'target_count': counts[target], 'target_rpm': abundance,
            'controls': comparisons, 'candidate': not reasons, 'rejection_reasons': reasons,
            'minimum_control_fold_enrichment': min(v['fold_enrichment'] for v in comparisons.values())}


def parent_landscape(features, references, samples, args, resolved=None, counts_by_id=None):
    groups = defaultdict(list)
    for feature in features:
        groups[(feature['parent_RNA'], feature['strand'])].append(feature)
    rows = []
    for (parent, strand), members in sorted(groups.items()):
        for sample in samples:
            values = [f['counts'][sample]['inclusive'] for f in members]
            total = sum(values)
            if not total:
                continue
            full = sum(value for f, value in zip(members, values)
                       if f['start'] <= args.endpoint_tolerance and len(references[parent]) - f['end'] <= args.endpoint_tolerance)
            if resolved is not None and counts_by_id is not None:
                # Full-parent evidence must come from this sample's actual reads,
                # not the pooled feature representative selected in another sample.
                full = sum(counts_by_id[sample].get(query, 0)
                           for feature in members for query in feature['member_ids']
                           if resolved[query]['locations'][0][1] <= args.endpoint_tolerance
                           and len(references[parent]) - resolved[query]['locations'][0][2] <= args.endpoint_tolerance)
            concentration = max(values) / total
            recurrent = sum(value for value in values if value >= args.min_species_reads)
            if total < args.min_landscape_reads:
                label = 'insufficient_support'
            elif full / total >= args.landscape_dominance:
                label = 'full_length_compatible'
            elif recurrent / total >= args.landscape_dominance:
                label = 'discrete_fragment_pattern'
            else:
                label = 'diffuse_or_multiple_fragments'
            rows.append({'parent_RNA': parent, 'strand': strand, 'sample': sample, 'eligible_read_count': total,
                         'feature_count': sum(v > 0 for v in values), 'largest_feature_fraction': concentration,
                         'recurrent_feature_support_fraction': recurrent / total,
                         'full_length_compatible_fraction': full / total, 'descriptive_pattern': label,
                         'interpretation': 'Diffuse support is smear-compatible, not proof of degradation; duplicate reads are not independent molecules.'})
    return rows


# %% Optional unresolved-read diagnostics and metadata

def sample_metadata(path, samples):
    if path is None:
        return [{'sample_id': sample, 'condition': sample, 'replicate': None, 'library_preparation': 'unspecified'} for sample in samples]
    with open_text(Path(path).expanduser()) as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or 'sample_id' not in reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Sample metadata needs unique headers and sample_id')
        rows = list(reader)
    if any(None in row or not row['sample_id'] for row in rows) or len({r['sample_id'] for r in rows}) != len(rows):
        raise ValueError('Malformed/duplicate sample metadata rows')
    if {r['sample_id'] for r in rows} != set(samples):
        raise ValueError('Sample metadata must cover exactly the input samples')
    by_name = {row['sample_id']: row for row in rows}
    return [by_name[name] for name in samples]


def diagnostic_categories(args):
    labels = [spec.partition('=')[0] for spec in args.diagnostic_reference]
    categories = {}
    for spec in args.diagnostic_reference_category:
        label, sep, category = spec.partition('=')
        if not sep or label not in labels or label in categories or category not in ('genome', 'microbial', 'other'):
            raise ValueError('Use unique --diagnostic-reference-category LABEL=genome|microbial|other for a supplied reference')
        categories[label] = category
    return {label: categories.get(label, 'other') for label in labels}


def unmapped_qc(diagnostics, counts_by_id, summaries, reference_categories, adapter_configured):
    """Non-exclusive evidence labels plus a disjoint read-count partition; no reassignment."""
    configured = set(reference_categories.values()) | ({'adapter'} if adapter_configured else set())
    read_labels = {}
    for query, evidence in diagnostics.items():
        labels = set()
        for item in evidence:
            if item['type'] == 'adapter_exact_substring':
                labels.add('adapter')
            elif item['type'] == 'secondary_reference_match':
                labels.add(reference_categories[item['reference_label']])
        read_labels[query] = sorted(labels)
    reports = {}
    for sample, summary in summaries.items():
        total, unmapped = summary['total_input_reads'], summary['unmapped_reads']
        partition, evidence_counts = Counter(), Counter()
        for query, labels in read_labels.items():
            count = counts_by_id[sample].get(query, 0)
            status = ('multiple_source_evidence' if len(labels) > 1 else labels[0] + '_evidence' if labels
                      else 'no_match_in_supplied_screens' if configured else 'not_tested')
            partition[status] += count
            for label in labels:
                evidence_counts[label] += count
        if sum(partition.values()) != unmapped:
            raise RuntimeError('Unmapped QC read accounting checkpoint failed')
        reports[sample] = {
            'input_scope': summary['input_scope'], 'primary_unmapped_reads': unmapped,
            'total_input_reads': total, 'primary_unmapped_fraction_of_input': unmapped / total,
            'below_min_read_length_reads': summary['below_min_read_length_reads'],
            'partition_counts': {k: v for k, v in sorted(partition.items()) if v},
            'screens': {label: {
                'status': ('not_tested' if label not in configured else 'no_unmapped_reads' if not unmapped else 'tested'),
                'matching_reads': evidence_counts[label] if label in configured else None,
                'fraction_of_unmapped': evidence_counts[label] / unmapped if label in configured and unmapped else None,
            } for label in ('genome', 'microbial', 'adapter', 'other')},
            'identity_status': 'no_primary_unmapped_reads' if not unmapped else 'unresolved_origin',
        }
    return reports, read_labels


def unresolved_diagnostics(sequences, mappings, args, tools, workdir):
    indices = [q for q, m in enumerate(mappings) if not m.hits]
    evidence = {q: [] for q in indices}
    sources = []
    adapters = [clean_sequence(sequence) for sequence in args.adapter_sequence]
    if any(len(a) < 8 for a in adapters):
        raise ValueError('Adapter diagnostic sequences must be at least 8 bases')
    for query in indices:
        for adapter in adapters:
            if adapter in sequences[query] or reverse_complement(adapter) in sequences[query]:
                evidence[query].append({'type': 'adapter_exact_substring', 'adapter': adapter})
    names = set()
    for index, specification in enumerate(args.diagnostic_reference):
        name, sep, raw_path = specification.partition('=')
        if not sep or not name or name in names:
            raise ValueError('Use unique --diagnostic-reference LABEL=FASTA')
        names.add(name)
        path = Path(raw_path).expanduser().resolve(strict=True)
        before = sha256(path)
        references, annotations = load_reference(path)
        with tempfile.TemporaryDirectory(prefix=f'diagnostic-{index}-', dir=workdir) as scratch:
            diagnostic_mappings = align_sequences([sequences[q] for q in indices], references, tools, args.threads, scratch)
        if sha256(path) != before:
            raise RuntimeError('Diagnostic reference changed during analysis')
        sources.append({'label': name, 'path': str(path), 'sha256': before})
        for query, mapping in zip(indices, diagnostic_mappings):
            if mapping.hits:
                evidence[query].append({'type': 'secondary_reference_match', 'reference_label': name,
                                        'best_alignment_score': mapping.score,
                                        'best_alignments': [h.__dict__.copy() for h in sorted(mapping.hits)]})
    return evidence, sources


# %% Orchestration, accounting checkpoints, and one atomic report

def category_ambiguity_qc(sample, retention, threshold, input_scope):
    """Read-weighted, descriptive QC only; never changes mapping or eligibility."""
    metrics, alerts = {}, []
    for category, counts in sorted(retention.items()):
        total = counts.get('mapped_reads', 0)
        ambiguous = counts.get('raw_locus_ambiguous_reads', 0)
        eligible = counts.get('feature_eligible_reads', 0)
        fraction = ambiguous / total if total else None
        triggered = ambiguous > 0 and fraction >= threshold
        metrics[category] = {
            'mapped_reads': total, 'raw_locus_ambiguous_reads': ambiguous,
            'raw_locus_ambiguous_fraction': fraction,
            'feature_eligible_reads': eligible,
            'feature_eligible_fraction': eligible / total if total else None,
            'warning_threshold': threshold, 'warning_triggered': triggered,
        }
        if triggered:
            scope = 'SELECTED SUBSET ONLY' if input_scope == 'selected_subset' else 'supplied library'
            message = (
                f'HIGH MAPPING AMBIGUITY — {sample}/{category} [{scope}]: '
                f'{ambiguous:,}/{total:,} mapped reads ({fraction:.1%}) have multiple raw best-score loci '
                f'(warning threshold {threshold:.1%}). '
                f'{eligible:,}/{total:,} ({eligible / total:.1%}) retain resolved endpoint support '
                'after equivalence and eligibility checks. Raw ambiguity can reflect redundant references; '
                'excluded support is not necessarily caused by ambiguity. Inspect category_retention '
                'and assignment_state_counts before interpreting fragment results.')
            alerts.append({'code': 'HIGH_CATEGORY_MAPPING_AMBIGUITY', 'severity': 'warning',
                           'sample': sample, 'category': category, 'input_scope': input_scope,
                           **metrics[category], 'message': message})
    return metrics, alerts


def support_reporting(features, resolved, counts_by_id, summaries):
    """Report disjoint support tracks without changing feature membership or counts."""
    labels = ('raw_unique_exact', 'raw_unique_inexact',
              'equivalence_resolved_exact', 'equivalence_resolved_inexact')
    sample_ids = list(counts_by_id)
    matrices = {label: {'feature_ids': [f['id'] for f in features], 'sample_ids': sample_ids,
                        'data': []} for label in labels}
    by_sample = {s: defaultdict(Counter) for s in sample_ids}
    for feature in features:
        rows = {label: [] for label in labels}
        for sample in sample_ids:
            counts = Counter()
            for query in feature['member_ids']:
                evidence = resolved[query]
                if not evidence['feature_eligible'] or evidence['status'] not in ('unique_locus', 'equivalent_parent_resolved'):
                    raise RuntimeError('Support reporting encountered an ineligible feature member')
                prefix = 'raw_unique' if evidence['status'] == 'unique_locus' else 'equivalence_resolved'
                counts[prefix + ('_exact' if evidence['exact'] else '_inexact')] += counts_by_id[sample].get(query, 0)
            exact = counts['raw_unique_exact'] + counts['equivalence_resolved_exact']
            inexact = counts['raw_unique_inexact'] + counts['equivalence_resolved_inexact']
            if feature['counts'][sample] != {'exact': exact, 'inexact': inexact, 'inclusive': exact + inexact}:
                raise RuntimeError('Support reporting feature-count checkpoint failed')
            for label in labels:
                rows[label].append(counts[label])
            by_sample[sample][feature['category']].update(counts)
        for label in labels:
            matrices[label]['data'].append(rows[label])

    def describe(counts, denominator):
        values = {label: counts[label] for label in labels}
        values['exact_only'] = values['raw_unique_exact'] + values['equivalence_resolved_exact']
        values['accepted_inexact'] = values['raw_unique_inexact'] + values['equivalence_resolved_inexact']
        values['inclusive'] = values['exact_only'] + values['accepted_inexact']
        return {'counts': values, 'rpm': {k: rpm(v, denominator) for k, v in values.items()},
                'inexact_fraction_of_eligible_support': values['accepted_inexact'] / values['inclusive'] if values['inclusive'] else None,
                'inclusive_to_exact_support_ratio': values['inclusive'] / values['exact_only'] if values['exact_only'] else None}

    output = {}
    for sample, categories in by_sample.items():
        # Include mapped RNA categories even when all support failed feature eligibility.
        for category in summaries[sample]['category_retention']:
            categories[category]
        total = Counter()
        for counts in categories.values():
            total.update(counts)
        denominator = summaries[sample]['rpm_denominator']
        output[sample] = {'input_scope': summaries[sample]['input_scope'], 'rpm_denominator': denominator,
                          **describe(total, denominator),
                          'by_category': {c: describe(n, denominator) for c, n in sorted(categories.items())}}
    return {'samples': output, 'assignment_count_matrices': matrices,
            'definitions': {
                'exact_only': 'Feature-eligible exact support; existing count_matrices.exact',
                'accepted_inexact': 'Feature-eligible non-exact support; existing count_matrices.inexact',
                'inclusive': 'Exact-only plus accepted-inexact; existing count_matrices.inclusive',
                'raw_unique': 'One raw best-score alignment, with existing endpoint/edit-rate eligibility checks',
                'equivalence_resolved': 'Multiple raw alignments resolved to one oriented parent interval by verified reference equivalence',
                'interpretation': 'Inexact includes accepted substitutions/indels; it is sensitivity evidence, not proof of RNA modification. All tracks use the same shared features and RPM denominator; overlapping totals must not be added together.'}}


def landscape_comparison(features, matrices, metadata, denominators, args):
    """Goal 2: aligned feature/metadata tables and descriptive sample contrasts."""
    tracks = ('exact', 'inexact', 'inclusive')
    ids = [f['id'] for f in features]
    samples = list(denominators)
    if len(set(ids)) != len(ids) or len({m['sample_id'] for m in metadata}) != len(metadata):
        raise RuntimeError('Landscape table contains duplicate feature/sample identifiers')
    metadata_by_id = {m['sample_id']: m for m in metadata}
    if set(metadata_by_id) != set(samples):
        raise RuntimeError('Landscape metadata does not match matrix samples')
    for track in tracks:
        matrix = matrices[track]
        if matrix['feature_ids'] != ids or matrix['sample_ids'] != samples or len(matrix['data']) != len(ids):
            raise RuntimeError('Landscape matrix axes are inconsistent')
        for feature, row in zip(features, matrix['data']):
            if len(row) != len(samples) or any(type(n) is not int or n < 0 for n in row):
                raise RuntimeError('Landscape counts must be nonnegative integers')
            if row != [feature['counts'][sample][track] for sample in samples]:
                raise RuntimeError('Landscape matrix disagrees with feature counts')
    if any(matrices['exact']['data'][i][j] + matrices['inexact']['data'][i][j] != matrices['inclusive']['data'][i][j]
           for i in range(len(ids)) for j in range(len(samples))):
        raise RuntimeError('Landscape exact/inexact conservation checkpoint failed')
    normalized = {track: {'feature_ids': ids, 'sample_ids': samples,
                         'data': [[rpm(n, denominators[sample]) for sample, n in zip(samples, row)]
                                  for row in matrices[track]['data']]} for track in tracks}
    annotations = [{key: f[key] for key in ('id', 'parent_RNA', 'category', 'strand', 'start', 'end',
                                           'start_cluster', 'end_cluster', 'endpoint_evidence')}
                   for f in features]
    categories = sorted({f['category'] for f in features})
    features_by_category = defaultdict(list)
    for feature in features:
        features_by_category[feature['category']].append(feature)
    category_counts = {track: {'category_ids': categories, 'sample_ids': samples,
                              'data': [[sum(f['counts'][s][track] for f in features_by_category[category])
                                        for s in samples] for category in categories]} for track in tracks}
    contrasts = []
    # Target is the numerator in its contrasts; other sample pairs follow input order.
    ordered = ([args.target] if args.target in samples else []) + [s for s in samples if s != args.target]
    defer_contrasts = getattr(args, 'analysis_mode', 'target_control') == 'exploratory'
    if args.input_scope == 'full_library' and not defer_contrasts:
        for i, numerator in enumerate(ordered):
            for denominator in ordered[i + 1:]:
                rankings = {}
                for track in tracks:
                    rows = []
                    for f in features:
                        a, b = f['counts'][numerator][track], f['counts'][denominator][track]
                        if not (a or b):
                            continue  # Both-zero rows remain in the shared matrices.
                        ar, br = rpm(a, denominators[numerator]), rpm(b, denominators[denominator])
                        ratio = (ar + args.pseudocount_rpm) / (br + args.pseudocount_rpm)
                        rows.append({'feature_id': f['id'], 'category': f['category'],
                                     'numerator_count': a, 'denominator_count': b,
                                     'numerator_rpm': ar, 'denominator_rpm': br,
                                     'fold_change': ratio, 'log2_fold_change': math.log2(ratio),
                                     'direction': 'higher_in_numerator' if ar > br else 'higher_in_denominator' if ar < br else 'equal',
                                     'detection': 'both_samples' if a and b else 'numerator_only' if a else 'denominator_only'})
                    rows.sort(key=lambda r: (-r['log2_fold_change'], -max(r['numerator_rpm'], r['denominator_rpm']), r['feature_id']))
                    rankings[track] = [{**row, 'rank': rank} for rank, row in enumerate(rows, 1)]
                contrasts.append({'numerator_sample': numerator, 'denominator_sample': denominator,
                                  'rankings': rankings})
    return {'status': 'descriptive_no_replicate_statistics' if args.input_scope == 'full_library' else 'disabled_selected_subset',
            'input_scope': args.input_scope,
            'feature_table': annotations,
            'raw_count_matrices_path': 'count_matrices',
            'sample_metadata': [dict(metadata_by_id[s]) for s in samples],
            'normalized_rpm_matrices': normalized,
            'normalization': {'method': args.denominator, 'denominators': dict(denominators), 'pseudocount_rpm': args.pseudocount_rpm},
            'category_feature_count_matrices': category_counts,
            'pairwise_contrasts': contrasts,
            'pairwise_contrasts_deferred_to_frontend': defer_contrasts,
            'interpretation': {
                'features': 'Shared observed-read endpoint features across all annotated categories; molecular fragment ends remain unverified',
                'category_counts': 'Feature-eligible support only, not all mapped reads or total RNA abundance; inspect category_retention for exclusions',
                'rankings': 'Descending signed log2 RPM fold change, including control-only features; both-zero pairs omitted; no significance or automatic candidate selection',
                'future_statistics': 'Use existing integer count_matrices and aligned sample_metadata with genuine biological replicates; do not pass RPM or pseudocount-adjusted values as raw counts',
                'scope': 'Subset matrices describe supplied reads only; fold rankings disabled for selected subsets'}}


def analyze(args, tools, workdir):
    stage_started = time.perf_counter()
    timings = {}
    def checkpoint(name):
        nonlocal stage_started
        now = time.perf_counter()
        timings[name] = now - stage_started
        stage_started = now
        print(f'PERFORMANCE — {name}: {timings[name]:.3f} seconds', file=sys.stderr)
    workdir = str(Path(workdir).resolve())
    paths = {}
    for specification in args.sample:
        name, separator, path = specification.partition('=')
        if not separator or not re.fullmatch(r'[A-Za-z0-9_-]+', name) or not path or name in paths:
            raise ValueError('Use unique --sample NAME=PATH entries')
        paths[name] = Path(path).expanduser().resolve(strict=True)
    manifest = read_manifest(args.manifest) if args.manifest else None
    if manifest:
        paths = {row['sample_id']: Path(row['r1']) for row in manifest}
    exploratory = args.analysis_mode == 'exploratory'
    controls = [] if exploratory else args.control if args.control is not None else ['BMDM', 'BMDC']
    if not exploratory and (not controls or len(set(controls)) != len(controls) or args.target in controls):
        raise ValueError('Controls must be distinct from one another and the target')
    if not exploratory and not {args.target, *controls} <= set(paths):
        raise ValueError('Supply the target and every declared control sample')
    if len(set(paths.values())) != len(paths):
        raise ValueError('A file cannot represent more than one sample')
    reference_path = Path(args.reference).expanduser().resolve(strict=True)
    tracked = set(paths.values()) | {reference_path}
    if manifest:
        tracked |= {Path(row['r2']) for row in manifest} | {Path(args.manifest).expanduser().resolve()}
    if args.sample_metadata:
        tracked.add(Path(args.sample_metadata).expanduser().resolve(strict=True))
    for specification in args.diagnostic_reference:
        if '=' not in specification:
            raise ValueError('Diagnostic reference must be LABEL=FASTA')
        tracked.add(Path(specification.split('=', 1)[1]).expanduser().resolve(strict=True))
    output_path = Path(args.output).expanduser().resolve()
    if output_path in tracked:
        raise ValueError('Output must not overwrite any input/reference/metadata')
    if output_path.exists() and not args.overwrite:
        raise FileExistsError('Output exists; choose another path or explicitly use --overwrite')
    hashes = {path: sha256(path) for path in tracked}
    metadata = [{k: v for k, v in row.items() if k not in ('r1', 'r2')} for row in manifest] if manifest else sample_metadata(args.sample_metadata, paths)
    input_paths = dict(paths)
    paired_audit = {}
    if manifest:
        paths, paired_audit = prepare_manifest_inputs(manifest, args, workdir)
    qualities = {name: {} for name in paths}
    libraries, union_sequences = {}, set()
    for name, path in paths.items():
        libraries[name] = load_counts(path, args.counts_complete, args.input_scope == 'selected_subset',
                                     qualities[name] if args.quality_detail == 'positional' else None, args.max_distinct_sequences)
        union_sequences.update(libraries[name])
        if args.max_total_distinct_sequences and len(union_sequences) > args.max_total_distinct_sequences:
            raise ValueError(f'Combined distinct-sequence limit ({args.max_total_distinct_sequences:,}) exceeded after {name}; '
                             'no sequences discarded and no report published. Assess whole-run memory before increasing --max-total-distinct-sequences.')
    if any(sum(libraries[s].values()) != row['analyzed_r1_reads'] for s, row in paired_audit.items()):
        raise RuntimeError('Paired-input/R1 count conservation checkpoint failed')
    checkpoint('input_hashing_validation_and_counting')
    references, annotations = load_reference(reference_path, args.reference_id)
    model = reference_model(references, annotations, args)
    all_sequences = sorted(union_sequences)
    del union_sequences
    sequences = [seq for seq in all_sequences if len(seq) >= args.min_read_length]
    ids = {seq: index for index, seq in enumerate(sequences)}
    print(f'Mapping {len(sequences):,} distinct sequences across {len(libraries)} libraries', file=sys.stderr)
    checkpoint('reference_model_and_query_union')
    mappings = align_sequences(sequences, references, tools, args.threads, workdir)
    checkpoint('index_build_and_alignment')
    resolved = [resolve_mapping(mapping, model, annotations, len(seq), args.max_edit_rate)
                for seq, mapping in zip(sequences, mappings)]
    counts_by_id = {name: {ids[seq]: count for seq, count in counts.items() if seq in ids} for name, counts in libraries.items()}
    denominators, summaries, warnings, qc_alerts = {}, {}, [], []
    if manifest:
        warnings.append('PAIRED INPUT / R1 ONLY: both mates validated; only R1 aligned and counted. No insert reconstruction or native-end inference.')
    if exploratory:
        warnings.append('EXPLORATORY: no negative controls designated; candidate selection disabled. Samples are never pooled automatically.')
    if args.quality_detail == 'counts-only':
        warnings.append('Positional quality aggregation disabled; FASTQ quality syntax still validated. Missing quality evidence is not high-quality evidence.')
    if args.input_scope == 'selected_subset':
        warnings.append('SELECTED SUBSET: counts and QC describe supplied sequences only; biological fold comparisons/candidates disabled.')
    warnings.append('Read counts are not independent molecules without UMI/PCR information; extraction protocols may confound comparisons.')
    for sample, counts in libraries.items():
        total, mapped = sum(counts.values()), 0
        categories, states, assignments, lengths = Counter(), Counter(), Counter(), Counter()
        retention = defaultdict(Counter)
        for seq, count in counts.items():
            lengths[len(seq)] += count
            if seq not in ids:
                category = state = assignment = 'below_min_read_length'
            else:
                query = ids[seq]
                mapping, evidence = mappings[query], resolved[query]
                category, assignment = evidence['category'], evidence['status']
                mapped += count if mapping.hits else 0
                state = ('unmapped' if not mapping.hits else 'ambiguous_locus' if len(mapping.hits) > 1
                         else 'unique_exact' if next(iter(mapping.hits)).exact else 'unique_inexact')
                if mapping.hits:
                    retention[category]['mapped_reads'] += count
                    if len(mapping.hits) > 1:
                        retention[category]['raw_locus_ambiguous_reads'] += count
                    if evidence['feature_eligible']:
                        retention[category]['feature_eligible_reads'] += count
                        retention[category]['exact_feature_reads' if evidence['exact'] else 'inexact_feature_reads'] += count
                    else:
                        retention[category]['excluded_from_features_reads'] += count
            categories[category] += count
            states[state] += count
            assignments[assignment] += count
        denominator = total if args.denominator == 'total_input' else mapped
        if denominator <= 0:
            raise ValueError(f'Zero mapped denominator for {sample}; use total_input or inspect reference/input')
        denominators[sample] = denominator
        ambiguity_qc, alerts = category_ambiguity_qc(sample, retention, args.ambiguity_warning_fraction, args.input_scope)
        qc_alerts.extend(alerts)
        warnings.extend(alert['message'] for alert in alerts)
        summaries[sample] = {'path': str(input_paths[sample]), 'sha256': hashes[input_paths[sample]],
                             'input_scope': args.input_scope, 'total_input_reads': total,
                             'distinct_sequences': len(counts), 'mapped_reads': mapped, 'rpm_denominator': denominator,
                             'unmapped_reads': categories['unmapped'], 'unmapped_fraction_of_input': categories['unmapped'] / total,
                             'below_min_read_length_reads': categories['below_min_read_length'],
                             'category_counts': dict(sorted(categories.items())),
                             'category_rpm_per_input_million': {c: rpm(n, total) for c, n in sorted(categories.items())},
                             'mapping_state_counts': dict(sorted(states.items())), 'assignment_state_counts': dict(sorted(assignments.items())),
                             'category_retention': {c: dict(v) for c, v in sorted(retention.items())},
                             'category_ambiguity_qc': ambiguity_qc,
                             'length_counts': dict(sorted(lengths.items())),
                             'top10_sequences_all_categories': [{'sequence': seq.replace('T', 'U'), 'count': count,
                                 'rpm_per_input_million': rpm(count, total),
                                 'category': resolved[ids[seq]]['category'] if seq in ids else 'below_min_read_length'}
                                 for seq, count in sorted(counts.items(), key=lambda p: (-p[1], p[0]))[:10]]}
        if sum(categories.values()) != total or sum(states.values()) != total or sum(assignments.values()) != total:
            raise RuntimeError('Read accounting checkpoint failed')
        scope_label = 'SELECTED SUBSET ONLY' if args.input_scope == 'selected_subset' else 'supplied library'
        print(f'HEADLINE QC — {sample} [{scope_label}]: {categories["unmapped"]}/{total} primary-unmapped '
              f'({categories["unmapped"]/total:.1%} of input); {categories["below_min_read_length"]} below minimum length '
              f'(not assessed for mapping); {sum(v["feature_eligible_reads"] for v in retention.values())} eligible for fragment features', file=sys.stderr)
    for warning in warnings:
        print('WARNING: ' + warning, file=sys.stderr)
    if qc_alerts:
        print(f'\n!!! MAPPING AMBIGUITY QC: {len(qc_alerts)} sample/category warning(s). '
              'Review the warnings above and qc_alerts in the JSON report. !!!\n', file=sys.stderr)
    features, endpoints = build_features(resolved, sequences, counts_by_id, references, annotations, args.endpoint_tolerance)
    feature_ids = [f['id'] for f in features]
    if len(set(feature_ids)) != len(feature_ids):
        raise RuntimeError('Feature identity collision')
    feature_members = [q for f in features for q in f['member_ids']]
    if len(set(feature_members)) != len(feature_members):
        raise RuntimeError('A sequence was assigned to multiple fragment features')
    projected = [(q, Hit(parent, start, end, strand, '', 0, False)) for q, evidence in enumerate(resolved)
                 for parent, start, end, strand in evidence['locations']]
    overlap_index = IntervalIndex(projected)
    matrices = {track: {'feature_ids': feature_ids, 'sample_ids': list(paths),
                        'data': [[f['counts'][sample][track] for sample in paths] for f in features]}
                for track in ('exact', 'inexact', 'inclusive')}
    for sample in paths:
        eligible_total = sum(counts_by_id[sample].get(q, 0) for q, e in enumerate(resolved) if e['feature_eligible'])
        if sum(f['counts'][sample]['inclusive'] for f in features) != eligible_total:
            raise RuntimeError('Feature count conservation checkpoint failed')
    support_report = support_reporting(features, resolved, counts_by_id, summaries)
    for sample, support in support_report['samples'].items():
        n = support['counts']
        print(f'SUPPORT TRACKS — {sample} [{support["input_scope"]}]: exact-only={n["exact_only"]:,}; '
              f'accepted-inexact={n["accepted_inexact"]:,} '
              f'(raw-unique={n["raw_unique_inexact"]:,}, equivalence-resolved={n["equivalence_resolved_inexact"]:,}); '
              f'inclusive={n["inclusive"]:,}', file=sys.stderr)
    for feature in features:
        feature['rpm'] = {sample: {track: rpm(n, denominators[sample]) for track, n in counts.items()}
                          for sample, counts in feature['counts'].items()}
        feature['comparison'] = {track: compare_counts({s: feature['counts'][s][track] for s in paths}, denominators,
                                                      args.target, controls, args) for track in ('exact', 'inclusive')}
        feature['sequence_metrics'] = sequence_metrics(feature['observed_sequence'])
        if args.max_homopolymer is not None and feature['sequence_metrics']['longest_homopolymer'] > args.max_homopolymer:
            for comparison in feature['comparison'].values():
                if comparison['candidate'] is not None:
                    comparison['candidate'] = False
                    comparison['rejection_reasons'].append('representative_homopolymer_limit')
        feature['candidate'] = feature['comparison']['inclusive']['candidate']
        # Use the full observed boundary envelope for a conservative regional screen.
        start, end = feature['start_cluster']['minimum'], feature['end_cluster']['maximum']
        overlap_ids = overlap_index.overlap(feature['parent_RNA'], feature['strand'], start, end)
        if args.control_strand == 'both':
            overlap_ids |= overlap_index.overlap(feature['parent_RNA'], '-' if feature['strand'] == '+' else '+', start, end)
        feature['control_overlap_screen'] = {control: {'count': sum(counts_by_id[control].get(q, 0) for q in overlap_ids)} for control in controls}
        for control, values in feature['control_overlap_screen'].items():
            values['rpm'] = rpm(values['count'], denominators[control])
        feature['control_overlap_screen_is_abundance_comparison'] = False
    # Select folding/display budgets within each RNA category, after all features are counted.
    top_feature_ids, top_ids = {}, {}
    for sample in paths:
        by_category = defaultdict(list)
        for feature in features:
            if feature['counts'][sample]['inclusive']:
                by_category[feature['category']].append(feature)
        top_feature_ids[sample] = [f['id'] for category in sorted(by_category)
                                  for f in sorted(by_category[category], key=lambda f: (-f['counts'][sample]['inclusive'], f['id']))[:args.top_fragments]]
        mapped_categories = defaultdict(list)
        for q, evidence in enumerate(resolved):
            if mappings[q].hits and counts_by_id[sample].get(q, 0):
                mapped_categories[evidence['category']].append(q)
        top_ids[sample] = [q for category in sorted(mapped_categories)
                          for q in sorted(mapped_categories[category], key=lambda q: (-counts_by_id[sample][q], q))[:args.top_fragments]]
    selected_features = set().union(*(set(v) for v in top_feature_ids.values()))
    folded_queries = set().union(*(set(v) for v in top_ids.values()))
    cache = {}
    def fold(sequence):
        return fold_sequence(sequence, tools, args.temperature, args.max_fold_length, args.fold_timeout, workdir, cache)
    for feature in features:
        if feature['id'] in selected_features:
            ranked = sorted(feature['member_ids'], key=lambda q: (-sum(c.get(q, 0) for c in counts_by_id.values()), q))
            chosen = ranked[:args.max_fold_variants]
            folded_queries.update(chosen)
            feature['observed_variant_folds'] = [{'read_id': q, 'sequence': sequences[q].replace('T', 'U'), 'fold': fold(sequences[q].replace('T', 'U'))} for q in chosen]
            feature['reference_fold'] = fold(feature['reference_sequence'])
            feature['unfolded_variant_count'] = len(ranked) - len(chosen)
        else:
            feature['observed_variant_folds'] = []
            feature['reference_fold'] = {'status': 'outside_folding_budget', 'structure': None, 'mfe_kcal_mol': None}
            feature['unfolded_variant_count'] = len(feature['member_ids'])
    checkpoint('assignment_features_and_folding')
    diagnostics, diagnostic_sources = unresolved_diagnostics(sequences, mappings, args, tools, workdir)
    checkpoint('unmapped_diagnostics')
    unmapped_reports, unmapped_labels = unmapped_qc(
        diagnostics, counts_by_id, summaries, diagnostic_categories(args), bool(args.adapter_sequence))
    for sample, qc in unmapped_reports.items():
        print(f'UNMAPPED INVESTIGATION — {sample}: ' + '; '.join(
            f'{label}={screen["status"]}' + (f' ({screen["matching_reads"]} reads)' if screen['matching_reads'] is not None else '')
            for label, screen in qc['screens'].items()) +
            '. Matches are non-exclusive evidence, not confirmed origin.', file=sys.stderr)
    fragments = []
    # Every eligible and unresolved sequence is retained; no majority is reduced to a counter.
    for q, sequence in (enumerate(sequences) if args.report_detail == 'complete' else []):
        evidence = resolved[q]
        hit = next(iter(mappings[q].hits)) if len(mappings[q].hits) == 1 else None
        quality = {sample: qualities[sample].get(sequence) for sample in paths}
        row = {'id': q, 'sequence': sequence.replace('T', 'U'), 'assignment': evidence,
               'counts': {s: counts_by_id[s].get(q, 0) for s in paths},
               'rpm': {s: rpm(counts_by_id[s].get(q, 0), denominators[s]) for s in paths},
               'quality_evidence_phred33': quality, 'best_alignment_score': mappings[q].score,
               'best_alignments': [{**h.__dict__, 'differences': alignment_differences(sequence, h, references)} for h in sorted(mappings[q].hits)],
               'parent_rRNA': hit.reference if hit else None, 'start': hit.start if hit else None,
               'parent_RNA': hit.reference if hit else None,
               'end': hit.end if hit else None, 'strand': hit.strand if hit else None,
               'mapping_status': ('unmapped' if not mappings[q].hits else 'ambiguous_locus' if not hit else 'unique_exact' if hit.exact else 'unique_inexact'),
               'assembly_eligible': evidence['feature_eligible'], 'sequence_metrics': sequence_metrics(sequence),
               'diagnostic_evidence': diagnostics.get(q, []),
               'unmapped_diagnostic_labels': unmapped_labels.get(q),
               'fold': fold(sequence.replace('T', 'U')) if q in folded_queries else {'status': 'outside_folding_budget', 'structure': None, 'mfe_kcal_mol': None}}
        fragments.append(row)
    short_sequences = [{'sequence': seq.replace('T', 'U'), 'counts': {s: c.get(seq, 0) for s, c in libraries.items()},
                        'reason': 'below_min_read_length'} for seq in all_sequences if seq not in ids and args.report_detail == 'complete']
    regions = []
    for sample in paths:
        entries = [(q, Hit(*e['locations'][0], '', 0, e['exact'])) for q, e in enumerate(resolved)
                   if e['feature_eligible'] and counts_by_id[sample].get(q, 0)]
        for index, region in enumerate(merge_intervals(entries), 1):
            region.update({'id': f'{sample}_coverage_{index}', 'sample': sample, 'candidate': False,
                           'role': 'coverage_summary_not_a_molecule', 'parent_RNA': region['parent_rRNA'],
                           'fragments_merged': len(region['member_ids']), 'contig_length': region['end'] - region['start']})
            region['seed_count'] = sum(counts_by_id[sample][q] for q in region['member_ids'])
            region['seed_rpm'] = rpm(region['seed_count'], denominators[sample])
            region['consensus_sequence'] = sequence_for_region(region, references)
            region['fold'] = fold(region['consensus_sequence']) if args.fold_coverage_regions else {'status': 'coverage_region_not_folded', 'structure': None, 'mfe_kcal_mol': None}
            regions.append(region)
    family_ids = sorted({e['family'] for e in resolved if e['family'] is not None})
    family_counts = {family: {s: 0 for s in paths} for family in family_ids}
    for q, evidence in enumerate(resolved):
        if evidence['family']:
            for sample in paths:
                family_counts[evidence['family']][sample] += counts_by_id[sample].get(q, 0)
    ranking = sorted((f for f in features if f['comparison']['inclusive'].get('target_count', 0)),
                     key=lambda f: (-f['comparison']['inclusive']['minimum_control_fold_enrichment'],
                                    -f['comparison']['inclusive']['target_rpm'], f['id']))
    if any(sha256(path) != checksum for path, checksum in hashes.items()):
        raise RuntimeError('Input changed during analysis; report not published')
    checkpoint('read_reporting_coverage_and_integrity_recheck')
    return {'schema_version': VERSION, 'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'complete',
            'performance': {'stage_seconds': timings, 'mapped_query_count': len(sequences),
                            'scope': 'Elapsed stage times through integrity recheck; excludes final report assembly/serialization. Not a full-library capacity estimate.'},
            'target': None if exploratory else args.target, 'controls': controls, 'parameters': vars(args).copy(), 'tools': tools,
            'paired_input_qc': paired_audit,
            'report_omissions': (['individual_read_alignments_and_quality', 'short_sequence_records', 'unused_reference_annotations', 'whole_reference_alias_and_family_maps']
                                 if args.report_detail == 'summary' else []),
            'reference': {'path': str(reference_path), 'sha256': hashes[reference_path],
                          'records': {n: {**annotations[n], 'length': len(references[n])} for n in
                                      (references if args.report_detail == 'complete' else sorted({f['parent_RNA'] for f in features}))},
                          'canonical_aliases': model['aliases'] if args.report_detail == 'complete' else {},
                          'verified_projections': model['projections'],
                          'families': model['families'] if args.report_detail == 'complete' else {}},
            'input_provenance': {str(path): checksum for path, checksum in sorted(hashes.items())},
            'sample_metadata': metadata, 'samples': summaries, 'warnings': warnings, 'qc_alerts': qc_alerts,
            'features': features, 'count_matrices': matrices, 'endpoint_profiles': endpoints,
            'support_reporting': support_report,
            'landscape_comparison': landscape_comparison(features, matrices, metadata, denominators, args),
            'parent_landscape': parent_landscape(features, references, paths, args, resolved, counts_by_id),
            'family_count_matrix': {'family_ids': family_ids, 'sample_ids': list(paths), 'data': [[family_counts[f][s] for s in paths] for f in family_ids]},
            'top_feature_ids': top_feature_ids, 'top_fragment_ids': top_ids, 'fragments': fragments,
            'unresolved_read_ids': [q for q, e in enumerate(resolved) if not e['feature_eligible']],
            'unmapped_read_ids': [q for q, m in enumerate(mappings) if not m.hits],
            'unmapped_qc': unmapped_reports, 'diagnostic_reference_categories': diagnostic_categories(args),
            'excluded_short_sequences': short_sequences, 'diagnostic_reference_provenance': diagnostic_sources,
            'regions': regions, 'candidate_region_ids': [],
            'candidate_feature_ids': [f['id'] for f in features if f['candidate']],
            'fold_change_ranking': [f['id'] for f in ranking],
            'interpretation': {
                'coordinates': '0-based half-open; strand is read-to-reference orientation, not automatically biological strand',
                'features': 'Shared paired-endpoint clusters supported by spanning reads; sequence variants retained; read-length censoring can hide molecular ends',
                'tracks': 'Exact and accepted-inexact counts are disjoint; inclusive is their sum. Eligibility is symmetric across samples.',
                'family_counts': 'Any resolved-family best hit counted once, even when endpoint/quality eligibility fails; not additive with feature counts',
                'quality': ('Phred+33 per-position evidence retained' if args.quality_detail == 'positional' else 'Positional quality evidence omitted by request') + '; deduplicated FASTA alignment uses quality-independent scoring',
                'equivalence': 'Identical-sequence/verified coordinate equivalence resolves RNA sequence evidence, not genomic copy of origin',
                'diagnostics': 'Optional adapter or secondary-reference matches are non-exclusive clues; not confirmed taxonomy or contamination',
                'coverage': 'Merged intervals are secondary summaries and never asserted to be intact molecules',
                'folding': 'Canonical-base observed variants and reference separately; no modeled glycosylation or proven accessibility',
                'statistics': 'Descriptive fold rankings only; no p-values, significance, or replicate inference',
                'limitations': ['Unmapped identity remains unresolved without suitable diagnostic references and experimental context',
                                'Genome diagnostic alignment is unspliced; failure is not proof of non-genomic origin',
                                'Bowtie2 best hits are heuristic search results; mismatches do not establish RNA modifications',
                                'Duplicate read counts do not prove independent endpoints; protocol/UMI information is required',
                                'Total-input and mapped RPM answer different conditional questions; neither fixes composition/extraction confounding']}}


def atomic_report(report, output, overwrite=False):
    output = Path(output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                         prefix=".csrna-", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(report, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary, output)
        else:
            # Atomic no-clobber publish, including concurrent writers.
            os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--sample", action="append", default=[], metavar="NAME=PATH")
    p.add_argument('--manifest', help='CSV with sample_id,r1,r2,condition,replicate; relative FASTQ paths resolve beside CSV')
    p.add_argument('--pair-handling', choices=['r1-only'], help='Explicit paired-manifest policy: validate both mates, analyze R1 only')
    p.add_argument('--analysis-mode', choices=['target_control', 'exploratory'], default='target_control',
                   help='Exploratory mode needs no target/controls and never calls control-dependent candidates')
    p.add_argument('--pilot-pairs', type=int, default=0, help='Uniform pairs sampled per manifest library; 0 uses all; requires selected_subset')
    p.add_argument('--sampling-seed', type=int, default=20261007)
    p.add_argument('--report-detail', choices=['complete', 'summary'], default='complete', help='Summary omits individual-read alignments and unused reference records; features/counts/QC are unchanged')
    p.add_argument('--quality-detail', choices=['positional', 'counts-only'], default='positional', help='Counts-only validates FASTQ but avoids per-base quality arrays; no quality evidence inferred')
    p.add_argument("--target", default="M1")
    p.add_argument("--control", action="append", help="Repeat for controls; default BMDM and BMDC")
    p.add_argument("--reference", help="Combined or curated RNA FASTA, optionally gzipped")
    p.add_argument("--reference-id", action="append", default=[], help="Optional exact FASTA IDs to include; repeat")
    p.add_argument('--reference-projection', action='append', default=[], help='SOURCE=TARGET:OFFSET:+/-; validated exact sequence projection')
    p.add_argument('--no-auto-rrna-projections', action='store_true', help='Disable annotated precursor-to-mature exact containment projections')
    p.add_argument('--parent-family', action='append', default=[], help='REFERENCE_ID=FAMILY; family membership does not assert coordinate equivalence')
    p.add_argument('--input-scope', choices=['full_library', 'selected_subset'], default='full_library', help='Selected subsets disable biological comparisons/candidates')
    p.add_argument('--sample-metadata', help='CSV with sample_id and optional condition, replicate, library_preparation')
    p.add_argument('--diagnostic-reference', action='append', default=[], help='LABEL=FASTA; optional secondary mapping of primary-unmapped reads')
    p.add_argument('--diagnostic-reference-category', action='append', default=[],
                   help='LABEL=genome|microbial|other; explicit diagnostic source category (default: other); repeat for each labeled reference')
    p.add_argument('--adapter-sequence', action='append', default=[], help='Exact adapter substring diagnostic; repeat, minimum 8 bases')
    p.add_argument('--endpoint-tolerance', type=int, default=2, help='Maximum nt distance to endpoint cluster center; diameter <=2*tolerance')
    p.add_argument('--max-edit-rate', type=float, default=0.10, help='NM/query length limit for endpoint features, after best-hit assignment')
    p.add_argument('--ambiguity-warning-fraction', type=float, default=0.5,
                   help='Warn when at least this fraction of mapped reads in a category has multiple raw best-score loci (default: 0.5); QC only')
    p.add_argument('--min-landscape-reads', type=int, default=10)
    p.add_argument('--min-species-reads', type=int, default=3, help='Read-count support for recurrent endpoint patterns, not independent molecules')
    p.add_argument('--landscape-dominance', type=float, default=0.5)
    p.add_argument('--max-fold-variants', type=int, default=5)
    p.add_argument('--fold-coverage-regions', action='store_true', help='Optional reference-region folds, explicitly not observed molecules')
    p.add_argument("--output", default="csrna_results.json")
    p.add_argument("--counts-complete", action="store_true", help="Attest that CSV inputs contain full-library sequence counts")
    p.add_argument('--max-distinct-sequences', type=int, default=250000,
                   help='Per-library distinct-sequence ingestion limit; fail explicitly rather than discard reads (0 disables; memory is not bounded during mapping)')
    p.add_argument('--max-total-distinct-sequences', type=int, default=500000,
                   help='Whole-run distinct-sequence union limit checked after each library; 0 disables. This is not a RAM guarantee.')
    p.add_argument("--denominator", choices=["total_input", "mapped"], default="total_input")
    p.add_argument("--min-rpm", type=float, default=100.0)
    p.add_argument("--min-enrichment", type=float, default=4.0)
    p.add_argument("--pseudocount-rpm", type=float, default=1.0)
    p.add_argument("--require-control-absence", action="store_true")
    p.add_argument("--control-strand", choices=["both", "same"], default="both",
                   help="Both is conservative across library protocols; same requires compatible stranded libraries")
    p.add_argument("--min-read-length", type=int, default=15)
    p.add_argument("--temperature", type=float, default=37.0)
    p.add_argument("--max-fold-length", type=int, default=2000)
    p.add_argument("--fold-timeout", type=float, default=120.0)
    p.add_argument("--top-fragments", type=int, default=100, help='Folding/report ranking budget per category per sample; never a count filter')
    p.add_argument("--max-homopolymer", type=int, help="Optional explicit screen; no universal biological cutoff assumed")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--temp-dir", help="Scratch parent; created run directory is always cleaned on normal/error exit")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--self-test", action="store_true", help="Run built-in deterministic checks without external tools")
    p.add_argument("--self-test-tools", action="store_true", help="Also run a tiny full workflow using real Bowtie2/RNAfold")
    return p


def validate_args(args):
    diagnostic_categories(args)
    if args.manifest and (args.sample or args.sample_metadata):
        raise ValueError('Use either --manifest or --sample/--sample-metadata, not both')
    if args.manifest and args.pair_handling != 'r1-only':
        raise ValueError('Manifest requires explicit --pair-handling r1-only; native paired alignment is not implemented')
    if args.pair_handling and not args.manifest:
        raise ValueError('--pair-handling requires --manifest')
    if args.analysis_mode == 'exploratory' and (args.control is not None or args.require_control_absence):
        raise ValueError('Exploratory mode cannot specify controls or control-dependent screens')
    if args.pilot_pairs < 0 or (args.pilot_pairs and (not args.manifest or args.input_scope != 'selected_subset')):
        raise ValueError('--pilot-pairs must be nonnegative and needs --manifest with --input-scope selected_subset')
    if args.max_distinct_sequences < 0 or args.max_total_distinct_sequences < 0:
        raise ValueError('Distinct-sequence limits must be nonnegative')
    if args.endpoint_tolerance < 0:
        raise ValueError('endpoint_tolerance must be nonnegative')
    for name in ('max_edit_rate', 'ambiguity_warning_fraction', 'landscape_dominance'):
        value = getattr(args, name)
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f'{name} must lie in [0,1]')
    if args.landscape_dominance <= 0:
        raise ValueError('landscape_dominance must be positive')
    for name in ("min_rpm", "min_enrichment", "pseudocount_rpm", "fold_timeout"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(args.temperature) or not 0 <= args.temperature <= 100:
        raise ValueError("temperature must be between 0 and 100 Celsius")
    for name in ("threads", "min_read_length", "max_fold_length", "top_fragments", "min_landscape_reads", "min_species_reads", "max_fold_variants"):
        if getattr(args, name) <= 0:
            raise ValueError(f"{name} must be positive")
    if args.min_read_length < 4:
        raise ValueError("min_read_length must be at least 4")
    if args.max_homopolymer is not None and args.max_homopolymer <= 0:
        raise ValueError("max_homopolymer must be positive")
    if not args.reference or not (args.sample or args.manifest):
        raise ValueError('--reference and either --sample or --manifest are required')


# %% Built-in regression tests (no production inputs are modified)
class PipelineTests(unittest.TestCase):
    def test_argument_validation(self):
        for option, value in (("--min-rpm", "nan"), ("--min-enrichment", "inf"),
                              ("--pseudocount-rpm", "0"), ("--min-read-length", "1"),
                              ("--temperature", "101"), ("--max-homopolymer", "0")):
            args = parser().parse_args(["--reference", "r.fa", "--sample", "M1=m.fa", option, value])
            with self.assertRaises(ValueError):
                validate_args(args)

    def test_alphabet(self):
        self.assertEqual(clean_sequence(" acgu "), "ACGT")
        self.assertEqual(reverse_complement("AACG"), "CGTT")
        for bad in ("", "AC-G", "AC GT", "AX"):
            with self.assertRaises(ValueError):
                clean_sequence(bad)

    def test_counts_and_fastq(self):
        import io
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.csv"
            path.write_text("sequence,count\nacgu,2\nACGT,3\nTTTT,0\n")
            self.assertEqual(load_counts(path, True), Counter({"ACGT": 5}))
            with self.assertRaises(ValueError):
                load_counts(path)
            path.write_text("sequence,count\nACGT,-1\n")
            with self.assertRaises(ValueError):
                load_counts(path, True)
            path = Path(directory) / "reads.fa.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(">one\nAC\nGT\n>two\nACGU\n")
            self.assertEqual(load_counts(path), Counter({"ACGT": 2}))
        self.assertEqual(list(fastq_records(io.StringIO("@x\nACGT\n+\nIIII\n"))), ["ACGT"])
        for text in ("@x\nACGT\n+\n", "@x\nACGT\n+\nIII\n", "@x\nACGT\n+y\nIIII\n"):
            with self.assertRaises(ValueError):
                list(fastq_records(io.StringIO(text)))

    def test_reference_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r.fa"
            path.write_text(">RIBO|a rRNA\nACGT\n>g protein_coding\nAACG\n")
            refs, annotations = load_reference(path)
            self.assertEqual(annotations["g"]["category"], "mRNA")
            self.assertEqual(load_reference(path, ["RIBO|a"])[0], {"RIBO|a": "ACGT"})
            with self.assertRaises(ValueError):
                load_reference(path, ["missing"])
            path.write_text(">RIBO|a\nACGT\n>RIBO|a\nAAAA\n")
            with self.assertRaises(ValueError):
                load_reference(path)

    def test_reject_legacy_and_empty_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.csv"
            for text in ("sequence,count,RPM_total_mapped\nACGT,1,100\n",
                         "sequence,count\nACGT,0\n", "sequence,count\nACGT,1.5\n",
                         "sequence,count,counts\nACGT,1,1\n"):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    load_counts(path, True)

    def test_sam_exact_and_balanced_indels(self):
        refs = {"r": "AACCGGTT"}
        def sam(cigar, nm=0, flag=0, query="AACCGGTT"):
            return f"q0\t{flag}\tr\t1\t42\t{cigar}\t*\t0\t0\t{query}\t*\tAS:i:0\tNM:i:{nm}"
        self.assertTrue(parse_sam_record(sam("8M"), ["AACCGGTT"], refs)[2].exact)
        self.assertFalse(parse_sam_record(sam("2M1D3M1I2M", 2), ["AACCGGTT"], refs)[2].exact)
        self.assertFalse(parse_sam_record(sam("8M", 1), ["AACCGGTA"], refs)[2].exact)
        refs = {"r": "AACCGGTA"}
        self.assertTrue(parse_sam_record(sam("8M", flag=16), ["TACCGGTT"], refs)[2].exact)

    def test_ambiguity_and_best_score(self):
        hit = Hit("r", 0, 8, "+", "8M", 0, True)
        other = Hit("r", 10, 18, "+", "8M", 0, True)
        mapping = Mapping()
        mapping.add(-2, other)
        mapping.add(0, hit)
        self.assertEqual(mapping.hits, {hit})
        mapping.add(0, other)
        self.assertIsNone(unique_exact_rrna(mapping, {"r": {"category": "rRNA"}}))

    def test_cross_category_ambiguity(self):
        mapping = Mapping(0, {Hit("r", 0, 8, "+", "8M", 0, True), Hit("g", 0, 8, "+", "8M", 0, True)})
        annotations = {"r": {"category": "rRNA"}, "g": {"category": "mRNA"}}
        self.assertEqual(mapping_category(mapping, annotations), "multi_category_ambiguous")
        self.assertIsNone(unique_exact_rrna(mapping, annotations))

    def test_interval_index_against_bruteforce(self):
        import random
        rng = random.Random(3)
        entries = []
        for query in range(80):
            start = rng.randrange(100)
            entries.append((query // 2, Hit("r", start, start + rng.randrange(1, 30), "+", "1M", 0, True)))
        index = IntervalIndex(entries)
        for _ in range(100):
            start = rng.randrange(100)
            end = start + rng.randrange(1, 30)
            expected = {query for query, hit in entries if hit.start < end and hit.end > start}
            self.assertEqual(index.overlap("r", "+", start, end), expected)

    def test_intervals_strand_adjacency_nesting(self):
        def hit(start, end, strand="+", ref="r"):
            return Hit(ref, start, end, strand, f"{end-start}M", 0, True)
        entries = [(0, hit(0, 10)), (1, hit(2, 4)), (2, hit(8, 15)),
                   (3, hit(15, 20)), (4, hit(0, 10, "-")), (5, hit(0, 10, ref="s"))]
        regions = merge_intervals(entries)
        self.assertEqual(len(regions), 4)
        self.assertEqual(regions[0]["member_ids"], [0, 1, 2])
        index = IntervalIndex(entries + [(0, hit(4, 12))])
        self.assertEqual(index.overlap("r", "+", 4, 8), {0})
        self.assertEqual(index.overlap("r", "+", 15, 16), {3})
        self.assertEqual(sequence_for_region({"parent_rRNA": "r", "start": 0, "end": 4, "strand": "-"}, {"r": "AACG"}), "CGUU")

    def test_fold_parsing(self):
        self.assertEqual(parse_fold("AAAA\n.... (  0.00)\n", "AAAA")["structure"], "....")
        self.assertEqual(parse_fold("GCGC\n(()) (-2.30)\n", "GCGC")["mfe_kcal_mol"], -2.3)
        for result in ("AAAA\n... (0.0)", "AAAA\n)(.. (0.0)", "CCCC\n.... (0.0)", "error"):
            with self.assertRaises(ValueError):
                parse_fold(result, "AAAA")

    def test_fold_failure_and_length_limit(self):
        from unittest.mock import patch
        tools = {"RNAfold": {"path": "RNAfold"}}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(subprocess, "run") as run:
                skipped = fold_sequence("AAAA", tools, 37, 3, 10, directory, {})
                self.assertIsNone(skipped["mfe_kcal_mol"])
                run.assert_not_called()
                run.return_value = subprocess.CompletedProcess([], 1, "", "synthetic failure")
                with self.assertRaises(RuntimeError):
                    fold_sequence("AAAA", tools, 37, 100, 10, directory, {})

    def test_region_controls_count_once(self):
        args = parser().parse_args([])
        hits = [(0, Hit("r", 0, 10, "+", "10M", 0, True)),
                (1, Hit("r", 5, 15, "+", "10M", 1, False)),
                (1, Hit("r", 7, 17, "+", "10M", 1, False))]
        region = {"parent_rRNA": "r", "strand": "+", "start": 0, "end": 10}
        result = compare_region(region, "M1", ["BMDM", "BMDC"],
                                {"M1": {0: 100}, "BMDM": {1: 50}, "BMDC": {}},
                                dict.fromkeys(["M1", "BMDM", "BMDC"], 1000),
                                IntervalIndex(hits[:1]), IntervalIndex(hits), args)
        self.assertEqual(result["controls"]["BMDM"]["overlap_count"], 50)
        self.assertFalse(result["abundance_control_pass"])
        self.assertTrue(result["controls"]["BMDC"]["observed_absent"])

    def test_opposite_strand_control_policy(self):
        args = parser().parse_args([])
        target_hit = Hit("r", 0, 10, "+", "10M", 0, True)
        control_hit = Hit("r", 0, 10, "-", "10M", 0, True)
        region = {"parent_rRNA": "r", "strand": "+", "start": 0, "end": 10}
        counts = {"M1": {0: 100}, "C": {1: 100}}
        denominators = {"M1": 100, "C": 100}
        exact, all_hits = IntervalIndex([(0, target_hit)]), IntervalIndex([(0, target_hit), (1, control_hit)])
        self.assertFalse(compare_region(region, "M1", ["C"], counts, denominators, exact, all_hits, args)["abundance_control_pass"])
        args.control_strand = "same"
        self.assertTrue(compare_region(region, "M1", ["C"], counts, denominators, exact, all_hits, args)["abundance_control_pass"])

    def test_absence_is_separate_from_enrichment(self):
        args = parser().parse_args([])
        hit = Hit("r", 0, 10, "+", "10M", 0, True)
        index = IntervalIndex([(0, hit)])
        region = {"parent_rRNA": "r", "strand": "+", "start": 0, "end": 10}
        counts, denominators = {"M1": {0: 100}, "C": {0: 1}}, {"M1": 1000, "C": 1000}
        self.assertTrue(compare_region(region, "M1", ["C"], counts, denominators, index, index, args)["abundance_control_pass"])
        args.require_control_absence = True
        result = compare_region(region, "M1", ["C"], counts, denominators, index, index, args)
        self.assertFalse(result["abundance_control_pass"])
        self.assertIn("observed_in_control:C", result["rejection_reasons"])

    def test_atomic_no_clobber_and_nan(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            atomic_report({"ok": 1}, path)
            with self.assertRaises(FileExistsError):
                atomic_report({"ok": 2}, path)
            with self.assertRaises(ValueError):
                atomic_report({"bad": float("nan")}, path, True)
            self.assertEqual(json.loads(path.read_text()), {"ok": 1})
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_empty_result_and_normalization(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            reference = directory / "r.fa"
            reference.write_text(">RIBO|r\n" + "A" * 60 + "\n")
            cli = ["--reference", str(reference), "--output", str(directory / "out.json"), "--counts-complete"]
            for name in ("M1", "BMDM", "BMDC"):
                path = directory / f"{name}.csv"
                path.write_text("sequence,count\n" + "A" * 30 + ",7\nACGT,3\n")
                cli.extend(["--sample", f"{name}={path}"])
            args = parser().parse_args(cli)
            with patch.dict(analyze.__globals__, {"align_sequences": lambda *a: [Mapping()]}):
                report = analyze(args, {}, directory)
                self.assertEqual(report["regions"], [])
                self.assertEqual(report["candidate_region_ids"], [])
                self.assertEqual(report["samples"]["M1"]["rpm_denominator"], 10)
                self.assertEqual(report["samples"]["M1"]["category_counts"], {"below_min_read_length": 3, "unmapped": 7})
                args.denominator = "mapped"
                with self.assertRaises(ValueError):
                    analyze(args, {}, directory)


# %% Fragment redesign regression tests
class RedesignTests(unittest.TestCase):
    @staticmethod
    def landscape_fixture():
        samples = ['M1', 'BMDM', 'BMDC']
        features = []
        for identifier, category, counts in [('r', 'rRNA', [(10, 2), (20, 4), (0, 0)]),
                                              ('s', 'snoRNA', [(0, 0), (0, 8), (0, 2)]),
                                              ('l', 'lncRNA', [(5, 0), (0, 0), (0, 0)])]:
            features.append({'id': identifier, 'parent_RNA': identifier, 'category': category,
                             'strand': '+', 'start': 0, 'end': 30, 'start_cluster': {}, 'end_cluster': {},
                             'endpoint_evidence': 'sequenced_read_boundaries',
                             'counts': {s: {'exact': a, 'inexact': b, 'inclusive': a + b}
                                        for s, (a, b) in zip(samples, counts)}})
        matrices = {t: {'feature_ids': [f['id'] for f in features], 'sample_ids': samples,
                        'data': [[f['counts'][s][t] for s in samples] for f in features]}
                    for t in ('exact', 'inexact', 'inclusive')}
        metadata = [{'sample_id': s, 'condition': s, 'replicate': None} for s in reversed(samples)]
        return features, matrices, metadata, {'M1': 100, 'BMDM': 200, 'BMDC': 100}, parser().parse_args([])

    def test_landscape_union_normalization_and_rankings(self):
        inputs = self.landscape_fixture()
        before = json.dumps(inputs[:4], sort_keys=True)
        result = landscape_comparison(*inputs)
        self.assertEqual([r['sample_id'] for r in result['sample_metadata']], ['M1', 'BMDM', 'BMDC'])
        self.assertEqual(len(result['pairwise_contrasts']), 3)
        rows = result['pairwise_contrasts'][0]['rankings']['inclusive']
        self.assertEqual([r['feature_id'] for r in rows], ['l', 'r', 's'])
        self.assertEqual(rows[1]['fold_change'], 1)
        self.assertEqual(rows[2]['detection'], 'denominator_only')
        self.assertEqual(rows[2]['direction'], 'higher_in_denominator')
        self.assertAlmostEqual(rows[2]['fold_change'], 1 / 40001)
        self.assertEqual(result['normalized_rpm_matrices']['inclusive']['data'][0], [120000, 120000, 0])
        self.assertEqual(result['category_feature_count_matrices']['inclusive']['data'], [[5, 0, 0], [12, 24, 0], [0, 8, 2]])
        self.assertEqual(before, json.dumps(inputs[:4], sort_keys=True))

    def test_landscape_subset_retains_tables_without_rankings(self):
        inputs = self.landscape_fixture()
        inputs[-1].input_scope = 'selected_subset'
        result = landscape_comparison(*inputs)
        self.assertEqual(result['status'], 'disabled_selected_subset')
        self.assertFalse(result['pairwise_contrasts'])
        self.assertEqual(len(result['feature_table']), 3)

    def test_landscape_empty_and_both_zero_tracks(self):
        f, m, meta, den, args = self.landscape_fixture()
        result = landscape_comparison(f, m, meta, den, args)
        self.assertEqual([r['feature_id'] for r in result['pairwise_contrasts'][0]['rankings']['exact']], ['l', 'r'])
        for matrix in m.values():
            matrix['feature_ids'], matrix['data'] = [], []
        result = landscape_comparison([], m, meta, den, args)
        self.assertFalse(result['feature_table'])
        self.assertTrue(all(not rows for pair in result['pairwise_contrasts'] for rows in pair['rankings'].values()))

    def test_landscape_rejects_misaligned_metadata_and_counts(self):
        f, m, meta, den, args = self.landscape_fixture()
        with self.assertRaisesRegex(RuntimeError, 'metadata'):
            landscape_comparison(f, m, meta[:2], den, args)
        m['exact']['data'][0][0] += 1
        with self.assertRaisesRegex(RuntimeError, 'disagrees'):
            landscape_comparison(f, m, meta, den, args)

    def test_distinct_limit_never_truncates_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'reads.fq'
            path.write_text('@a\nACGT\n+\nIIII\n@b\nACGT\n+\n!!!!\n@c\nAAAA\n+\nIIII\n')
            with self.assertRaisesRegex(ValueError, 'Distinct-sequence limit'):
                load_counts(path, max_distinct_sequences=1)
            self.assertEqual(load_counts(path, max_distinct_sequences=2), {'ACGT': 2, 'AAAA': 1})
            self.assertEqual(load_counts(path, max_distinct_sequences=0), {'ACGT': 2, 'AAAA': 1})
            for suffix, contents in (('.fa', '>a\nACGT\n>b\nAAAA\n'),
                                     ('.csv', 'sequence,count\nACGT,2\nAAAA,1\n')):
                other = Path(directory) / ('reads' + suffix)
                other.write_text(contents)
                with self.assertRaisesRegex(ValueError, 'Distinct-sequence limit'):
                    load_counts(other, counts_complete=True, max_distinct_sequences=1)

    def test_landscape_full_length_uses_sample_reads(self):
        args = parser().parse_args([])
        feature = {'parent_RNA': 'p', 'strand': '+', 'start': 0, 'end': 100,
                   'member_ids': [0, 1], 'counts': {'A': {'inclusive': 30}, 'B': {'inclusive': 20}}}
        resolved = [{'locations': [('p', 0, 100, '+')]}, {'locations': [('p', 4, 96, '+')]}]
        rows = parent_landscape([feature], {'p': 'A' * 100}, ['A', 'B'], args,
                                resolved, {'A': {0: 30}, 'B': {1: 20}})
        self.assertEqual(rows[0]['full_length_compatible_fraction'], 1)
        self.assertEqual(rows[1]['full_length_compatible_fraction'], 0)
        self.assertNotEqual(rows[1]['descriptive_pattern'], 'full_length_compatible')

    def test_support_report_disjoint_tracks_and_normalization(self):
        resolved = [{'feature_eligible': True, 'status': status, 'exact': exact}
                    for status, exact in [('unique_locus', True), ('unique_locus', False),
                                          ('equivalent_parent_resolved', True), ('equivalent_parent_resolved', False)]]
        features = [{'id': 'f', 'category': 'rRNA', 'member_ids': [0, 1, 2, 3],
                     'counts': {'S': {'exact': 10, 'inexact': 20, 'inclusive': 30}}}]
        summaries = {'S': {'rpm_denominator': 100, 'input_scope': 'selected_subset', 'category_retention': {'rRNA': {}}}}
        before = json.dumps(features, sort_keys=True)
        result = support_reporting(features, resolved, {'S': {0: 3, 1: 7, 2: 7, 3: 13}}, summaries)
        row = result['samples']['S']
        self.assertEqual(row['counts']['raw_unique_inexact'], 7)
        self.assertEqual(row['counts']['equivalence_resolved_inexact'], 13)
        self.assertEqual(row['rpm']['accepted_inexact'], 200000)
        self.assertEqual(row['inclusive_to_exact_support_ratio'], 3)
        self.assertAlmostEqual(row['inexact_fraction_of_eligible_support'], 2 / 3)
        self.assertEqual(row['input_scope'], 'selected_subset')
        self.assertEqual(result['assignment_count_matrices']['raw_unique_exact']['data'], [[3]])
        self.assertEqual(before, json.dumps(features, sort_keys=True))

    def test_support_report_zero_support(self):
        summaries = {'S': {'rpm_denominator': 100, 'input_scope': 'complete', 'category_retention': {'rRNA': {}}}}
        result = support_reporting([], [], {'S': {}}, summaries)
        row = result['samples']['S']
        self.assertEqual(row['counts']['inclusive'], 0)
        self.assertIsNone(row['inclusive_to_exact_support_ratio'])
        self.assertIsNone(row['inexact_fraction_of_eligible_support'])
        self.assertEqual(row['by_category']['rRNA']['counts']['accepted_inexact'], 0)
        self.assertEqual(result['assignment_count_matrices']['raw_unique_inexact']['data'], [])

    def test_support_report_rejects_bad_accounting(self):
        summaries = {'S': {'rpm_denominator': 10, 'input_scope': 'complete', 'category_retention': {}}}
        feature = {'id': 'f', 'category': 'rRNA', 'member_ids': [0],
                   'counts': {'S': {'exact': 0, 'inexact': 0, 'inclusive': 0}}}
        for evidence in ({'feature_eligible': False, 'status': 'unresolved_parents', 'exact': False},
                         {'feature_eligible': True, 'status': 'unique_locus', 'exact': False}):
            with self.assertRaises(RuntimeError):
                support_reporting([feature], [evidence], {'S': {0: 5}}, summaries)

    def test_unmapped_qc_overlap_and_conservation(self):
        evidence = {0: [{'type': 'secondary_reference_match', 'reference_label': 'host'},
                        {'type': 'secondary_reference_match', 'reference_label': 'bugs'},
                        {'type': 'adapter_exact_substring'}], 1: []}
        summary = {'S': {'total_input_reads': 20, 'unmapped_reads': 10,
                         'below_min_read_length_reads': 2, 'input_scope': 'selected_subset'}}
        qc, labels = unmapped_qc(evidence, {'S': {0: 7, 1: 3}}, summary,
                                 {'host': 'genome', 'bugs': 'microbial'}, True)
        self.assertEqual(labels[0], ['adapter', 'genome', 'microbial'])
        self.assertEqual(qc['S']['partition_counts'], {'multiple_source_evidence': 7, 'no_match_in_supplied_screens': 3})
        self.assertEqual(qc['S']['primary_unmapped_fraction_of_input'], .5)
        self.assertEqual(qc['S']['screens']['genome']['fraction_of_unmapped'], .7)
        self.assertEqual(qc['S']['screens']['other']['status'], 'not_tested')
        self.assertEqual(qc['S']['input_scope'], 'selected_subset')
        with self.assertRaises(RuntimeError):
            unmapped_qc(evidence, {'S': {0: 6, 1: 3}}, summary, {'host': 'genome', 'bugs': 'microbial'}, True)

    def test_unmapped_qc_missing_screens_and_empty(self):
        summary = {'S': {'total_input_reads': 8, 'unmapped_reads': 8,
                         'below_min_read_length_reads': 0, 'input_scope': 'complete'}}
        qc, _ = unmapped_qc({0: []}, {'S': {0: 8}}, summary, {}, False)
        self.assertEqual(qc['S']['partition_counts'], {'not_tested': 8})
        self.assertIsNone(qc['S']['screens']['genome']['matching_reads'])
        summary['S']['unmapped_reads'] = 0
        qc, _ = unmapped_qc({}, {'S': {}}, summary, {'host': 'genome'}, False)
        self.assertEqual(qc['S']['screens']['genome']['status'], 'no_unmapped_reads')
        self.assertIsNone(qc['S']['screens']['genome']['fraction_of_unmapped'])

    def test_diagnostic_source_categories_are_explicit(self):
        args = parser().parse_args(['--diagnostic-reference', 'mouse=ref.fa'])
        self.assertEqual(diagnostic_categories(args), {'mouse': 'other'})
        for spec in ('unknown=genome', 'mouse=virus', 'mouse'):
            args.diagnostic_reference_category = [spec]
            with self.assertRaises(ValueError):
                validate_args(args)
        args.diagnostic_reference_category = ['mouse=genome', 'mouse=genome']
        with self.assertRaises(ValueError):
            validate_args(args)

    def test_ambiguity_qc_weighted_counts_and_boundary(self):
        retention = {'rRNA': {'mapped_reads': 1000, 'raw_locus_ambiguous_reads': 999,
                              'feature_eligible_reads': 1000},
                     'snoRNA': {'mapped_reads': 100, 'raw_locus_ambiguous_reads': 50}}
        metrics, alerts = category_ambiguity_qc('M1', retention, 0.5, 'complete')
        self.assertEqual(len(alerts), 2)
        self.assertEqual(metrics['rRNA']['raw_locus_ambiguous_fraction'], 0.999)
        self.assertEqual(metrics['rRNA']['feature_eligible_fraction'], 1)
        self.assertTrue(metrics['snoRNA']['warning_triggered'])
        self.assertEqual(len(category_ambiguity_qc('M1', retention, 0.501, 'complete')[1]), 1)

    def test_ambiguity_qc_empty_and_unambiguous(self):
        self.assertEqual(category_ambiguity_qc('M1', {}, 0.5, 'complete'), ({}, []))
        metrics, alerts = category_ambiguity_qc('M1', {'rRNA': {}, 'snoRNA': {'mapped_reads': 8}}, 0, 'complete')
        self.assertFalse(alerts)
        self.assertIsNone(metrics['rRNA']['raw_locus_ambiguous_fraction'])
        self.assertEqual(metrics['snoRNA']['raw_locus_ambiguous_fraction'], 0)

    def test_ambiguity_qc_subset_and_cross_category(self):
        retention = {'multi_category_ambiguous': {'mapped_reads': 9, 'raw_locus_ambiguous_reads': 9}}
        before = json.dumps(retention, sort_keys=True)
        metrics, alerts = category_ambiguity_qc('M1', retention, 1, 'selected_subset')
        self.assertEqual(list(metrics), ['multi_category_ambiguous'])
        self.assertEqual(alerts[0]['code'], 'HIGH_CATEGORY_MAPPING_AMBIGUITY')
        self.assertIn('SELECTED SUBSET ONLY', alerts[0]['message'])
        self.assertEqual(alerts[0]['input_scope'], 'selected_subset')
        self.assertEqual(before, json.dumps(retention, sort_keys=True))

    @staticmethod
    def model(refs, headers=None, **kwargs):
        args = parser().parse_args([])
        for key, value in kwargs.items():
            setattr(args, key, value)
        annotations = {name: {'header': (headers or {}).get(name, name + ' snoRNA'),
                              'category': classify_reference((headers or {}).get(name, name + ' snoRNA'))} for name in refs}
        return reference_model(refs, annotations, args), annotations, args

    def test_non_rrna_reference_allowed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'r.fa'
            path.write_text('>s snoRNA\nACGTACGT\n')
            self.assertEqual(load_reference(path)[1]['s']['category'], 'snoRNA')

    def test_reference_category_not_gene_name_substring(self):
        self.assertEqual(classify_reference('GENCODE|a|RRNAD1|protein_coding|'), 'mRNA')
        self.assertEqual(classify_reference('GENCODE|a|rRNA_pseudogene|'), 'pseudogene')
        self.assertEqual(classify_reference('RNAcentral|x rRNA'), 'rRNA')

    def test_palindromic_projection_retains_orientation_ambiguity(self):
        refs = {'m': 'ACGTACGT', 'p': 'TTACGTACGTGG'}
        headers = {'m': 'RIBO|m 28S ribosomal RNA', 'p': 'RIBO|p complete repeating unit'}
        model, annotations, _ = self.model(refs, headers)
        mapping = Mapping(0, {Hit('p', 2, 10, '+', '8M', 0, True)})
        self.assertFalse(resolve_mapping(mapping, model, annotations, 8, .1)['feature_eligible'])

    def test_checked_tool_failure(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(subprocess, 'run', return_value=subprocess.CompletedProcess([], 2)):
                with self.assertRaises(RuntimeError):
                    checked_to_log(['mock-aligner'], directory)

    def test_fastq_whitespace_rejected(self):
        import io
        with self.assertRaises(ValueError):
            list(fastq_records(io.StringIO('@a\n ACGT\n+\nIIIII\n')))

    def test_identical_references_resolve_without_double_count(self):
        refs = {'a': 'ACGTACGT', 'b': 'ACGTACGT'}
        model, annotations, args = self.model(refs)
        mapping = Mapping(0, {Hit('a', 0, 8, '+', '8M', 0, True), Hit('b', 0, 8, '+', '8M', 0, True)})
        evidence = resolve_mapping(mapping, model, annotations, 8, .1)
        self.assertEqual(evidence['status'], 'equivalent_parent_resolved')
        self.assertTrue(evidence['feature_eligible'])
        features, _ = build_features([evidence], ['ACGTACGT'], {'M1': {0: 5}, 'C': {0: 2}}, refs, annotations, 2)
        self.assertEqual(features[0]['counts']['M1']['inclusive'], 5)

    def test_automatic_mature_projection(self):
        refs = {'m': 'AACGTTAC', 'p': 'TT' + 'AACGTTAC' + 'GG'}
        headers = {'m': 'RIBO|m 28S ribosomal RNA', 'p': 'RIBO|p complete repeating unit'}
        model, annotations, args = self.model(refs, headers)
        self.assertEqual(len(model['projections']), 1)
        mapping = Mapping(0, {Hit('m', 1, 7, '+', '6M', 0, True), Hit('p', 3, 9, '+', '6M', 0, True)})
        self.assertEqual(resolve_mapping(mapping, model, annotations, 6, .1)['locations'], [('m', 1, 7, '+')])
        self.assertEqual(projected_hits(Hit('p', 0, 4, '+', '4M', 0, True), model), {('p', 0, 4, '+')})

    def test_reverse_projection_and_invalid_projection(self):
        refs = {'m': 'AACGTTAC', 'p': 'GG' + reverse_complement('AACGTTAC') + 'TT'}
        model, _, _ = self.model(refs, reference_projection=['p=m:2:-'])
        self.assertEqual(projected_hits(Hit('p', 3, 7, '+', '4M', 0, True), model), {('m', 3, 7, '-')})
        with self.assertRaises(ValueError):
            self.model(refs, reference_projection=['p=m:1:+'])

    def test_family_resolution_is_not_coordinate_resolution(self):
        refs = {'a': 'AAAACCCC', 'b': 'AAAAGGGG'}
        model, annotations, _ = self.model(refs, parent_family=['a=family', 'b=family'])
        mapping = Mapping(0, {Hit('a', 0, 4, '+', '4M', 0, True), Hit('b', 0, 4, '+', '4M', 0, True)})
        evidence = resolve_mapping(mapping, model, annotations, 4, .1)
        self.assertEqual(evidence['status'], 'family_resolved_coordinates_ambiguous')
        self.assertFalse(evidence['feature_eligible'])

    def test_inexact_eligible_and_terminal_indel_excluded(self):
        refs = {'a': 'A' * 40}
        model, annotations, _ = self.model(refs)
        for cigar, length, edits, expected in [('40M', 40, 1, True), ('20M1D19M', 39, 1, True),
                                              ('1I39M', 40, 1, False), ('40M', 40, 10, False), ('1S39M', 40, 1, False)]:
            evidence = resolve_mapping(Mapping(-6, {Hit('a', 0, 40, '+', cigar, edits, False)}), model, annotations, length, .1)
            self.assertEqual(evidence['feature_eligible'], expected)

    def test_cluster_tolerance_does_not_chain(self):
        assignment, clusters = endpoint_clusters([0, 2, 4, 6, 8, 10], 2)
        self.assertEqual([(c['minimum'], c['maximum']) for c in clusters], [(0, 4), (6, 10)])
        for value, group in assignment.items():
            self.assertLessEqual(abs(value - clusters[group]['center']), 2)
        self.assertEqual(len(endpoint_clusters([0, 1, 2], 0)[1]), 3)

    def test_paired_features_do_not_merge_overlaps(self):
        refs = {'a': 'A' * 150}
        annotations = {'a': {'category': 'snoRNA'}}
        resolved = [{'feature_eligible': True, 'locations': [('a', s, e, '+')], 'exact': exact}
                    for s, e, exact in [(0, 60, True), (1, 61, False), (40, 100, True)]]
        features, _ = build_features(resolved, ['A' * 60, 'C' + 'A' * 59, 'A' * 59 + 'C'],
                                     {'M1': {0: 4, 1: 5, 2: 6}, 'C': {0: 2, 1: 3}}, refs, annotations, 2)
        self.assertEqual(len(features), 2)
        paired = next(f for f in features if len(f['member_ids']) == 2)
        self.assertEqual(paired['counts']['M1'], {'exact': 4, 'inexact': 5, 'inclusive': 9})
        self.assertEqual(paired['counts']['C']['inclusive'], 5)
        self.assertEqual(sum(f['counts']['M1']['inclusive'] for f in features), 15)
        self.assertIn((paired['start'], paired['end']), [(0, 60), (1, 61)])

    def test_threshold_after_aggregation_and_symmetry(self):
        args = parser().parse_args([])
        args.min_enrichment = 1
        comparison = compare_counts({'M1': 120, 'C': 120}, {'M1': 1_000_000, 'C': 1_000_000}, 'M1', ['C'], args)
        self.assertTrue(comparison['candidate'])
        self.assertEqual(comparison['controls']['C']['fold_enrichment'], 1)
        args.input_scope = 'selected_subset'
        self.assertIsNone(compare_counts({'M1': 120, 'C': 0}, {'M1': 120, 'C': 1}, 'M1', ['C'], args)['candidate'])

    def test_phred_retained_and_normalized(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'q.fq'
            path.write_text('@a\nACGU\n+\nIIII\n@b\nACGT\n+\n!!!!\n')
            quality = {}
            self.assertEqual(load_counts(path, qualities=quality), {'ACGT': 2})
            self.assertEqual(quality['ACGT']['phred_sum'], [40] * 4)
            self.assertEqual(quality['ACGT']['below_q20'], [1] * 4)
            self.assertEqual(quality['ACGT']['records'], 2)

    def test_cigar_mismatch_reverse_positions(self):
        refs = {'a': 'AACCGGTT'}
        seq = 'AACCGGTA'
        events = alignment_differences(seq, Hit('a', 0, 8, '+', '8M', 1, False), refs)
        self.assertEqual(events[0]['read_position'], 7)
        events = alignment_differences(reverse_complement(seq), Hit('a', 0, 8, '-', '8M', 1, False), refs)
        self.assertEqual(events[0]['read_position'], 0)
        events = alignment_differences('AACGGGTT', Hit('a', 0, 8, '+', '2M1D3M1I2M', 2, False), refs)
        self.assertEqual([e['type'] for e in events], ['deletion', 'insertion'])

    def test_selected_tables_only_in_explicit_scope(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'selected.csv'
            path.write_text('sequence,count,RPM_total_mapped\nACGT,10,300\n')
            self.assertEqual(load_counts(path, allow_selected=True), {'ACGT': 10})
            with self.assertRaises(ValueError):
                load_counts(path, counts_complete=True)

    def test_sample_metadata_exact_membership(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'm.csv'
            path.write_text('sample_id,condition,replicate\nb,control,1\na,target,1\n')
            self.assertEqual([r['sample_id'] for r in sample_metadata(path, ['a', 'b'])], ['a', 'b'])
            with self.assertRaises(ValueError):
                sample_metadata(path, ['a'])

    def test_landscape_requires_support_and_is_descriptive(self):
        args = parser().parse_args([])
        feature = {'parent_RNA': 'p', 'strand': '+', 'start': 0, 'end': 100,
                   'counts': {'M1': {'inclusive': 20}, 'C': {'inclusive': 1}}}
        rows = parent_landscape([feature], {'p': 'A' * 100}, ['M1', 'C'], args)
        self.assertEqual(rows[0]['descriptive_pattern'], 'full_length_compatible')
        self.assertEqual(rows[1]['descriptive_pattern'], 'insufficient_support')


# %% Real-tool integration fixtures
def test_with_tools():
    """Four complete real-tool fixtures: biology units, controls, duplication, subsets."""
    import random
    randomizer = random.Random(1701)
    random_dna = lambda n: ''.join(randomizer.choice('ACGT') for _ in range(n))
    reference, sno, unknown = random_dna(500), random_dna(130), random_dna(60)
    variant = reference[20:80]
    variant = variant[:30] + next(b for b in 'ACGT' if b != variant[30]) + variant[31:]
    tools = require_tools()
    with tempfile.TemporaryDirectory(prefix='csrna-integration-') as directory:
        directory = Path(directory)
        ref, diagnostic = directory / 'ref.fa', directory / 'diagnostic.fa'
        ref_text = '>RIBO|mature 28S ribosomal RNA\n' + reference + '\n>RIBO|precursor complete repeating unit\n' + 'G' * 40 + reference + 'C' * 40 + '\n>sno snoRNA\n' + sno + '\n'
        ref.write_text(ref_text)
        diagnostic.write_text('>diagnostic mock sequence\n' + unknown + '\n')
        reads = {'M1': [reference[20:80]] * 8 + [variant] * 8 + [reference[60:120]] * 8 +
                       [reverse_complement(reference[250:310])] * 4 + [sno[15:75]] * 20 + [unknown] * 4,
                 'BMDM': [reference[20:80]] * 4 + [reference[200:260]] * 48,
                 'BMDC': [reference[350:410]] * 52}
        cli = ['--reference', str(ref), '--output', str(directory / 'result.json'), '--threads', '1',
               '--min-enrichment', '2', '--diagnostic-reference', 'mock_genome=' + str(diagnostic),
               '--diagnostic-reference-category', 'mock_genome=genome', '--adapter-sequence', unknown[:8]]
        for name, sample_reads in reads.items():
            path = directory / f'{name}.fq'
            path.write_text(''.join(f'@{i}\n{seq}\n+\n{"I" * len(seq)}\n' for i, seq in enumerate(sample_reads)))
            cli += ['--sample', f'{name}={path}']
        args = parser().parse_args(cli)
        validate_args(args)
        def run():
            with tempfile.TemporaryDirectory(dir=directory) as scratch:
                return analyze(args, tools, scratch)
        report = run()
        assert {a['sample'] for a in report['qc_alerts']} == {'M1', 'BMDM', 'BMDC'}
        assert all(a['category'] == 'rRNA' and a['raw_locus_ambiguous_fraction'] == 1
                   and a['feature_eligible_fraction'] == 1 for a in report['qc_alerts'])
        assert not report['samples']['M1']['category_ambiguity_qc']['snoRNA']['warning_triggered']
        atomic_report(report, args.output)
        feature = next(f for f in report['features'] if f['parent_RNA'] == 'RIBO|mature' and f['start'] == 20)
        assert feature['counts']['M1'] == {'exact': 8, 'inexact': 8, 'inclusive': 16}
        landscape = report['landscape_comparison']
        assert len(landscape['feature_table']) == len(report['features'])
        assert any(f['category'] == 'snoRNA' for f in landscape['feature_table'])
        assert len(landscape['pairwise_contrasts']) == 3
        rows = landscape['pairwise_contrasts'][0]['rankings']['inclusive']
        pair_feature = next(row for row in rows if row['feature_id'] == feature['id'])
        assert pair_feature['numerator_count'] == 16 and pair_feature['denominator_count'] == 4
        assert any(row['detection'] == 'denominator_only' for row in rows)
        support = report['support_reporting']['samples']['M1']['counts']
        assert support['exact_only'] == 40 and support['accepted_inexact'] == 8 and support['inclusive'] == 48
        assert support['raw_unique_inexact'] == 0 and support['equivalence_resolved_inexact'] == 8
        assert support['raw_unique_exact'] == 20 and support['equivalence_resolved_exact'] == 20
        assert feature['counts']['BMDM']['inclusive'] == 4
        assert feature['candidate'] and len(feature['observed_variant_folds']) == 2
        assert any(f['category'] == 'snoRNA' and f['counts']['M1']['inclusive'] == 20 for f in report['features'])
        assert len([f for f in report['features'] if f['parent_RNA'] == 'RIBO|mature' and f['strand'] == '+' and f['counts']['M1']['inclusive']]) == 2
        assert all(sum(row[i] for row in report['count_matrices']['inclusive']['data']) == (48 if sample == 'M1' else 52)
                   for i, sample in enumerate(report['count_matrices']['inclusive']['sample_ids']))
        unresolved = next(r for r in report['fragments'] if r['sequence'].replace('U', 'T') == unknown)
        assert unresolved['counts']['M1'] == 4 and len(unresolved['diagnostic_evidence']) == 2
        assert unresolved['unmapped_diagnostic_labels'] == ['adapter', 'genome']
        assert report['unmapped_qc']['M1']['partition_counts'] == {'multiple_source_evidence': 4}
        assert report['unmapped_qc']['M1']['screens']['genome']['matching_reads'] == 4
        assert report['unmapped_qc']['M1']['screens']['microbial']['status'] == 'not_tested'
        assert report['samples']['M1']['unmapped_reads'] == 4
        assert report['samples']['M1']['total_input_reads'] == 52
        assert not list(directory.glob('*.ps'))
        # Shifted control fragments are different species, but remain a control-overlap concern.
        shifted = reference[40:100]
        (directory / 'BMDM.fq').write_text(''.join(f'@{i}\n{shifted}\n+\n{"I" * 60}\n' for i in range(52)))
        args.output = str(directory / 'second.json')
        shifted_report = run()
        f = next(f for f in shifted_report['features'] if f['id'] == feature['id'])
        assert f['counts']['BMDM']['inclusive'] == 0 and f['control_overlap_screen']['BMDM']['count'] == 52
        assert f['candidate']  # Descriptive species enrichment, explicitly not probe specificity.
        # Adding a duplicate reference must not change the feature counts or identity.
        ref.write_text(ref_text + '>RIBO|mature_copy 28S ribosomal RNA\n' + reference + '\n')
        duplicate_report = run()
        assert duplicate_report['count_matrices'] == shifted_report['count_matrices']
        args.input_scope = 'selected_subset'
        subset_report = run()
        assert all(a['input_scope'] == 'selected_subset' and 'SELECTED SUBSET ONLY' in a['message']
                   for a in subset_report['qc_alerts'])
        assert not subset_report['candidate_feature_ids'] and not subset_report['fold_change_ranking']
        assert all(f['candidate'] is None for f in subset_report['features'])
        assert subset_report['count_matrices'] == duplicate_report['count_matrices']
        assert not subset_report['landscape_comparison']['pairwise_contrasts']
        assert subset_report['landscape_comparison']['feature_table'] == duplicate_report['landscape_comparison']['feature_table']
        # No-control paired manifest reproduces exactly the single-end R1 counts.
        manifest_path = directory / 'paired.csv'
        with manifest_path.open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['sample_id', 'r1', 'r2', 'condition', 'replicate', 'replicate_type', 'pcr_cycles'])
            for i, name in enumerate(reads, 1):
                mate = directory / f'{name}.r2.fq'
                with open_text(directory / f'{name}.fq') as source:
                    records = list(fastq_records(source, with_name=True))
                mate.write_text(''.join(f'@{header}\n{reverse_complement(seq)}\n+\n{quality[::-1]}\n' for header, seq, quality in records))
                writer.writerow([name, f'{name}.fq', mate.name, 'example', str(i), 'unknown', 15])
        args.sample = []
        args.manifest = str(manifest_path)
        args.pair_handling = 'r1-only'
        args.analysis_mode = 'exploratory'
        args.input_scope = 'full_library'
        args.report_detail = 'summary'
        args.quality_detail = 'counts-only'
        validate_args(args)
        paired_report = run()
        assert paired_report['count_matrices'] == subset_report['count_matrices']
        assert paired_report['endpoint_profiles'] == subset_report['endpoint_profiles']
        assert paired_report['parent_landscape'] == subset_report['parent_landscape']
        assert not paired_report['controls'] and paired_report['target'] is None
        assert not paired_report['candidate_feature_ids'] and not paired_report['fold_change_ranking']
        assert not paired_report['fragments']
        assert not paired_report['landscape_comparison']['pairwise_contrasts']
        assert all(row['validated_pairs'] == 52 and row['analyzed_r1_reads'] == 52 for row in paired_report['paired_input_qc'].values())
        for feature in paired_report['features']:
            assert feature['candidate'] is None and feature['comparison']['inclusive']['status'] == 'disabled_no_controls'
        args.pilot_pairs = 7
        args.input_scope = 'selected_subset'
        sampled_report = run()
        assert all(row['total_input_reads'] == 7 for row in sampled_report['samples'].values())
        assert all(row['validated_pairs'] == 52 for row in sampled_report['paired_input_qc'].values())
    print('Six real-tool workflows passed; all synthetic inputs/results removed.')


class AccuracyPerformanceTests(unittest.TestCase):
    def test_fast_validation_equivalence(self):
        import io
        # Exhaustive byte-range and selected Unicode cases, plus valid long sequences.
        for ref in (False, True):
            alphabet = 'ACGTRYSWKMBDHVN' if ref else 'ACGTN'
            for value in [chr(i) for i in range(512)] + ['acgu', ' ACGT ', 'ACGT' * 1000, 'A G', 'A\nG']:
                expected = value.strip().upper().replace('U', 'T')
                valid = bool(expected) and all(c in alphabet for c in expected)
                if valid:
                    self.assertEqual(clean_sequence(value, ref), expected)
                else:
                    with self.assertRaises(ValueError): clean_sequence(value, ref)
        for code in list(range(256)) + [1000, 65535]:
            quality = chr(code) * 4
            text = '@q\nACGT\n+\n' + quality + '\n'
            if 33 <= code <= 126:
                self.assertEqual(list(fastq_records(io.StringIO(text), True)), [('ACGT', quality)])
            else:
                with self.assertRaises(ValueError): list(fastq_records(io.StringIO(text), True))

    def test_combined_limit_stops_before_mapping(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref=root/'ref.fa';ref.write_text('>r\nACGTACGT\n')
            a=root/'a.fa';a.write_text('>q\nAAAA\n')
            b=root/'b.fa';b.write_text('>q\nCCCC\n')
            args=parser().parse_args(['--sample',f'A={a}','--sample',f'B={b}', '--reference',str(ref),
                                     '--analysis-mode','exploratory','--max-total-distinct-sequences','1'])
            validate_args(args)
            with patch.dict(analyze.__globals__, {'align_sequences': lambda *a: self.fail('Mapping must not start')}):
                with self.assertRaisesRegex(ValueError,'Combined distinct-sequence limit'):
                    analyze(args,{},d)
            self.assertFalse((root/'csrna_results.json').exists())
        args=parser().parse_args(['--sample','A=a.fa','--reference','r.fa','--max-total-distinct-sequences','-1'])
        with self.assertRaises(ValueError):validate_args(args)


class ReplicateInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.r1, self.r2 = self.root / 'one.fastq.gz', self.root / 'two.fastq.gz'
        self.write_pair(100)
        self.manifest = self.root / 'samples.csv'
        self.manifest.write_text('sample_id,r1,r2,condition,replicate,replicate_type,pcr_cycles\nS,one.fastq.gz,two.fastq.gz,group,1,unknown,15\n')

    def tearDown(self):
        self.temp.cleanup()

    def write_pair(self, n):
        for mate, path in enumerate((self.r1, self.r2), 1):
            with gzip.open(path, 'wt') as h:
                for i in range(n):
                    h.write(f'@read{i}/{mate}\nACGTACGT\n+\nIIIIIIII\n')

    def args(self, extra=()):
        return parser().parse_args(['--manifest', str(self.manifest), '--reference', 'unused.fa',
                                   '--pair-handling', 'r1-only', '--analysis-mode', 'exploratory', *extra])

    def test_manifest_identity_and_metadata(self):
        row = read_manifest(self.manifest)[0]
        self.assertEqual(row['r1'], str(self.r1))
        self.assertEqual(row['replicate_type'], 'unknown')
        self.assertEqual(len(list(paired_r1_records(self.r1, self.r2))), 100)
        self.manifest.write_text(self.manifest.read_text().replace('two.fastq.gz', 'one.fastq.gz'))
        with self.assertRaises(ValueError): read_manifest(self.manifest)

    def test_pair_corruption_and_mate_identity(self):
        for text in ('@wrong/2\nACGT\n+\nIIII\n', '@read0/1\nACGT\n+\nIIII\n',
                     '@read0/2\nAXGT\n+\nIIII\n', '@read0/2\nACGT\n+\nIII\n'):
            with gzip.open(self.r2, 'wt') as h: h.write(text)
            with self.assertRaises(ValueError): list(paired_r1_records(self.r1, self.r2))

    def test_count_mismatch_and_illumina_headers(self):
        with gzip.open(self.r2, 'wt') as h: h.write('@read0 2:N:0:INDEX\nACGT\n+\nIIII\n')
        with self.assertRaisesRegex(ValueError, 'counts differ'): list(paired_r1_records(self.r1, self.r2))
        with gzip.open(self.r1, 'wt') as h: h.write('@read0 1:N:0:INDEX\nACGT\n+\nIIII\n')
        self.assertEqual(list(paired_r1_records(self.r1, self.r2)), [('ACGT', 'IIII')])

    def test_uniform_sampling_reproducible_and_complete_validation(self):
        args = self.args(['--pilot-pairs', '10', '--input-scope', 'selected_subset'])
        validate_args(args)
        rows = read_manifest(self.manifest)
        paths, audit = prepare_manifest_inputs(rows, args, self.root)
        first = paths['S'].read_text()
        self.assertEqual(audit['S']['validated_pairs'], 100)
        self.assertEqual(audit['S']['analyzed_r1_reads'], 10)
        ids = [int(line.split('_')[-1]) for line in first.splitlines() if line.startswith('@')]
        self.assertTrue(any(i > 10 for i in ids))
        prepare_manifest_inputs(rows, args, self.root)
        self.assertEqual(first, paths['S'].read_text())
        args.sampling_seed += 1
        prepare_manifest_inputs(rows, args, self.root)
        self.assertNotEqual(first, paths['S'].read_text())
        with gzip.open(self.r2, 'at') as h: h.write('@extra/2\nACGT\n+\nIIII\n')
        with self.assertRaises(ValueError): prepare_manifest_inputs(rows, args, self.root)

    def test_small_and_full_libraries(self):
        args = self.args(['--pilot-pairs', '200', '--input-scope', 'selected_subset'])
        rows = read_manifest(self.manifest)
        paths, audit = prepare_manifest_inputs(rows, args, self.root)
        self.assertEqual(audit['S']['analyzed_r1_reads'], 100)
        self.assertEqual(sum(load_counts(paths['S']).values()), 100)
        args.pilot_pairs = 0
        paths, audit = prepare_manifest_inputs(rows, args, self.root)
        self.assertEqual(paths['S'], self.r1)
        self.assertEqual(audit['S']['sampling'], 'all_pairs')

    def test_unsafe_cli_combinations(self):
        for extra in (['--pilot-pairs', '2'], ['--pilot-pairs', '-1'], ['--control', 'S'],
                      ['--sample', 'S=file.fastq'], ['--sample-metadata', 'metadata.csv']):
            with self.assertRaises(ValueError): validate_args(self.args(extra))
        args = self.args(); args.pair_handling = None
        with self.assertRaises(ValueError): validate_args(args)
        validate_args(self.args())

    def test_no_control_comparison(self):
        result = compare_counts({'S': 10}, {'S': 10}, None, [], self.args())
        self.assertIsNone(result['candidate'])
        self.assertEqual(result['status'], 'disabled_no_controls')


# %% Command-line entry point
def main(argv=None):
    args = parser().parse_args(argv)
    if args.self_test or args.self_test_tools:
        suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(cls) for cls in (PipelineTests, RedesignTests, ReplicateInputTests, AccuracyPerformanceTests))
        if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
            return 1
        if args.self_test_tools:
            test_with_tools()
        return 0
    validate_args(args)
    tools = require_tools()
    with tempfile.TemporaryDirectory(prefix="csrna-", dir=args.temp_dir) as workdir:
        report = analyze(args, tools, workdir)
    atomic_report(report, args.output, args.overwrite)
    print(f"Complete: {len(report['features'])} fragment features; {len(report['candidate_feature_ids'])} descriptive candidates; {Path(args.output).resolve()}")
    return 0


if __name__ == "__main__" and "ipykernel" not in sys.modules:
    try:
        raise SystemExit(main())
    except (ValueError, OSError, EOFError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"Pipeline failed; no new report published: {error}", file=sys.stderr)
        raise SystemExit(1)
