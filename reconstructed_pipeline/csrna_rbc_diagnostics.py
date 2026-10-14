#!/usr/bin/env python3
"""RBC-only post-analysis. Run beside csrna_pipeline.py; never modifies its report.

SQLite holds full read-pair membership and sequence counts on disk. R1 and R2 are
aligned independently in both modes. No insert reconstruction, trimming, molecule
count, contamination verdict, or native endpoint inference is performed.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time

from csrna_pipeline import fastq_records, open_text, fasta_records, clean_sequence, classify_reference, sha256, reverse_complement

VERSION = '1.0.0'
SCHEMA = 'csrna-rbc-diagnostics/1'
MODES = ('end_to_end', 'local')


def atomic_json(path, value):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(value, handle, allow_nan=False, separators=(',', ':'))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@contextmanager
def lock(directory):
    directory.mkdir(parents=True, exist_ok=True)
    file = directory / '.lock'
    try:
        fd = os.open(file, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise ValueError(f'Workspace/index is locked: {file}. Check for an active process before removing a stale lock.')
    try:
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)
        yield
    finally:
        file.unlink()


def complexity(sequence):
    n = len(sequence)
    bases = Counter(sequence)
    entropy = -sum((c/n) * math.log2(c/n) for b, c in bases.items() if b in 'ACGT')
    runs = {b: max((len(m.group()) for m in re.finditer(b + '+', sequence)), default=0) for b in 'ACGT'}
    terminal = {end: {b: len(sequence) - len(sequence.lstrip(b) if end == 'five_prime' else sequence.rstrip(b))
                      for b in 'ACGT'} for end in ('five_prime', 'three_prime')}
    # Fraction of distinct 3-mers: always include the length, because this statistic is length-dependent.
    kmers = [sequence[i:i+3] for i in range(max(0, n-2))]
    kmer_fraction = len(set(kmers))/len(kmers) if kmers else None
    return {'length': n, 'gc_fraction': (bases['G']+bases['C'])/n, 'n_fraction': bases['N']/n,
            'entropy_bits': entropy, 'longest_runs': runs, 'terminal_runs': terminal,
            'longest_run_fraction': max(runs.values())/n, 'distinct_3mer_fraction': kmer_fraction}


def complexity_labels(metrics, args):
    labels = []
    if metrics['entropy_bits'] < args.entropy_threshold:
        labels.append('low_base_entropy')
    if metrics['distinct_3mer_fraction'] is not None and metrics['distinct_3mer_fraction'] < args.kmer_threshold:
        labels.append('low_3mer_complexity')
    if max(metrics['longest_runs'].values()) >= args.homopolymer_length:
        labels.append('homopolymer_rich')
    for base in 'ACGT':
        if any(v[base] >= args.homopolymer_length for v in metrics['terminal_runs'].values()):
            labels.append('terminal_poly_' + base)
    if metrics['n_fraction'] > 0:
        labels.append('contains_N')
    return labels


def pair_records(r1, r2):
    def identity(header, mate):
        parts = header.split()
        token = parts[0]
        if token.endswith(('/1', '/2')):
            if token[-1] != str(mate):
                raise ValueError('Incorrect FASTQ mate suffix')
            token = token[:-2]
        if len(parts) > 1 and re.match(r'^[12]:', parts[1]) and parts[1][0] != str(mate):
            raise ValueError('Incorrect Illumina mate field')
        return token
    with open_text(r1) as a, open_text(r2) as b:
        for ordinal, pair in enumerate(zip_longest(fastq_records(a, with_name=True), fastq_records(b, with_name=True)), 1):
            left, right = pair
            if left is None or right is None or identity(left[0], 1) != identity(right[0], 2):
                raise ValueError(f'Mate count/identity mismatch at pair {ordinal}')
            yield ordinal, left, right


def selected_manifest(path, report):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle)
        if len(reader.fieldnames or []) != len(set(reader.fieldnames or [])) or not {'sample_id', 'condition', 'r1', 'r2'} <= set(reader.fieldnames or []):
            raise ValueError('Manifest requires sample_id,condition,r1,r2')
        rows = list(reader)
    seen = set()
    files = set()
    selected = []
    meta = {x['sample_id']: x for x in report.get('sample_metadata', [])}
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError('Malformed manifest row')
        row = {k:v.strip() for k,v in row.items()}
        name = row['sample_id']
        if name in seen or not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise ValueError('Manifest sample IDs must be unique and filesystem-safe')
        seen.add(name)
        if row['condition'] != 'RBC':
            continue
        if name not in report['samples'] or meta.get(name, {}).get('condition') != 'RBC':
            raise ValueError(f'{name} must be labeled RBC in both manifest and main report')
        for mate in ('r1', 'r2'):
            p = Path(row[mate]).expanduser()
            p = (p if p.is_absolute() else Path(path).resolve().parent / p).resolve(strict=True)
            if p in files:
                raise ValueError('A FASTQ cannot represent two mates or samples')
            files.add(p)
            row[mate] = str(p)
        selected.append(row)
    if not selected:
        raise ValueError('No RBC samples: manifest condition must be exactly RBC')
    return selected


def connect(path):
    db = sqlite3.connect(path)
    db.execute('PRAGMA cache_size=-32768')
    db.execute('PRAGMA temp_store=FILE')
    db.executescript('''
      CREATE TABLE IF NOT EXISTS seq(id INTEGER PRIMARY KEY, sequence TEXT UNIQUE NOT NULL);
      CREATE TABLE IF NOT EXISTS counts(sample TEXT, mate INTEGER, qid INTEGER, n INTEGER, qsum REAL,
        PRIMARY KEY(sample,mate,qid));
      CREATE TABLE IF NOT EXISTS pairs(sample TEXT, ordinal INTEGER, r1 INTEGER, r2 INTEGER, PRIMARY KEY(sample,ordinal));
      CREATE TABLE IF NOT EXISTS done(stage TEXT PRIMARY KEY, data TEXT);
      CREATE TABLE IF NOT EXISTS align(qid INTEGER, label TEXT, mode TEXT, score INTEGER,
        PRIMARY KEY(qid,label,mode));
      CREATE TABLE IF NOT EXISTS hit(qid INTEGER,label TEXT,mode TEXT,ref TEXT,start INTEGER,strand TEXT,
        cigar TEXT,aligned INTEGER,edits INTEGER,left_clip INTEGER,right_clip INTEGER,passed INTEGER,
        PRIMARY KEY(qid,label,mode,ref,start,strand,cigar));
      CREATE TABLE IF NOT EXISTS refs(label TEXT,ref TEXT,category TEXT,length INTEGER,PRIMARY KEY(label,ref));
    ''')
    return db


def done(db, stage):
    row = db.execute('SELECT data FROM done WHERE stage=?', (stage,)).fetchone()
    return json.loads(row[0]) if row else None


def finish(db, stage, data):
    db.execute('INSERT OR REPLACE INTO done VALUES (?,?)', (stage,json.dumps(data)))
    db.commit()


def ingest(db, row, args, work):
    stage = 'input:' + row['sample_id']
    if done(db, stage):
        previous = done(db, stage)
        for retained in previous['retained_pilot_fastqs']:
            if not Path(retained['path']).is_file() or sha256(retained['path']) != retained['sha256']:
                raise ValueError('Retained pilot file changed or is missing; use a new workspace')
        return previous
    sample = row['sample_id']
    for table in ('pairs', 'counts'):
        db.execute(f'DELETE FROM {table} WHERE sample=?', (sample,))
    db.commit()
    seed = int(hashlib.sha256(f'{args.sampling_seed}:{sample}'.encode()).hexdigest(),16)
    rng, reservoir, total = random.Random(seed), [], 0
    def add(record):
        ordinal,left,right = record
        ids=[]
        for mate,read in enumerate((left,right),1):
            sequence,quality = read[1:]
            db.execute('INSERT OR IGNORE INTO seq(sequence) VALUES (?)',(sequence,))
            qid=db.execute('SELECT id FROM seq WHERE sequence=?',(sequence,)).fetchone()[0]
            ids.append(qid)
            meanq=sum(ord(c)-33 for c in quality)/len(quality)
            db.execute('INSERT INTO counts VALUES (?,?,?,1,?) ON CONFLICT(sample,mate,qid) DO UPDATE SET n=n+1,qsum=qsum+excluded.qsum',
                       (sample,mate,qid,meanq))
        db.execute('INSERT INTO pairs VALUES (?,?,?,?)',(sample,ordinal,*ids))
    for record in pair_records(row['r1'],row['r2']):
        total=record[0]
        if args.pilot_pairs:
            if total <= args.pilot_pairs:
                reservoir.append(record)
            else:
                index=rng.randrange(total)
                if index < args.pilot_pairs: reservoir[index]=record
        else:
            add(record)
        if total % 1000000 == 0:
            db.commit()
            print(f'INPUT {sample}: validated {total:,} pairs',file=sys.stderr,flush=True)
    if not total: raise ValueError(f'Empty paired library: {sample}')
    retained = []
    if args.pilot_pairs:
        reservoir.sort(key=lambda x:x[0])
        for mate in (1,2):
            p=work/f'{sample}.pilot.R{mate}.fastq.gz'
            temporary=Path(str(p)+'.partial')
            with gzip.open(temporary,'wt') as h:
                for record in reservoir:
                    header,sequence,quality=record[mate]
                    h.write(f'@{header}\n{sequence}\n+\n{quality}\n')
            os.replace(temporary,p)
            retained.append({'path':str(p),'sha256':sha256(p)})
        for record in reservoir: add(record)
    used=min(total,args.pilot_pairs) if args.pilot_pairs else total
    assert db.execute('SELECT count(*) FROM pairs WHERE sample=?',(sample,)).fetchone()[0]==used
    for mate in (1,2):
        assert db.execute('SELECT sum(n) FROM counts WHERE sample=? AND mate=?',(sample,mate)).fetchone()[0]==used
    result={'validated_pairs':total,'analyzed_pairs':used,'retained_pilot_fastqs':retained}
    finish(db,stage,result)
    return result


def checked(command, log):
    with Path(log).open('w') as h:
        p=subprocess.run(command,stdout=h,stderr=subprocess.STDOUT)
    if p.returncode: raise RuntimeError(f'Command failed ({p.returncode}); see {log}')


def index_reference(db,label,path,args,tools,work):
    digest=sha256(path)
    key=hashlib.sha256((digest+tools['bowtie2-build']['version']).encode()).hexdigest()[:32]
    cache=Path(args.index_cache).expanduser().resolve()/key
    with lock(cache):
        fasta=cache/'reference.fa'
        prefix=cache/'index'
        marker=cache/'complete.json'
        complete=json.loads(marker.read_text()) if marker.exists() else None
        # Check every expected index file, not just an empty success marker.
        suffixes=['1','2','3','4','rev.1','rev.2']
        valid_index=any(all(Path(str(prefix)+'.'+s+'.'+ext).is_file() and Path(str(prefix)+'.'+s+'.'+ext).stat().st_size>0 for s in suffixes)
                        for ext in ('bt2','bt2l'))
        if not complete or not valid_index:
            if marker.exists(): marker.unlink()
            with open_text(path) as source, fasta.open('w') as dest:
                # A disk table enforces unique IDs without retaining the reference in RAM.
                db.execute('DELETE FROM refs WHERE label=?',(label,))
                for header,sequence in fasta_records(source):
                    name=header.split()[0];sequence=clean_sequence(sequence,reference=True)
                    db.execute('INSERT INTO refs VALUES (?,?,?,?)',(label,name,classify_reference(header),len(sequence)))
                    dest.write(f'>{name}\n{sequence}\n')
            db.commit()
            checked([tools['bowtie2-build']['path'],'--threads',str(args.threads),str(fasta),str(prefix)],cache/'build.log')
            if sha256(path)!=digest: raise ValueError('Reference changed during indexing')
            atomic_json(marker,{'sha256':digest,'builder':tools['bowtie2-build']})
        else:
            db.execute('DELETE FROM refs WHERE label=?',(label,))
            with open_text(path) as source:
                for header,sequence in fasta_records(source):
                    sequence=clean_sequence(sequence,reference=True)
                    db.execute('INSERT INTO refs VALUES (?,?,?,?)',(label,header.split()[0],classify_reference(header),len(sequence)))
            db.commit()
    return str(prefix)


def parse_hit(fields,sequence,args):
    cigar=fields[5]
    pieces=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in pieces)!=cigar: raise ValueError('Invalid alignment CIGAR')
    ops=[(int(n),op) for n,op in pieces]
    if not ops or any(n<=0 or op in 'NHP' for n,op in ops): raise ValueError('Unsupported alignment CIGAR')
    if any(op=='S' and i not in (0,len(ops)-1) for i,(_,op) in enumerate(ops)): raise ValueError('Internal soft clip')
    if sum(n for n,op in ops if op in 'MIS=X')!=len(sequence): raise ValueError('Alignment/query length mismatch')
    tags={t.split(':',2)[0]:t.split(':',2)[2] for t in fields[11:] if t.count(':')>=2}
    if 'NM' not in tags or 'AS' not in tags: raise ValueError('Missing Bowtie2 alignment evidence')
    left=ops[0][0] if ops[0][1]=='S' else 0
    right=ops[-1][0] if ops[-1][1]=='S' else 0
    # Convert clipping to original sequenced-read orientation.
    reverse=bool(int(fields[1])&16)
    if reverse: left,right=right,left
    aligned=len(sequence)-left-right
    edits=int(tags['NM'])
    passed=aligned>=args.min_aligned_length and aligned/len(sequence)>=args.min_aligned_fraction and edits/max(1,aligned)<=args.max_edit_rate
    return int(tags['AS']),fields[2],int(fields[3])-1,'-' if reverse else '+',cigar,aligned,edits,left,right,int(passed)


def align(db,label,mode,prefix,args,tools,work,query):
    stage=f'align:{label}:{mode}'
    if done(db,stage): return done(db,stage)
    for table in ('align','hit'):
        db.execute(f'DELETE FROM {table} WHERE label=? AND mode=?',(label,mode))
    db.commit()
    command=[tools['bowtie2']['path'],'-x',prefix,'-f','-U',str(query),'-a','--seed','0','--reorder','-p',str(args.threads)]
    command+=['--end-to-end','--very-sensitive'] if mode=='end_to_end' else ['--local','--very-sensitive-local']
    logfile=work/f'{label}.{mode}.log'
    count=0;started=time.monotonic()
    with logfile.open('w') as log:
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=log,text=True)
        try:
            for line in p.stdout:
                if line.startswith('@'): continue
                fields=line.rstrip('\n').split('\t');qid=int(fields[0][1:]);flag=int(fields[1])
                if flag & 2048: raise ValueError('Unexpected supplementary alignment')
                row=db.execute('SELECT sequence FROM seq WHERE id=?',(qid,)).fetchone()
                if row is None: raise ValueError('Unknown alignment query')
                db.execute('INSERT OR IGNORE INTO align VALUES (?,?,?,NULL)',(qid,label,mode))
                if not flag&4:
                    score,*hit=parse_hit(fields,row[0],args)
                    old=db.execute('SELECT score FROM align WHERE qid=? AND label=? AND mode=?',(qid,label,mode)).fetchone()[0]
                    if old is None or score>old:
                        db.execute('DELETE FROM hit WHERE qid=? AND label=? AND mode=?',(qid,label,mode))
                        db.execute('UPDATE align SET score=? WHERE qid=? AND label=? AND mode=?',(score,qid,label,mode))
                    if old is None or score>=old:
                        db.execute('INSERT OR IGNORE INTO hit VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',(qid,label,mode,*hit))
                count+=1
                if count%10000==0: db.commit()
            if p.wait(): raise RuntimeError(f'Bowtie2 failed; see {logfile}')
        finally:
            p.stdout.close()
            if p.poll() is None: p.kill()
            p.wait()
    expected=db.execute('SELECT count(*) FROM seq WHERE id IN (SELECT qid FROM counts)').fetchone()[0]
    if db.execute('SELECT count(*) FROM align WHERE label=? AND mode=?',(label,mode)).fetchone()[0]!=expected:
        raise ValueError('Missing alignment query records')
    result={'seconds':time.monotonic()-started,'queries':expected,'command':command,'log':str(logfile)}
    finish(db,stage,result)
    return result


def evidence(db,qid,label,mode):
    rows=db.execute('SELECT ref,start,strand,cigar,aligned,edits,left_clip,right_clip,passed FROM hit WHERE qid=? AND label=? AND mode=?',
                    (qid,label,mode))
    n=0;passed=0;parent=None;category=None;mincoverage=None;maxcoverage=None;example=None;reported_example=None
    for ref,start,strand,cigar,aligned,edits,left,right,ok in rows:
        n+=1
        if reported_example is None:
            reported_example={'reference':ref,'start':start,'strand':strand,'cigar':cigar,'aligned_bases':aligned,'edits':edits,'five_prime_clipped_bases':left,'three_prime_clipped_bases':right}
        if not ok: continue
        passed+=1
        parent=ref if passed==1 else parent if parent==ref else '__ambiguous__'
        cat=db.execute('SELECT category FROM refs WHERE label=? AND ref=?',(label,ref)).fetchone()[0]
        category=cat if passed==1 else category if category==cat else 'multi_category_ambiguous'
        mincoverage=aligned if mincoverage is None else min(mincoverage,aligned)
        maxcoverage=aligned if maxcoverage is None else max(maxcoverage,aligned)
        if example is None:example={'reference':ref,'start':start,'strand':strand,'cigar':cigar,'aligned_bases':aligned,'edits':edits,'five_prime_clipped_bases':left,'three_prime_clipped_bases':right}
    return {'reported_best_loci':n,'accepted_best_loci':passed,'parent':None if parent=='__ambiguous__' else parent,
            'category':category,'min_aligned_bases':mincoverage,'max_aligned_bases':maxcoverage,'example_alignment':example,
            'reported_example_alignment':reported_example}


def mode_state(end,local):
    a,b=end['accepted_best_loci']>0,local['accepted_best_loci']>0
    if a and b:
        if end['parent'] and local['parent']:
            return 'both_same_parent' if end['parent']==local['parent'] else 'both_different_parent'
        return 'both_parent_ambiguous'
    return 'end_to_end_only' if a else 'local_only' if b else 'neither_accepted'


def summarize(db,args,rows,work):
    # Small summary JSON; full per-sequence evidence streams to compressed TSV on disk.
    summaries={};pair_summary={};started=time.monotonic()
    screens=[spec.split('=',1)[0] for spec in args.screen_reference]
    adapters=[clean_sequence(s) for s in args.adapter_sequence]
    db.execute('CREATE TABLE IF NOT EXISTS accepted(qid INTEGER,mode TEXT,parent TEXT,ok INTEGER,PRIMARY KEY(qid,mode))')
    db.execute('DELETE FROM accepted')
    # Compute sequence evidence once and store it on disk for reused R1/R2 observations.
    db.execute('CREATE TABLE IF NOT EXISTS evidence(qid INTEGER PRIMARY KEY,data TEXT)')
    db.execute('DELETE FROM evidence')
    for qid,sequence in db.execute('SELECT id,sequence FROM seq WHERE id IN (SELECT qid FROM counts)'):
        metrics=complexity(sequence);labels=complexity_labels(metrics,args)
        modes={m:evidence(db,qid,'primary',m) for m in MODES}
        for mode,v in modes.items():
            db.execute('INSERT INTO accepted VALUES (?,?,?,?)',(qid,mode,v['parent'],int(v['accepted_best_loci']>0)))
        origin=[]
        for label in screens:
            hit=evidence(db,qid,label,'local')
            if hit['accepted_best_loci']:origin.append({'label':label,'evidence':hit})
        if any(a in sequence or reverse_complement(a) in sequence for a in adapters):labels.append('adapter_exact_substring')
        data={'complexity':metrics,'labels':labels,'alignment':modes,'mode_state':mode_state(modes['end_to_end'],modes['local']),
              'origin_screens':origin}
        db.execute('INSERT INTO evidence VALUES (?,?)',(qid,json.dumps(data,separators=(',',':'))))
    db.commit()
    details=work/'sequence_evidence.tsv.gz'
    temporary=Path(str(details)+'.partial')
    with gzip.open(temporary,'wt',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t')
        writer.writerow(['sample','mate','sequence_id','sequence','read_count','mean_read_phred','evidence_json'])
        for sample,mate,qid,n,qsum,sequence,data in db.execute('SELECT c.sample,c.mate,c.qid,c.n,c.qsum,s.sequence,e.data FROM counts c JOIN seq s ON s.id=c.qid JOIN evidence e ON e.qid=c.qid ORDER BY c.sample,c.mate,c.qid'):
            key=f'{sample}/R{mate}'
            out=summaries.setdefault(key,{'sample':sample,'mate':mate,'reads':0,'distinct_sequences':0,'states':Counter(),
                                         'labels':{},'complexity_by_state':{},'categories_by_mode':{m:Counter() for m in MODES},
                                         'reported_match_reads':Counter(),'accepted_match_reads':Counter(),'origin_screen_matches':Counter(),
                                         'origin_screen_matches_neither':Counter(),'local_clipped_bases':Counter(),'top_sequences':[]})
            ev=json.loads(data);state=ev['mode_state'];out['reads']+=n;out['distinct_sequences']+=1;out['states'][state]+=n
            for label in ev['labels']:
                tally=out['labels'].setdefault(label,{'reads':0,'distinct_sequences':0});tally['reads']+=n;tally['distinct_sequences']+=1
            c=out['complexity_by_state'].setdefault(state,{'reads':0,'distinct_sequences':0,'length_sum':0,'gc_sum':0,'entropy_sum':0,'mean_phred_sum':0,'labels':Counter(),'length_histogram':Counter(),'entropy_histogram':Counter(),'longest_run_histogram':Counter()})
            c['reads']+=n;c['distinct_sequences']+=1;c['length_sum']+=len(sequence)*n;c['gc_sum']+=ev['complexity']['gc_fraction']*n;c['entropy_sum']+=ev['complexity']['entropy_bits']*n;c['mean_phred_sum']+=qsum
            for label in ev['labels']:c['labels'][label]+=n
            c['length_histogram'][str(len(sequence))]+=n
            c['entropy_histogram'][str(round(ev['complexity']['entropy_bits'],1))]+=n
            c['longest_run_histogram'][str(max(ev['complexity']['longest_runs'].values()))]+=n
            for m,v in ev['alignment'].items():
                if v['reported_best_loci']:out['reported_match_reads'][m]+=n
                if v['accepted_best_loci']:
                    out['accepted_match_reads'][m]+=n;out['categories_by_mode'][m][v['category']]+=n
                    # Parent-specific counts are aggregated on disk below.
            example=ev['alignment']['local']['example_alignment']
            if example:
                # One representative accepted best alignment per read; ties are not extra reads.
                out['local_clipped_bases'][str(example['five_prime_clipped_bases']+example['three_prime_clipped_bases'])]+=n
            for origin in ev['origin_screens']:
                out['origin_screen_matches'][origin['label']]+=n
                if state=='neither_accepted':out['origin_screen_matches_neither'][origin['label']]+=n
            item={'sequence':sequence,'reads':n,'mode_state':state,'labels':ev['labels']}
            top=out['top_sequences'];top.append(item);top.sort(key=lambda x:(-x['reads'],x['sequence']));del top[20:]
            writer.writerow([sample,f'R{mate}',qid,sequence,n,qsum/n,data])
    # Publish detailed tables only after the final input-integrity checkpoint.
    for row in rows:
        sample=row['sample_id'];pair_summary[sample]={}
        for mode in MODES:
            counter=Counter()
            query='''SELECT a.ok,a.parent,b.ok,b.parent FROM pairs p
                JOIN accepted a ON a.qid=p.r1 AND a.mode=? JOIN accepted b ON b.qid=p.r2 AND b.mode=? WHERE p.sample=?'''
            for a,pa,b,pb in db.execute(query,(mode,mode,sample)):
                state=('same_parent' if pa and pb and pa==pb else 'different_parents' if pa and pb else 'parent_ambiguous') if a and b else 'one_mate_accepted' if a or b else 'neither_mate_accepted'
                counter[state]+=1
            expected=done(db,'input:'+sample)['analyzed_pairs']
            if sum(counter.values())!=expected:raise ValueError('Pair accounting failed')
            pair_summary[sample][mode]=dict(counter)
    parentfile=work/'parent_counts.tsv'
    with Path(str(parentfile)+'.partial').open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t');w.writerow(['sample','mate','mode','parent','reads','rpm_per_analyzed_mate_reads'])
        query="""SELECT c.sample,c.mate,a.mode,a.parent,sum(c.n) FROM counts c
            JOIN accepted a ON a.qid=c.qid WHERE a.ok=1 AND a.parent IS NOT NULL
            GROUP BY c.sample,c.mate,a.mode,a.parent ORDER BY c.sample,c.mate,a.mode,a.parent"""
        for sample,mate,mode,parent,n in db.execute(query):
            w.writerow([sample,f'R{mate}',mode,parent,n,1e6*n/summaries[f'{sample}/R{mate}']['reads']])
    for out in summaries.values():
        out['top_parents_by_mode']={}
        for mode in MODES:
            query='''SELECT a.parent,sum(c.n) AS n FROM counts c JOIN accepted a ON a.qid=c.qid
                WHERE c.sample=? AND c.mate=? AND a.mode=? AND a.ok=1 AND a.parent IS NOT NULL
                GROUP BY a.parent ORDER BY n DESC,a.parent LIMIT 20'''
            out['top_parents_by_mode'][mode]=[{'parent':p,'reads':n,'rpm':1e6*n/out['reads']}
                for p,n in db.execute(query,(out['sample'],out['mate'],mode))]
        assert sum(out['states'].values())==out['reads']
        for mode in MODES:assert sum(out['categories_by_mode'][mode].values())==out['accepted_match_reads'][mode]
    return {'samples':summaries,'pair_agreement':pair_summary,'sequence_evidence_path':str(details),
            'parent_counts_path':str(parentfile),'summary_seconds':time.monotonic()-started}


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-report',required=True,help='Completed main-pipeline JSON; never modified')
    p.add_argument('--manifest',required=True,help='Paired manifest; only condition RBC is analyzed')
    p.add_argument('--reference',required=True,help='Same combined reference as main analysis')
    p.add_argument('--output-dir',required=True,help='Dedicated resumable workspace and output directory')
    p.add_argument('--index-cache',required=True,help='Reusable reference index cache')
    p.add_argument('--pilot-pairs',type=int,default=10000,help='Uniform sample per RBC library; 0 means all pairs')
    p.add_argument('--sampling-seed',type=int,default=20261007)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--min-aligned-length',type=int,default=30)
    p.add_argument('--min-aligned-fraction',type=float,default=.5)
    p.add_argument('--max-edit-rate',type=float,default=.1)
    p.add_argument('--homopolymer-length',type=int,default=10)
    p.add_argument('--entropy-threshold',type=float,default=1.2)
    p.add_argument('--kmer-threshold',type=float,default=.1)
    p.add_argument('--screen-reference',action='append',default=[],help='genome=FASTA, microbial=FASTA, technical=FASTA, or other=FASTA; labels remain nonexclusive evidence')
    p.add_argument('--adapter-sequence',action='append',default=[],help='Exact substring diagnostic, >=8 bases; never trims')
    return p


def run(args):
    if args.pilot_pairs<0 or min(args.threads,args.min_aligned_length,args.homopolymer_length)<=0:
        raise ValueError('Invalid sampling/thread/length settings')
    if not 0<args.min_aligned_fraction<=1 or not 0<=args.max_edit_rate<=1 or not 0<=args.entropy_threshold<=2 or not 0<=args.kmer_threshold<=1:
        raise ValueError('Invalid complexity/alignment thresholds')
    for s in args.adapter_sequence:
        if len(clean_sequence(s))<8:raise ValueError('Adapter diagnostics require at least 8 bases')
    report_path=Path(args.analysis_report).expanduser().resolve(strict=True)
    report=json.loads(report_path.read_text())
    if report.get('status')!='complete' or not str(report.get('schema_version','')).startswith('2.'):
        raise ValueError('Expected a completed main-pipeline schema 2 report')
    rows=selected_manifest(args.manifest,report)
    refs={'primary':str(Path(args.reference).expanduser().resolve(strict=True))}
    for spec in args.screen_reference:
        label,sep,path=spec.partition('=')
        if not sep or label not in ('genome','microbial','technical','other') or label in refs:
            raise ValueError('Screen labels must be unique genome, microbial, technical, or other')
        refs[label]=str(Path(path).expanduser().resolve(strict=True))
    tools={}
    for name in ('bowtie2','bowtie2-build'):
        binary=shutil.which(name)
        if not binary:raise ValueError(f'Required executable not on PATH: {name}')
        version=subprocess.check_output([binary,'--version'],text=True).splitlines()[0]
        tools[name]={'path':binary,'version':version}
    work=Path(args.output_dir).expanduser().resolve()
    with lock(work):
        started=time.monotonic()
        if not (work/'run_config.json').exists() and any(p.name!='.lock' for p in work.iterdir()):
            raise ValueError('Choose an empty dedicated output directory; refusing to overwrite unrelated files')
        inputs={str(report_path),str(Path(args.manifest).resolve()),*refs.values(),*(row[m] for row in rows for m in ('r1','r2'))}
        reserved = {'diagnostics.sqlite','queries.fa','run_config.json','rbc_diagnostics.json','sequence_evidence.tsv.gz','parent_counts.tsv'}
        reserved |= {f'{row["sample_id"]}.pilot.R{m}.fastq.gz' for row in rows for m in (1,2)}
        if any((work/name).resolve() in {Path(p).resolve() for p in inputs} for name in reserved):
            raise ValueError('Output must not overwrite an input')
        hashes={p:sha256(p) for p in sorted(inputs)}
        for row in rows:
            original=report.get('paired_input_qc',{}).get(row['sample_id'],{})
            for mate in ('r1','r2'):
                digest=report.get('input_provenance',{}).get(original.get(mate,''))
                if digest and hashes[row[mate]]!=digest:
                    raise ValueError(f'{row["sample_id"]}/{mate}: checksum differs from main-analysis input')
        if hashes[refs['primary']]!=report['reference']['sha256']:
            raise ValueError('Primary reference checksum differs from main analysis')
        config={'version':VERSION,'parameters':vars(args),'sha256':hashes,'tools':tools,
                'code_sha256':{'csrna_rbc_diagnostics.py':sha256(Path(__file__)),
                               'csrna_pipeline.py':sha256(Path(sys.modules['csrna_pipeline'].__file__))}}
        config_path=work/'run_config.json'
        if config_path.exists() and json.loads(config_path.read_text())!=config:
            raise ValueError('Workspace configuration or input checksums changed. Choose a new output directory.')
        atomic_json(config_path,config)
        db=connect(work/'diagnostics.sqlite')
        try:
            audits={row['sample_id']:ingest(db,row,args,work) for row in rows}
            for sample,a in audits.items():
                original=report.get('paired_input_qc',{}).get(sample,{})
                if original.get('validated_pairs') is not None and original['validated_pairs']!=a['validated_pairs']:
                    raise ValueError(f'{sample}: pair count differs from main report')
            ingestion_seconds=time.monotonic()-started
            query=work/'queries.fa'
            with query.open('w') as h:
                for qid,seq in db.execute('SELECT id,sequence FROM seq WHERE id IN (SELECT qid FROM counts) ORDER BY id'):
                    h.write(f'>q{qid}\n{seq}\n')
            stages={}
            for label,path in refs.items():
                start=time.monotonic();prefix=index_reference(db,label,path,args,tools,work)
                stages[f'index:{label}']={'seconds':time.monotonic()-start}
                for mode in MODES if label=='primary' else ('local',):
                    print(f'ALIGN {label} {mode}',file=sys.stderr,flush=True)
                    stages[f'{label}:{mode}']=align(db,label,mode,prefix,args,tools,work,query)
            result=summarize(db,args,rows,work)
            if any(sha256(p)!=digest for p,digest in hashes.items()):raise ValueError('Input changed during analysis; no new report published')
            result.update({'schema_version':SCHEMA,'status':'complete','version':VERSION,
                'created_utc':datetime.now(timezone.utc).isoformat(),'analysis_report_sha256':hashes[str(report_path)],
                'parameters':vars(args),'input_scope':'full_library' if not args.pilot_pairs else 'selected_subset',
                'sample_metadata':[{**next(m for m in report['sample_metadata'] if m['sample_id']==row['sample_id']), **{k:v for k,v in row.items() if k not in ('r1','r2')}} for row in rows],
                'main_analysis_sample_metadata':[m for m in report['sample_metadata'] if m['sample_id'] in audits],
                'metadata_policy':'Current manifest fields take precedence for diagnostic annotations only; original main-report metadata is retained separately. Counts and identity checks are unchanged.',
                'input_audit':audits,'provenance':config,'stages':stages,'input_stage_seconds_this_invocation':ingestion_seconds,
                'screen_status':{label:'tested' if label in refs else 'not_tested' for label in ('genome','microbial','technical','other')},
                'interpretation':['RBC-associated preparation; surface localization is not established.',
                'Mates aligned independently, not native paired alignment. Pair agreement counts pairs once.',
                'Library fragmentation and cDNA synthesis can create boundaries; even reconstructed library inserts would not establish native RNA ends. Fragmentation status remains as supplied.',
                'Strand-specific preparation does not by itself verify mate-to-RNA orientation. UDI indexes are not evidence of UMIs. Conflicting PCR records require review, not numerical count correction.',
                'FASTA sequence-level alignment ignores base qualities for scoring; mean quality is diagnostic only.',
                'Local boundaries are alignment boundaries, not native fragment ends. Raw scores across modes are not comparable.',
                'Best reported loci retained; no equivalence projection or arbitrary hit cap. Ambiguous parents are excluded from parent-specific counts.',
                'Coverage/edit thresholds apply after best-score selection. A rejected best match is not no sequence similarity.',
                'Origin screens are overlapping matches, not contamination verdicts. Genome matching is unspliced evidence only.',
                'No trimming, UMI inference, deduplication of abundance, or automatic sample exclusion.',
                'Full-library capacity depends on disk and all-hit alignment volume; test resources before full runs.'],
                'cross_cell_comparison':'Use consistently processed data; do not combine these raw-reference mate counts with the main endpoint feature counts.'})
            for key in ('sequence_evidence_path','parent_counts_path'):
                target=Path(result[key]);os.replace(str(target)+'.partial',target)
            atomic_json(work/'rbc_diagnostics.json',result)
            print(f'Complete: {work / "rbc_diagnostics.json"}')
            return result
        finally: db.close()


if __name__=='__main__':
    try: run(parser().parse_args())
    except (ValueError,OSError,RuntimeError,sqlite3.Error,KeyError) as error:
        print(f'RBC diagnostics failed: {error}',file=sys.stderr)
        sys.exit(1)
