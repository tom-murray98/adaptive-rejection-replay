#!/usr/bin/env python3
"""Paired, fixed-time adaptive-sampling replay; inputs are never modified."""
import argparse
import collections
import csv
import gzip
import json
import math
import sqlite3
import sys
import time
from pathlib import Path


class Progress:
    """Time-throttled stderr updates, visible even when stdout is redirected."""
    def __init__(self, label, total=None):
        self.label, self.total = label, total
        self.started = self.last = time.monotonic()
        self.update(0, force=True)

    def update(self, count, force=False):
        now = time.monotonic()
        if not force and now - self.last < 5:
            return
        elapsed = now - self.started
        detail = f'{count:,}'
        if self.total:
            detail += f'/{self.total:,} ({100 * count / self.total:.1f}%)'
        print(f'[{time.strftime("%H:%M:%S")}] {self.label}: {detail}; '
              f'elapsed {elapsed:.0f}s', file=sys.stderr, flush=True)
        self.last = now


def op(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def replay(events, window, fraction, overhead):
    """Events are (id, start, duration, length, rejected), ordered per channel.

    Fixed overhead is the part of a rejected read's duration that is not sped up.
    Inter-event gaps remain fixed. Only completed molecules enter the BAM.
    """
    previous_end = 0.0
    saved = 0.0
    for rid, start, duration, length, rejected in events:
        if start < previous_end - 0.01:
            raise ValueError(f'Overlapping summary events at {rid}; child-read/channel reconciliation needed')
        if duration <= 0 or length < 0 or not all(math.isfinite(x) for x in (start, duration)):
            raise ValueError(f'Invalid timing/length at {rid}')
        original_end = start + duration
        new_duration = duration
        if rejected:
            fixed = min(overhead, duration)
            new_duration = fixed + fraction * (duration - fixed)
        fast_end = start - saved + new_duration
        baseline = original_end <= window
        faster = fast_end <= window
        yield rid, baseline, faster, original_end, fast_end
        saved += duration - new_duration
        previous_end = original_end


def prepare(a):
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    dbpath = out / 'selection.sqlite'
    if dbpath.exists():
        raise ValueError('Selection database exists; choose a new output directory')
    db = sqlite3.connect(dbpath)
    db.executescript('''CREATE TABLE events(id TEXT PRIMARY KEY, channel INTEGER, start REAL,
        duration REAL, length INTEGER, rejected INTEGER);
        CREATE TABLE selection(id TEXT PRIMARY KEY, baseline INTEGER, faster INTEGER, rejected INTEGER);
        CREATE TABLE matched(id TEXT PRIMARY KEY);''')
    ids = None
    if a.ids:
        progress = Progress('Loading unblock IDs')
        ids = set()
        count = 0
        with op(a.ids) as f:
            for count, line in enumerate(f, 1):
                if line.strip():
                    ids.add(line.strip())
                if count % 20000 == 0:
                    progress.update(count)
        progress.update(count, force=True)
        if not ids:
            raise ValueError('Empty unblock-ID list')
    reasons = collections.Counter()
    batch = []
    progress = Progress('Reading summary events')
    count = 0
    with op(a.summary) as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'read_id', 'channel', 'start_time', 'duration', 'sequence_length_template', 'end_reason'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Missing summary fields: {required - set(reader.fieldnames or [])}')
        for r in reader:
            count += 1
            confirmed = r['end_reason'] == a.reject_reason
            listed = ids is None or r['read_id'] in ids or r.get('parent_read_id') in ids
            rejected = int(confirmed and listed)
            reasons[r['end_reason']] += 1
            batch.append((r['read_id'], int(r['channel']), float(r['start_time']),
                          float(r['duration']), int(r['sequence_length_template']), rejected))
            if len(batch) >= 20000:
                db.executemany('INSERT INTO events VALUES (?,?,?,?,?,?)', batch)
                db.commit(); batch.clear()
                progress.update(count)
        db.executemany('INSERT INTO events VALUES (?,?,?,?,?,?)', batch)
    progress.update(count, force=True)
    print('Building channel/time index; this may take a while...', file=sys.stderr, flush=True)
    # SQLite invokes this during long SQL operations, including index creation.
    progress = Progress('Building index (SQLite steps)')
    steps = 0
    def sql_progress():
        nonlocal steps
        steps += 100000
        progress.update(steps)
        return 0
    db.set_progress_handler(sql_progress, 100000)
    db.execute('CREATE INDEX event_time ON events(channel,start)')
    db.set_progress_handler(None, 0)
    print('Channel/time index ready.', file=sys.stderr, flush=True)
    if not db.execute('SELECT count(*) FROM events WHERE rejected=1').fetchone()[0]:
        raise ValueError('No confirmed rejection events found')
    stats = {'parameters': vars(a), 'end_reasons': dict(reasons), 'channels': 0,
             'baseline_events': 0, 'faster_events': 0, 'reserve_exhausted_channels': []}
    channels = db.execute('SELECT DISTINCT channel FROM events ORDER BY channel').fetchall()
    progress = Progress('Planning channels', len(channels))
    for (channel,) in channels:
        stats['channels'] += 1
        flags = []
        events = db.execute('SELECT id,start,duration,length,rejected FROM events WHERE channel=? ORDER BY start', (channel,))
        last = None
        for item in replay(events, a.seconds, a.fraction, a.fixed_overhead):
            rid, base, fast, original_end, fast_end = item
            last = item
            if base or fast:
                flags.append((rid, int(base), int(fast), db.execute('SELECT rejected FROM events WHERE id=?',(rid,)).fetchone()[0]))
                stats['baseline_events'] += int(base)
                stats['faster_events'] += int(fast)
            if len(flags) >= 10000:
                db.executemany('INSERT INTO selection VALUES (?,?,?,?)', flags); flags.clear()
            if fast_end > a.seconds:
                break
        db.executemany('INSERT INTO selection VALUES (?,?,?,?)', flags)
        if last and last[4] <= a.seconds:
            stats['reserve_exhausted_channels'].append(channel)
        db.commit()
        progress.update(stats['channels'])
    progress.update(stats['channels'], force=True)
    stats['reserve_warning'] = 'Exhaustion may indicate a naturally inactive pore or insufficient held-out data; inspect before interpretation.'
    (out / 'plan.json').write_text(json.dumps(stats, indent=2) + '\n')
    db.close()
    print(f'Plan saved: {out / "plan.json"}', file=sys.stderr, flush=True)
    print(json.dumps(stats, indent=2))


def export(a):
    import pysam
    out = Path(a.outdir)
    plan = json.loads((out / 'plan.json').read_text())
    db = sqlite3.connect(out / 'selection.sqlite')
    if db.execute('SELECT count(*) FROM matched').fetchone()[0]:
        raise ValueError('Export already started; use a fresh plan directory')
    fraction = plan['parameters']['fraction']
    stats = {c: collections.Counter() for c in ('baseline', 'trim_only', 'faster')}
    progress = Progress('Scanning BAM records')
    scanned = 0
    with pysam.AlignmentFile(a.bam, 'rb', check_sq=False) as bam:
        rg = bam.header.to_dict().get('RG', [])
        if len(rg) != 1:
            raise ValueError('Exactly one BAM read group is required')
        # Use FASTQ + original subset per condition, then realign and restore tags.
        writers = {}
        for c in stats:
            d = out / c; d.mkdir()
            (d / 'read_group.txt').write_text('\\t'.join(['@RG']+[f'{k}:{v}' for k,v in rg[0].items()])+'\n')
            writers[c] = (gzip.open(d/'reads.fastq.gz','wt',compresslevel=1),
                          pysam.AlignmentFile(d/'original.bam','wb',template=bam))
        try:
            for r in bam.fetch(until_eof=True):
                scanned += 1
                if scanned % 20000 == 0:
                    progress.update(scanned)
                if r.is_secondary or r.is_supplementary:
                    continue
                flag = db.execute('SELECT baseline,faster,rejected FROM selection WHERE id=?',(r.query_name,)).fetchone()
                if flag is None:
                    continue
                db.execute('INSERT INTO matched VALUES (?)',(r.query_name,))
                seq, qual = r.get_forward_sequence(), r.get_forward_qualities()
                if not seq or qual is None or any(op == 5 for op,n in (r.cigartuples or [])):
                    raise ValueError(f'Missing sequence/quality or hard-clipped primary: {r.query_name}')
                if r.has_tag('MM') != r.has_tag('ML') or not r.has_tag('MM'):
                    raise ValueError(f'Missing MM/ML methylation tags: {r.query_name}')
                for c, selected in [('baseline',flag[0]),('trim_only',flag[0]),('faster',flag[1])]:
                    if not selected: continue
                    n = max(1,int(len(seq)*fraction)) if flag[2] and c!='baseline' else len(seq)
                    fq, original = writers[c]
                    fq.write('@'+r.query_name+'\n'+seq[:n]+'\n+\n'+''.join(chr(q+33) for q in qual[:n])+'\n')
                    original.write(r)
                    stats[c]['reads'] += 1; stats[c]['bases'] += n
                    stats[c]['rejected_reads'] += flag[2]
        finally:
            for fq,b in writers.values(): fq.close(); b.close()
    db.commit()
    progress.update(scanned, force=True)
    for c in stats:
        if not stats[c]['reads']: raise ValueError(f'Empty condition: {c}')
    missing = db.execute('SELECT count(*) FROM selection LEFT JOIN matched USING(id) WHERE matched.id IS NULL').fetchone()[0]
    result = {'conditions':stats, 'selected_summary_events_missing_from_bam':missing,
              'note':'All summary events consume time; only events present in the input BAM supply classifier observations.'}
    (out/'export.json').write_text(json.dumps(result,indent=2)+'\n')
    db.close()
    print(json.dumps(result,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare')
    q.add_argument('--summary',required=True); q.add_argument('--ids')
    q.add_argument('--seconds',type=float,required=True)
    q.add_argument('--fraction',type=float,default=0.5)
    q.add_argument('--fixed-overhead',type=float,default=0)
    q.add_argument('--reject-reason',default='data_service_unblock_mux_change')
    q.add_argument('--outdir',required=True)
    q=sub.add_parser('export');q.add_argument('--bam',required=True);q.add_argument('--outdir',required=True)
    a=p.parse_args()
    if a.command=='prepare':
        if a.seconds<=0 or a.fraction!=0.5 or a.fixed_overhead<0:
            p.error('Require seconds>0, fraction=0.5 (tag restoration supports half-length), overhead>=0')
        prepare(a)
    else: export(a)

if __name__=='__main__': main()
