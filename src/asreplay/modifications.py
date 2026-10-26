#!/usr/bin/env python3
"""Restore clipped MM/ML/MN from original reads to name-sorted realignments.

Both BAMs must have the same read-name ordering (samtools sort -n). The
shortened BAM must retain either the full read or its first floor(L/2) bases.
Modification tags always refer to the as-sequenced orientation (SAMtags spec).
"""
import argparse
from array import array
from collections import Counter
from itertools import groupby, zip_longest
import json
import re

import pysam


GROUP = re.compile(r'([ACGTUN])([+-])([a-z]+|[0-9]+|[ACGTUN])([.?]?)(.*)')


def clip_modifications(sequence, mm, ml, start, end):
    """Re-encode MM deltas and matching interleaved ML for [start, end)."""
    if not 0 <= start <= end <= len(sequence):
        raise ValueError('Invalid modification clipping interval')
    if mm and not mm.endswith(';'):
        raise ValueError('Malformed MM tag: missing trailing semicolon')
    encoded, probabilities = [], array('B')
    offset = 0
    for group in mm.split(';')[:-1]:
        match = GROUP.fullmatch(group)
        if match is None:
            raise ValueError(f'Unsupported MM group: {group}')
        base, strand, codes, mode, tail = match.groups()
        if tail and not re.fullmatch(r'(,[0-9]+)+', tail):
            raise ValueError(f'Malformed MM deltas: {group}')
        deltas = [int(x) for x in tail.split(',')[1:]] if tail else []
        width = 1 if codes.isdigit() else len(codes)
        positions = [i for i, b in enumerate(sequence.upper()) if base == 'N' or b == base]
        before = sum(i < start for i in positions)
        previous, index = -1, -1
        kept = []
        for delta in deltas:
            index += delta + 1
            if index >= len(positions):
                raise ValueError('MM position extends beyond original sequence')
            values = ml[offset:offset + width]
            if len(values) != width:
                raise ValueError('ML tag is shorter than MM requires')
            offset += width
            if start <= positions[index] < end:
                local = index - before
                kept.append(local - previous - 1)
                previous = local
                probabilities.extend(values)
        encoded.append(base + strand + codes + mode + ''.join(f',{d}' for d in kept) + ';')
    if offset != len(ml):
        raise ValueError('ML tag has extra probabilities')
    return ''.join(encoded), probabilities


def restore_group(original, records, stats):
    seq = original.get_forward_sequence()
    if not seq or any(op == 5 for op, _ in (original.cigartuples or [])):
        raise ValueError(f'Missing or hard-clipped original sequence: {original.query_name}')
    primary = [r for r in records if not r.is_secondary and not r.is_supplementary]
    if len(primary) != 1:
        raise ValueError(f'Expected one shortened primary: {original.query_name}')
    length = primary[0].query_length
    if length not in (len(seq), max(1, len(seq) // 2)):
        raise ValueError(f'Unexpected shortened read length: {original.query_name}')
    stats['primary_reads'] += 1
    mm_tag = 'MM' if original.has_tag('MM') else 'Mm'
    ml_tag = 'ML' if original.has_tag('ML') else 'Ml'
    has_mm, has_ml = original.has_tag(mm_tag), original.has_tag(ml_tag)
    if has_mm != has_ml:
        raise ValueError(f'Incomplete MM/ML pair: {original.query_name}')
    if has_mm:
        if original.has_tag('MN') and original.get_tag('MN') != len(seq):
            raise ValueError(f'Stale original MN tag: {original.query_name}')
        mm, ml = clip_modifications(seq, original.get_tag(mm_tag), original.get_tag(ml_tag), 0, length)
        stats['original_primary_reads_with_modifications'] += 1
    for r in records:
        for tag in ('MM', 'ML', 'MN', 'Mm', 'Ml'):
            r.set_tag(tag, None)
        stats['alignment_records'] += 1
        if r.query_sequence is None:
            stats['records_without_sequence'] += 1
            continue
        cigar = r.cigartuples or []
        left = cigar[0][1] if cigar and cigar[0][0] == 5 else 0
        right = cigar[-1][1] if cigar and cigar[-1][0] == 5 else 0
        if left + r.query_length + right != length:
            raise ValueError(f'Inconsistent hard clipping: {r.query_name}')
        start = right if r.is_reverse else left
        end = start + r.query_length
        if r.get_forward_sequence() != seq[start:end]:
            raise ValueError(f'Original/shortened sequence mismatch: {r.query_name}')
        if has_mm:
            clipped_mm, clipped_ml = clip_modifications(seq[:length], mm, ml, start, end)
            r.set_tag('MM', clipped_mm, 'Z')
            r.set_tag('ML', clipped_ml)
            r.set_tag('MN', r.query_length, 'i')
            # Independently decode with HTSlib to catch invalid encodings.
            if r.modified_bases_forward is None:
                raise ValueError(f'HTSlib could not decode modifications: {r.query_name}')
            stats['records_with_restored_tags'] += 1
            stats['retained_probability_values'] += len(clipped_ml)
            if not r.is_secondary and not r.is_supplementary:
                stats['primary_reads_with_restored_tags'] += 1
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', required=True)
    p.add_argument('--shortened', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--stats', required=True)
    a = p.parse_args()
    stats = Counter()
    with pysam.AlignmentFile(a.original, 'rb') as source, pysam.AlignmentFile(a.shortened, 'rb') as target:
        original_groups = groupby((r for r in source if not r.is_secondary and not r.is_supplementary), lambda r: r.query_name)
        target_groups = groupby(target, lambda r: r.query_name)
        with pysam.AlignmentFile(a.output, 'wb', template=target) as output:
            for src, dst in zip_longest(original_groups, target_groups):
                if src is None or dst is None or src[0] != dst[0]:
                    raise ValueError('Read names differ between name-sorted original and shortened BAMs')
                originals = list(src[1])
                if len(originals) != 1:
                    raise ValueError(f'Duplicate original primary: {src[0]}')
                for record in restore_group(originals[0], list(dst[1]), stats):
                    output.write(record)
                if stats['primary_reads'] % 100000 == 0:
                    print(f"Processed {stats['primary_reads']:,} reads", flush=True)
    if not stats['primary_reads_with_restored_tags']:
        raise ValueError('No primary reads received modification tags')
    with open(a.stats, 'w') as f:
        json.dump(dict(stats), f, indent=2)
        f.write('\n')
    print(json.dumps(dict(stats), indent=2))


if __name__ == '__main__':
    main()
