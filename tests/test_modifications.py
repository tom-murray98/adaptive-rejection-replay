"""Regression tests for sequence-relative methylation tag restoration."""
from array import array
from collections import Counter
import unittest

import pysam

from asreplay.modifications import clip_modifications, restore_group


class RestoreTests(unittest.TestCase):
    sequence = 'ACCGTCCGCC'
    mm = 'C+mh?,0,1,2;C+h.;'
    ml = array('B', [10, 20, 30, 40, 50, 60])

    def record(self, seq, reverse=False, cigar=None):
        header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 100}]})
        r = pysam.AlignedSegment(header)
        r.query_name = 'read'
        r.flag = 16 if reverse else 0
        r.reference_id = 0
        r.reference_start = 0
        r.query_sequence = seq
        r.query_qualities = [30] * len(seq)
        r.cigarstring = cigar or f'{len(seq)}M'
        return r

    def source(self):
        r = self.record('GGCGGACGGT', True)
        r.set_tag('MM', self.mm)
        r.set_tag('ML', self.ml)
        r.set_tag('MN', 10)
        return r

    def test_interleaved_probabilities_and_clipping(self):
        mm, ml = clip_modifications(self.sequence, self.mm, self.ml, 2, 8)
        self.assertEqual(mm, 'C+mh?,1;C+h.;')
        self.assertEqual(list(ml), [30, 40])

    def test_orientation_and_supplementary_hard_clipping(self):
        primary = self.record('ACCGT')
        reverse = self.record('ACGGT', True)
        reverse.flag |= 2048
        clipped = self.record('CCGT', cigar='1H4M')
        clipped.flag |= 2048
        reverse_clipped = self.record('ACGG', True, '4M1H')
        reverse_clipped.flag |= 2048
        restore_group(self.source(), [primary, reverse, clipped, reverse_clipped], Counter())
        self.assertEqual(primary.modified_bases_forward, reverse.modified_bases_forward)
        self.assertEqual(clipped.modified_bases_forward, reverse_clipped.modified_bases_forward)
        self.assertEqual(primary.get_tag('MM'), 'C+mh?,0;C+h.;')
        self.assertEqual(list(primary.get_tag('ML')), [10, 20])
        self.assertEqual(clipped.get_tag('MN'), 4)

    def test_empty_calls_one_base_and_missing_secondary_sequence(self):
        primary = self.record('A')
        secondary = pysam.AlignedSegment(primary.header)
        secondary.query_name = 'read'
        secondary.flag = 256
        source = self.record('AC')
        source.set_tag('MM', 'C+m?,0;')
        source.set_tag('ML', array('B', [200]))
        restore_group(source, [primary, secondary], Counter())
        self.assertEqual(primary.get_tag('MM'), 'C+m?;')
        self.assertEqual(list(primary.get_tag('ML')), [])
        self.assertEqual(primary.get_tag('MN'), 1)
        self.assertFalse(secondary.has_tag('MM'))

    def test_rejects_mismatched_sequence_and_stale_tags(self):
        with self.assertRaisesRegex(ValueError, 'sequence mismatch'):
            restore_group(self.source(), [self.record('AAAAA')], Counter())
        stale = self.source()
        stale.set_tag('MN', 11)
        with self.assertRaisesRegex(ValueError, 'Stale original MN'):
            restore_group(stale, [self.record('ACCGT')], Counter())


if __name__ == '__main__':
    unittest.main()
