import argparse
from array import array
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
import pysam
from asreplay.simulation import replay, prepare, export


class SimulationTests(unittest.TestCase):
    def test_saved_time_reaches_new_molecule(self):
        events=[('reject',0,4,100,1),('accept',5,4,100,0),('later',10,4,100,1)]
        result=list(replay(events,7,0.5,0))
        self.assertEqual([(r[1],r[2]) for r in result],[(True,True),(False,True),(False,False)])
        self.assertEqual(result[1][4],7)  # 1-second capture gap is retained
        self.assertEqual(result[2][4],10)  # additional rejected reads save time too

    def test_overhead_not_halved(self):
        result=list(replay([('r',0,4,100,1)],4,0.5,2))
        self.assertEqual(result[0][4],3)

    def test_overlap_fails(self):
        with self.assertRaises(ValueError):
            list(replay([('a',0,4,100,1),('b',3,4,100,0)],10,0.5,0))

    def test_plan_export_paired_and_reverse(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);summary=d/'summary.tsv';out=d/'experiment'
            with summary.open('w') as f:
                w=csv.writer(f,delimiter='\t')
                w.writerow(['read_id','channel','start_time','duration','sequence_length_template','end_reason'])
                w.writerow(['r',1,0,4,10,'data_service_unblock_mux_change'])
                w.writerow(['a',1,5,4,10,'signal_positive'])
                w.writerow(['z',1,10,4,10,'signal_positive'])
            prepare(argparse.Namespace(outdir=str(out),summary=str(summary),ids=None,
                seconds=7,fraction=0.5,fixed_overhead=0,reject_reason='data_service_unblock_mux_change',command='prepare'))
            bam=d/'input.bam'
            header={'SQ':[{'SN':'chr1','LN':100}], 'RG':[{'ID':'test','PL':'ONT'}]}
            with pysam.AlignmentFile(bam,'wb',header=header) as f:
                for name in ['r','a','z']:
                    r=pysam.AlignedSegment(f.header);r.query_name=name;r.flag=16
                    r.reference_id=0;r.reference_start=0;r.cigarstring='10M'
                    r.query_sequence='GGCGGACGGT';r.query_qualities=[30]*10
                    r.set_tag('MM','C+m?,0;');r.set_tag('ML',array('B',[200]));f.write(r)
            export(argparse.Namespace(bam=str(bam),outdir=str(out)))
            stats=json.loads((out/'export.json').read_text())['conditions']
            self.assertEqual(stats['baseline']['reads'],1)
            self.assertEqual(stats['trim_only']['bases'],5)
            self.assertEqual(stats['faster']['reads'],2)
            with gzip.open(out/'trim_only/reads.fastq.gz','rt') as f:
                lines=f.read().splitlines()
            self.assertEqual(lines[1],'ACCGT')  # first bases in pore order, not BAM order
            with self.assertRaises(ValueError):
                export(argparse.Namespace(bam=str(bam),outdir=str(out)))

if __name__=='__main__':unittest.main()
