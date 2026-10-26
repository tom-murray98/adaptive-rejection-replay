"""Real BAM/aligner orchestration with stub classifiers; no models or Docker needed."""
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
from array import array
import pysam
from asreplay.runner import run


@unittest.skipUnless(os.environ.get('ASREPLAY_INTEGRATION_TOOLCHAIN')=='1', 'Set ASREPLAY_INTEGRATION_TOOLCHAIN=1 with samtools/minimap2 on PATH')
class PipelineTests(unittest.TestCase):
    def test_complete_run_and_cleanup_with_real_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            reference=root/'ref.fa'; reference.write_text('>chr1\n'+('ACGT'*1000)+'\n')
            pysam.faidx(str(reference))
            bam=root/'original.bam'
            header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':4000}],'RG':[{'ID':'test','SM':'test'}]}
            with pysam.AlignmentFile(bam,'wb',header=header) as output:
                for i in range(3):
                    r=pysam.AlignedSegment(output.header)
                    r.query_name=f'r{i}';r.query_sequence='ACGT'*50;r.query_qualities=array('B',[30]*200)
                    r.flag=4;r.set_tag('RG','test');r.set_tag('MM','C+m?,0;');r.set_tag('ML',array('B',[200]))
                    output.write(r)
            summary=root/'summary.tsv'
            summary.write_text('read_id\tchannel\tstart_time\tduration\tsequence_length_template\tend_reason\n'
              'r0\t1\t0\t4\t200\tdata_service_unblock_mux_change\n'
              'r1\t1\t5\t4\t200\tsignal_positive\n'
              'r2\t1\t10\t4\t200\tsignal_positive\n')
            ids=root/'ids.txt';ids.write_text('r0\n')
            workflow=root/'workflow';workflow.mkdir()
            for name in ('bin','wf-human-variation'):(workflow/name).mkdir()
            (workflow/'main.nf').write_text('// synthetic classifier orchestration stub\n')
            (workflow/'rapidcns2.R').write_text('meth <- read.delim(opt$in_file,header=FALSE)\n')
            tools=root/'tools';tools.mkdir()
            docker=tools/'docker';docker.write_text('#!/bin/sh\necho test-container\n');docker.chmod(0o755)
            nf=tools/'nextflow';nf.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
if '--version' in sys.argv:
 print('synthetic-nextflow-stub');sys.exit(0)
a=sys.argv
work=Path(a[a.index('-work-dir')+1])/'aa'/'task';work.mkdir(parents=True)
(work/'.command.sh').write_text('synthetic classifier fixture')
(work/'test_nanoDx_votes.txt').write_text('class\\tscore\\tnum_features\\nGBM RTK II\\t0.5\\t10\\n')
(work/'merged_CAPPER_MODEL.csv').write_text('number_probes,GBM RTK II\\n12,0.8\\n')
(work/'test_calibrated_classification.tsv').write_text('Number of features: 5\\nPredicted Class: GBM_RTK_II\\nInitial Score: 0.06\\nCalibrated Score: 0.9\\n')
''');nf.chmod(0o755)
            out=root/'output'
            config=root/'config.json'
            config.write_text(json.dumps(dict(sample='test',bam=str(bam),summary=str(summary),ids=str(ids),reference=str(reference),workflow=str(workflow),outdir=str(out),seconds=7,threads=4,memory_gb=8)))
            with patch.dict(os.environ,{'PATH':str(tools)+os.pathsep+os.environ['PATH']}):
                run(config)
            self.assertTrue((out/'EXPERIMENT_SUCCESS').exists())
            yields=json.loads((out/'export.json').read_text())['conditions']
            self.assertEqual(yields['baseline']['reads'],1)
            self.assertEqual(yields['trim_only']['reads'],1)
            self.assertEqual(yields['faster']['reads'],2)
            for condition in ('baseline','trim_only','faster'):
                d=out/condition
                self.assertTrue((d/'input.bam.bai').is_file())
                self.assertFalse((d/'work').exists())
                self.assertFalse((d/'reads.fastq.gz').exists())
                self.assertTrue((d/'classifications/test_calibrated_classification.tsv').exists())
