from pathlib import Path
import contextlib
import io
import json
import tempfile
import unittest
import pysam
from asreplay.cleanup import clean
from asreplay.cli import main, slurm_script
from asreplay.runner import load_config, report
from asreplay.rapidcns2 import repair


class PackageTests(unittest.TestCase):
    def test_modification_fix_idempotent_and_rejects_unknown_script(self):
        text = 'meth <- read.delim(opt$in_file,header=FALSE)\n'
        patched = repair(text)
        self.assertIn('as.character(meth[[4]]) == "m"', patched)
        self.assertEqual(repair(patched), patched)
        with self.assertRaises(ValueError):
            repair('unrecognised input parser')

    def test_cleanup_refuses_unowned_and_incomplete_and_preserves_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError):
                clean(root, 'baseline', 'build')
            (root/'.asreplay-owned.json').write_text(json.dumps({'outdir': str(root)}))
            d = root/'baseline'; d.mkdir()
            with pysam.AlignmentFile(d/'input.bam', 'wb', header={'SQ': [{'SN':'chr1','LN':1000}]}):
                pass
            pysam.index(str(d/'input.bam'))
            (d/'BUILD_SUCCESS').touch()
            (d/'original.bam').write_bytes(b'redundant')
            clean(root, 'baseline', 'build')
            self.assertFalse((d/'original.bam').exists())
            self.assertTrue((d/'input.bam').is_file())
            clean(root, 'baseline', 'build')
            with self.assertRaises(ValueError):
                clean(root, 'baseline', 'classification')
            work=d/'work'/'ab'/'task'; work.mkdir(parents=True)
            (work/'.command.sh').write_text('echo command\n')
            (work/'prediction.txt').write_text('preserved')
            results=d/'classifications'; results.mkdir()
            (results/'prediction.txt').symlink_to(work/'prediction.txt')
            for name in ('x_nanoDx_votes.txt','x_calibrated_classification.tsv','x_CAPPER_MODEL.csv'):
                (results/name).write_text('result')
            (d/'CLASSIFICATION_SUCCESS').touch()
            clean(root, 'baseline', 'classification')
            self.assertEqual((results/'prediction.txt').read_text(), 'preserved')
            self.assertFalse((d/'work').exists())
            self.assertTrue((d/'task_provenance.tar.gz').is_file())
            self.assertTrue((d/'input.bam.bai').is_file())

    def test_cleanup_never_follows_condition_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'out'; root.mkdir()
            (root/'.asreplay-owned.json').write_text(json.dumps({'outdir': str(root)}))
            outside=Path(tmp)/'outside'; outside.mkdir()
            (root/'baseline').symlink_to(outside, target_is_directory=True)
            with self.assertRaises(ValueError):
                clean(root,'baseline','build')
            self.assertTrue(outside.is_dir())

    def test_config_paths_relative_to_config_and_bad_values_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name in ('bam','summary','ids','reference'):(root/name).touch()
            wf=root/'workflow'; wf.mkdir()
            for name in ('main.nf','rapidcns2.R','bin','wf-human-variation'):(wf/name).touch()
            c=dict(sample='sample',bam='bam',summary='summary',ids='ids',reference='reference',workflow='workflow',outdir='output',seconds=3600)
            path=root/'config.json'; path.write_text(json.dumps(c))
            self.assertEqual(load_config(path)['bam'], str(root/'bam'))
            for change in [dict(seconds=float('nan')),dict(sample='a\n#SBATCH --bad'),dict(threads=0),dict(memory_gb=0),dict(outdir='.')]:
                path.write_text(json.dumps({**c,**change}))
                with self.assertRaises(ValueError):load_config(path)

    def test_slurm_quotes_shell_arguments(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'config;unexpected.json'
            text=slurm_script(path,dict(sample='sample',threads=4,memory_gb=32))
            self.assertIn("'"+str(path)+"'",text)
            self.assertIn('#SBATCH --mem=32G',text)
            import subprocess
            subprocess.run(['bash','-n'], input=text, text=True, check=True)

    def test_summary_uses_initial_not_calibrated_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); d=root/'baseline'/'classifications';d.mkdir(parents=True)
            (d/'x_calibrated_classification.tsv').write_text('Number of features: 10\nPredicted Class: GBM_RTK_II\nInitial Score: 0.06\nCalibrated Score: 0.9\n')
            report(root)
            import csv
            row=next(csv.DictReader((root/'classification_summary.csv').open()))
            self.assertEqual(float(row['score_0_to_1']),0.06)

    def test_status_missing_output_is_not_complete(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()) as output:
            main(['status', str(Path(tmp)/'missing')])
            self.assertFalse(json.loads(output.getvalue())['complete'])


if __name__ == '__main__':
    unittest.main()
