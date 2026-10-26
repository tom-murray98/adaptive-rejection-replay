"""Conservative cleanup: owned outputs, validated final BAMs, archived provenance."""
import json
import os
from pathlib import Path
import shutil
import tarfile
import pysam

CONDITIONS = ('baseline', 'trim_only', 'faster')
BUILD_FILES = ('original.bam', 'original.names.bam', 'aligned.names.bam',
               'restored.names.bam', 'reads.fastq.gz', 'tmp')
PATTERNS = ('*nanoDx_votes.txt', '*calibrated_classification.tsv', '*CAPPER_MODEL.csv')


def clean(root, condition, stage):
    root = Path(root).resolve()
    if condition not in CONDITIONS or stage not in ('build', 'classification'):
        raise ValueError('Unknown cleanup condition or stage')
    marker = root / '.asreplay-owned.json'
    if not marker.is_file() or marker.is_symlink():
        raise ValueError('Cleanup requires a package-owned experiment')
    if json.loads(marker.read_text())['outdir'] != str(root):
        raise ValueError('Ownership marker does not match output location')
    d = root / condition
    if d.is_symlink():
        raise ValueError('Condition directory must not be a symlink')
    bam = d / 'input.bam'
    if bam.is_symlink() or not (d / 'BUILD_SUCCESS').is_file() or not (d / 'input.bam.bai').is_file():
        raise ValueError('Final indexed BAM not validated')
    pysam.quickcheck(str(bam))
    audit = d / f'{stage}_cleanup.json'
    if audit.exists():
        return  # Idempotent cleanup of an already completed stage.
    targets = BUILD_FILES if stage == 'build' else ('work',)
    if stage == 'classification':
        if not (d / 'CLASSIFICATION_SUCCESS').is_file():
            raise ValueError('Classification has not completed')
        if not all(any(f.is_file() and f.stat().st_size for f in (d / 'classifications').glob(p)) for p in PATTERNS):
            raise ValueError('Missing published classifier results')
        # Published symlinks would become dangling when work is removed.
        for f in (d / 'classifications').rglob('*'):
            if f.is_symlink():
                if not f.is_file():
                    raise ValueError(f'Unsupported or broken published symlink: {f}')
                tmp = f.with_name(f.name + '.materializing')
                shutil.copy2(f, tmp)
                os.replace(tmp, f)
        with tarfile.open(d / 'task_provenance.tar.gz', 'x:gz') as archive:
            for directory, dirs, files in os.walk(d / 'work', followlinks=False):
                for name in files:
                    f = Path(directory) / name
                    if not f.is_symlink() and (name.startswith(('.command.', '.nextflow.log')) or name in ('.exitcode', 'params.json') or name.endswith('_probes_for_training.csv')):
                        archive.add(f, arcname=str(f.relative_to(d)), recursive=False)
    removed = []
    for name in targets:
        f = d / name
        if f.is_symlink():
            f.unlink()  # Never follow a symlink during deletion.
        elif f.is_file():
            removed.append({'path': name, 'bytes': f.stat().st_size})
            f.unlink()
        elif f.is_dir():
            size = sum(Path(dp, fn).stat().st_size for dp, _, names in os.walk(f, followlinks=False) for fn in names if not Path(dp, fn).is_symlink())
            removed.append({'path': name, 'bytes': size})
            shutil.rmtree(f)
    pysam.quickcheck(str(bam))
    audit.write_text(json.dumps({'stage': stage, 'removed': removed, 'final_bam_quickcheck': 'passed'}, indent=2) + '\n')
