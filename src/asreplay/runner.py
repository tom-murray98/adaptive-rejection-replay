"""Portable orchestration of the validated replay and methylation workflow."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import uuid

from . import __version__
from .simulation import prepare, export
from .cleanup import clean, CONDITIONS, PATTERNS


def load_config(path):
    config = json.loads(Path(path).read_text())
    required = ('sample', 'bam', 'summary', 'ids', 'reference', 'workflow', 'outdir', 'seconds')
    for key in required:
        if key not in config:
            raise ValueError(f'Missing config field: {key}')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', config['sample']):
        raise ValueError('Sample name must contain letters, numbers, dots, underscores or hyphens')
    config.setdefault('threads', 16)
    config.setdefault('memory_gb', 128)
    config.setdefault('container', 'graefox/scarlet:latest')
    config.setdefault('fixed_overhead', 0)
    config.setdefault('cleanup', True)
    config.setdefault('reject_reason', 'data_service_unblock_mux_change')
    import math
    if not math.isfinite(config['seconds']) or config['seconds'] <= 0:
        raise ValueError('seconds must be finite and positive')
    if not math.isfinite(config['fixed_overhead']) or config['fixed_overhead'] < 0:
        raise ValueError('fixed_overhead must be finite and nonnegative')
    if not isinstance(config['threads'], int) or config['threads'] < 4:
        raise ValueError('threads must be an integer >= 4')
    if not isinstance(config['memory_gb'], int) or config['memory_gb'] < 4:
        raise ValueError('memory_gb must be an integer >= 4')
    if not isinstance(config['cleanup'], bool):
        raise ValueError('cleanup must be a JSON boolean')
    for key in ('bam', 'summary', 'ids', 'reference', 'workflow', 'outdir'):
        p = Path(config[key]).expanduser()
        if not p.is_absolute():
            p = Path(path).resolve().parent / p
        config[key] = str(p.resolve())
    out = Path(config['outdir'])
    for key in ('bam', 'summary', 'ids', 'reference', 'workflow'):
        p = Path(config[key])
        if not p.exists():
            raise ValueError(f'{key} not found: {p}')
        if p == out or out in p.parents:
            raise ValueError('Inputs must be outside the output directory')
    workflow = Path(config['workflow'])
    for name in ('main.nf', 'rapidcns2.R', 'bin', 'wf-human-variation'):
        if not (workflow / name).exists():
            raise ValueError(f'Workflow missing {name}; run prepare-workflow first')
    return config


def run_command(args, log):
    args = [str(a) for a in args]
    with log.open('a') as f:
        f.write(json.dumps(args) + '\n')
    subprocess.run(args, check=True)


def collect(root, condition):
    d = root / condition
    dest = d / 'classifications'
    dest.mkdir(exist_ok=True)
    patterns = PATTERNS + ('*nanoDx_output.txt', '*rf_details.tsv', '*votes.tsv', '*CAPPER_MODEL.pdf')
    for pattern in patterns:
        found = {}
        for f in (d / 'work').rglob(pattern):
            if not f.is_file():
                continue
            if f.name in found and f.read_bytes() != found[f.name].read_bytes():
                raise ValueError(f'Conflicting classifier outputs: {f.name}')
            found[f.name] = f
        for name, source in found.items():
            target = dest / name
            if target.is_symlink():
                target.unlink()
            shutil.copy2(source, target)
    if not all(any(f.is_file() and f.stat().st_size for f in dest.glob(p)) for p in PATTERNS):
        raise ValueError('Incomplete classifier outputs')


def report(root):
    root = Path(root)
    rows = []
    for condition in CONDITIONS:
        d = root / condition / 'classifications'
        files = list(d.glob('*nanoDx_votes.txt'))
        if files:
            for row in csv.DictReader(files[0].open(), delimiter='\t'):
                rows.append({'condition': condition, 'class': row['class'], 'score': row['score'], 'num_features': row['num_features']})
    with (root / 'nanodx_comparison.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=('condition', 'class', 'score', 'num_features'), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    combined = []
    for condition in CONDITIONS:
        d = root / condition / 'classifications'
        for pattern, classifier in (('*nanoDx_votes.txt', 'nanoDx'), ('*CAPPER_MODEL.csv', 'Sturgeon'), ('*calibrated_classification.tsv', 'RapidCNS2 initial')):
            files = list(d.glob(pattern))
            if not files:
                continue
            source = files[0]
            if classifier == 'nanoDx':
                r = max(csv.DictReader(source.open(), delimiter='\t'), key=lambda r: float(r['score']))
                label, score, features = r['class'], float(r['score']), int(r['num_features'])
            elif classifier == 'Sturgeon':
                r = next(csv.DictReader(source.open()))
                label, score = max(((k, float(v)) for k, v in r.items() if k != 'number_probes'), key=lambda x: x[1])
                features = int(r['number_probes'])
            else:
                r = dict(line.split(': ', 1) for line in source.read_text().splitlines() if ': ' in line)
                label, score, features = r['Predicted Class'], float(r['Initial Score']), int(r['Number of features'])
            combined.append({'condition': condition, 'classifier': classifier, 'class': label, 'score_0_to_1': score, 'features': features, 'source': str(source.relative_to(root))})
    with (root / 'classification_summary.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=('condition', 'classifier', 'class', 'score_0_to_1', 'features', 'source'))
        writer.writeheader()
        writer.writerows(combined)


def run(config_path):
    c = load_config(config_path)
    for tool in ('samtools', 'minimap2', 'nextflow', 'docker'):
        if not shutil.which(tool):
            raise ValueError(f'Required executable not on PATH: {tool}')
    root = Path(c['outdir'])
    root.mkdir(parents=True, exist_ok=False)
    token = str(uuid.uuid4())
    (root / '.asreplay-owned.json').write_text(json.dumps({'outdir': str(root), 'token': token}))
    (root / 'config.json').write_text(json.dumps(c, indent=2) + '\n')
    log = root / 'commands.jsonl'
    try:
        run_command(['docker', 'info'], log)
        run_command(['samtools', 'quickcheck', c['bam']], log)
        workflow = root / 'workflow'
        shutil.copytree(c['workflow'], workflow, symlinks=True)
        from .rapidcns2 import repair
        script = workflow / 'rapidcns2.R'
        if script.is_symlink():
            text = script.read_text(); script.unlink(); script.write_text(text)
        script.write_text(repair(script.read_text()))
        provenance = {'package_version': __version__, 'inputs': {}, 'workflow_files': {}}
        for key in ('bam', 'summary', 'ids', 'reference'):
            f = Path(c[key]); provenance['inputs'][key] = {'path': str(f), 'bytes': f.stat().st_size, 'mtime_ns': f.stat().st_mtime_ns}
        for f in workflow.glob('*'):
            if f.is_file():
                provenance['workflow_files'][f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
        for tool in ('samtools', 'minimap2', 'nextflow'):
            result = subprocess.run([tool, '--version'], capture_output=True, text=True, errors='replace', check=True)
            provenance[tool] = (result.stdout + result.stderr).strip()
        image = subprocess.run(['docker', 'image', 'inspect', c['container'], '--format', '{{json .RepoDigests}}'], capture_output=True, text=True)
        provenance['container_repo_digests'] = image.stdout.strip()
        (root / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
        prepare(argparse.Namespace(command='prepare', summary=c['summary'], ids=c['ids'], seconds=c['seconds'], fraction=0.5, fixed_overhead=c['fixed_overhead'], reject_reason=c['reject_reason'], outdir=str(root)))
        export(argparse.Namespace(bam=c['bam'], outdir=str(root)))
        config = root / 'resources.config'
        # JSON strings are valid quoted Groovy strings; forbid Groovy interpolation.
        if '$' in c['container']:
            raise ValueError('Container name must not contain interpolation')
        config.write_text(f'''params.threads = {c['threads']}
process.cpus = {c['threads']}
process.memory = '{c['memory_gb']} GB'
process.container = {json.dumps(c['container'])}
process {{ withName: human_variation_mods {{ container = null }} }}
executor.queueSize = 1
executor.cpus = {c['threads']}
executor.memory = '{c['memory_gb']} GB'
docker.runOptions = '-u {__import__('os').getuid()}:{__import__('os').getgid()} --cpus={c['threads']} --memory={c['memory_gb']-1}g'
''')
        (workflow / 'nested.config').write_text(f"""process.cpus = {c['threads']}
process.memory = '{min(32,c['memory_gb'])} GB'
executor.queueSize = 1
executor.cpus = {c['threads']}
executor.memory = '{c['memory_gb']} GB'
docker.runOptions = '-u {__import__('os').getuid()}:{__import__('os').getgid()} --cpus={c['threads']} --memory={min(32,c['memory_gb']-1)}g'
env {{ OMP_NUM_THREADS='1'; OPENBLAS_NUM_THREADS='1'; MKL_NUM_THREADS='1' }}
""")
        for condition in CONDITIONS:
            d = root / condition
            (d / 'tmp').mkdir()
            (d / 'work').mkdir()
            (d / 'classifications').mkdir()
            # Argument arrays and explicit pipe handles avoid shell interpolation.
            mapping = ['minimap2', '-t', str(c['threads']-2), '-ax', 'map-ont', '-L', '-R', (d / 'read_group.txt').read_text().strip(), c['reference'], str(d / 'reads.fastq.gz')]
            sort = ['samtools', 'sort', '-n', '-@', '1', '-m', '1G', '-T', str(d/'tmp/align'), '-o', str(d/'aligned.names.bam'), '-']
            with log.open('a') as f:
                f.write(json.dumps({'pipe': [mapping, sort]})+'\n')
            with subprocess.Popen(mapping, stdout=subprocess.PIPE) as process:
                result = subprocess.run(sort, stdin=process.stdout)
                process.stdout.close()
                code = process.wait()
                if result.returncode or code:
                    raise RuntimeError(f'Alignment pipeline failed: minimap2={code}, samtools={result.returncode}')
            run_command(['samtools', 'sort', '-n', '-@', '1', '-m', '1G', '-T', d/'tmp/original', '-o', d/'original.names.bam', d/'original.bam'], log)
            run_command([sys.executable, '-m', 'asreplay.modifications', '--original', d/'original.names.bam', '--shortened', d/'aligned.names.bam', '--output', d/'restored.names.bam', '--stats', d/'modification_stats.json'], log)
            run_command(['samtools', 'sort', '-@', '1', '-m', '1G', '-T', d/'tmp/final', '-o', d/'input.bam', d/'restored.names.bam'], log)
            run_command(['samtools', 'quickcheck', d/'input.bam'], log)
            run_command(['samtools', 'index', d/'input.bam'], log)
            (d / 'BUILD_SUCCESS').touch()
            if c['cleanup']:
                clean(root, condition, 'build')
            run_command(['nextflow', '-c', config, '-log', d/'nextflow.log', 'run', workflow/'main.nf', '-work-dir', d/'work', '-with-report', d/'execution-report.html', '--bam', d/'input.bam', '--reference', c['reference'], '--sample', c['sample']+'_'+condition, '--threads', c['threads'], '--outdir', d/'classifications'], log)
            collect(root, condition)
            (d / 'CLASSIFICATION_SUCCESS').touch()
            report(root)
            if c['cleanup']:
                clean(root, condition, 'classification')
        (root / 'EXPERIMENT_SUCCESS').touch()
    except BaseException as error:
        (root / 'failure.json').write_text(json.dumps({'error_type': type(error).__name__, 'error': str(error)}, indent=2)+'\n')
        raise
