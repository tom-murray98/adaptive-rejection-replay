"""Command line interface; configuration paths are relative to the config file."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from . import __version__
from .runner import load_config, run, report


def slurm_script(config_path, c, hours=48):
    config_path = Path(config_path).resolve()
    logs = config_path.parent / 'logs'
    # Slurm directives do not use shell quoting; reject ambiguous paths.
    if any(ch.isspace() or ch in '\r\n' for ch in str(logs)) or '%' in str(logs):
        raise ValueError('Slurm log directory must have no whitespace or percent signs')
    if not 1 <= hours <= 720:
        raise ValueError('Slurm hours must be between 1 and 720')
    logs.mkdir(exist_ok=True)
    return f'''#!/usr/bin/env bash
#SBATCH --job-name=asreplay_{c['sample']}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={c['threads']}
#SBATCH --mem={c['memory_gb']}G
#SBATCH --time={hours}:00:00
#SBATCH --output={logs}/experiment-%j.out
#SBATCH --error={logs}/experiment-%j.err
set -euo pipefail
export NXF_SYNTAX_PARSER=v1
{shlex.quote(sys.executable)} -m asreplay run {shlex.quote(str(config_path))}
'''


def main(argv=None):
    parser = argparse.ArgumentParser(description='Paired adaptive sampling rejection replay')
    parser.add_argument('--version', action='version', version=__version__)
    subs = parser.add_subparsers(dest='command', required=True)
    init = subs.add_parser('init', help='Write a validated JSON config without starting a run')
    init.add_argument('config')
    for name in ('sample', 'bam', 'summary', 'ids', 'reference', 'workflow', 'outdir'):
        init.add_argument('--'+name, required=True)
    init.add_argument('--seconds', required=True, type=float)
    init.add_argument('--threads', type=int, default=16)
    init.add_argument('--memory-gb', type=int, default=128)
    init.add_argument('--fixed-overhead', type=float, default=0)
    init.add_argument('--container', default='graefox/scarlet:latest')
    init.add_argument('--keep-intermediates', action='store_true')
    for name in ('run', 'validate', 'submit'):
        p = subs.add_parser(name)
        p.add_argument('config')
        if name == 'submit':
            p.add_argument('--hours', type=int, default=48)
            p.add_argument('--dry-run', action='store_true', help='Write the Slurm script without submitting')
    p = subs.add_parser('prepare-workflow', help='Extract SCARLET methylation processes and patch RapidCNS2')
    p.add_argument('--scarlet', required=True)
    p.add_argument('--human-variation', required=True)
    p.add_argument('--outdir', required=True)
    for name in ('status', 'report'):
        p = subs.add_parser(name)
        p.add_argument('outdir')
    args = parser.parse_args(argv)
    try:
        if args.command == 'init':
            dest = Path(args.config).resolve()
            if dest.exists():
                raise ValueError('Config already exists; refusing to overwrite')
            c = {k: v for k, v in vars(args).items() if k not in ('command', 'config', 'keep_intermediates')}
            for key in ('bam', 'summary', 'ids', 'reference', 'workflow', 'outdir'):
                c[key] = str(Path(c[key]).expanduser().resolve())
            c['cleanup'] = not args.keep_intermediates
            if Path(c['outdir']) in dest.parents:
                raise ValueError('Keep configuration outside the new experiment output directory')
            # Validate in memory through a temporary sibling, without overwriting files.
            dest.parent.mkdir(parents=True, exist_ok=True)
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', dir=dest.parent) as f:
                json.dump(c, f); f.flush()
                load_config(f.name)
            with dest.open('x') as f:
                json.dump(c, f, indent=2); f.write('\n')
            print(dest)
        elif args.command in ('run', 'validate', 'submit'):
            c = load_config(args.config)
            if Path(c['outdir']).exists():
                raise ValueError('Output directory already exists; choose a fresh output directory')
            if args.command == 'run':
                os.environ.setdefault('NXF_SYNTAX_PARSER', 'v1')
                run(args.config)
            elif args.command == 'validate':
                print(json.dumps(c, indent=2))
            else:
                script = Path(args.config).resolve().with_suffix('.sbatch')
                text = slurm_script(args.config, c, args.hours)
                if script.exists() and script.read_text() != text:
                    raise ValueError('Slurm script exists with different contents')
                script.write_text(text)
                subprocess.run(['bash', '-n', str(script)], check=True)
                print(script)
                if not args.dry_run:
                    result = subprocess.run(['sbatch', '--parsable', str(script)], capture_output=True, text=True, check=True)
                    job_id = result.stdout.strip()
                    script.with_suffix('.submission.json').write_text(json.dumps({'job_id': job_id, 'config': str(Path(args.config).resolve())}, indent=2)+'\n')
                    print('Submitted job '+job_id)
        elif args.command == 'prepare-workflow':
            from .workflow import create
            print(create(args.scarlet, args.human_variation, args.outdir))
        elif args.command == 'report':
            report(Path(args.outdir))
        else:
            root = Path(args.outdir)
            status = {'complete': (root/'EXPERIMENT_SUCCESS').exists(), 'conditions': {}}
            if (root/'failure.json').exists():
                status['failure'] = json.loads((root/'failure.json').read_text())
            for c in ('baseline', 'trim_only', 'faster'):
                status['conditions'][c] = {k: (root/c/f).exists() for k, f in [('built','BUILD_SUCCESS'),('classified','CLASSIFICATION_SUCCESS'),('build_cleaned','build_cleanup.json'),('work_cleaned','classification_cleanup.json')]}
            print(json.dumps(status, indent=2))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Error: {error}\n')
