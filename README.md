# Adaptive rejection replay

Simulate earlier adaptive sampling rejection from an existing Oxford Nanopore run, then compare methylation classification within the same sequencing time.

The package constructs **baseline**, **trim-only**, and **faster-rejection** BAMs, realigns them, restores retained modification tags, and runs SCARLET's Sturgeon, nanoDx and RapidCNS2 methylation branches. It supports Slurm submission and cleans validated intermediate files condition by condition.

**Status:** research software, version 0.1.0. The underlying scripts have been used on real experiments. This packaged interface has automated regression tests, but has not yet completed a new full real-data run as a package. It is not a diagnostic device. No sample data, clinical tables, pretrained models or upstream repositories are distributed here.

## What is simulated?

| Condition | Molecules | Confirmed rejected reads |
| --- | --- | --- |
| Baseline | Finish within the observed time window | Original sequences and durations |
| Trim only | Exactly the baseline read IDs | First half of each sequence; no additional molecules |
| Faster rejection | Finish within the accelerated channel timeline | First half of each sequence and reduced duration; later recorded molecules can enter the window |

A confirmed rejection must have a read or parent ID in the matching unblock list **and** sequencing-summary `end_reason=data_service_unblock_mux_change`. This is a summary field value, not a BAM optional tag. Rejected reads are a proxy for off-target molecules; target status is not recomputed. `mux_change`, `unblock_mux_change`, and pass/fail folders alone are not the rejection rule.

Events are replayed in order per channel, preserving observed gaps. With duration `d` and fixed overhead `o`, a rejected event takes:

```text
min(o, d) + 0.5 * (d - min(o, d))
```

`fixed_overhead=0` halves the entire rejected duration. This is an idealised intervention, not a measurement of a twofold computational speedup. Sequences retain `max(1, floor(length/2))` bases in sequencing orientation. Other events are unchanged. Only completed molecules are exported. Later reads are neither duplicated nor sampled with replacement.

## Requirements

- Linux, Python 3.10 or newer, and `pysam` (installed with the package).
- `samtools`, `minimap2`, Nextflow and a working Docker daemon on `PATH` for actual experiments.
- Docker access for the submitting user **on the compute node**. The package does not modify permissions or automatically invoke `sudo`/`sg`.
- An existing SCARLET checkout with its assets and a Git checkout of EPI2ME `wf-human-variation`. Models, R packages and classifier executables must be available in the chosen container and SCARLET assets. Reference/annotation/model use is subject to upstream licences.
- Matching original BAM, sequencing summary, unblock-ID list and indexed hg38 reference. The current exporter requires one BAM read group and full primary sequences, quality values, and `MM`/`ML` modification tags.
- Slurm `sbatch` for submission; local execution does not require Slurm.

The original experiments used samtools 1.21, minimap2 2.31, Nextflow 26.04.4 and `graefox/scarlet:latest`. Other combinations need validation. `latest` is a convenient default, **not a reproducibility lock**: use a validated image digest for final analyses, and pin compatible upstream revisions.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
as-replay --help
```

For development: `python -m pip install -e .`.

## Prepare the methylation workflow

Extraction reads your local upstream checkouts and writes a fresh workflow directory:

```bash
as-replay prepare-workflow \
  --scarlet /path/to/SCARLET \
  --human-variation /path/to/wf-human-variation \
  --outdir /path/to/methylation-workflow
```

This extracts the supported process definitions, snapshots the human-variation source and records its Git revision. SCARLET's `bin` assets remain referenced through a symlink; keep that checkout and its models available. The corrected RapidCNS2 script filters modification-code `m` before discarding the code column, preventing duplicate `h` and `m` observations from competing for the first probe match. The workflow also limits R prediction threads.

Extraction intentionally fails if expected SCARLET processes are missing; it does not guess how to adapt a changed upstream interface. When reusing an existing extracted workflow, `run` applies the same checked RapidCNS2 patch to a private snapshot. The package does not change the supplied workflow or upstream checkouts.

## Configure an experiment

```bash
as-replay init experiment-1h.json \
  --sample example_1h \
  --bam /path/to/original.sorted.bam \
  --summary /path/to/sequencing_summary.txt.gz \
  --ids /path/to/unblocked_read_ids.txt.gz \
  --reference /path/to/hg38.fa \
  --workflow /path/to/methylation-workflow \
  --outdir /path/to/new-output-1h \
  --seconds 3600 \
  --threads 16 --memory-gb 128

as-replay validate experiment-1h.json
```

Use `1800`, `21600`, or `43200` seconds for 30 minutes, six hours or 12 hours. Give every window a fresh config and output directory. Relative paths in a JSON config resolve against that config's directory; `init` writes absolute paths. An editable example is in `examples/config.example.json`.

`validate` checks configuration and input existence; it does not scan the full BAM/summary or guarantee compatibility. During preparation, duplicate read IDs and overlapping channel events cause a deliberate stop. The output directory must not already exist. Inputs are read-only and must be outside the output directory.

## Run locally or submit to Slurm

```bash
# Local execution
as-replay run experiment-1h.json

# Inspect the generated batch script without submitting
as-replay submit experiment-1h.json --dry-run

# Submit; prints the job ID and saves submission metadata
as-replay submit experiment-1h.json --hours 48

# Inspect markers and recorded failures
as-replay status /path/to/new-output-1h
```

The batch script uses the active Python interpreter and inherited `PATH`; ensure those tools are available on the compute node. It requests the configured cores and memory. The example 128 GB allocation is deliberately conservative, not a measured minimum. Logs are written beside the configuration under `logs/`. Slurm log directory paths cannot contain whitespace or percent signs.

## Storage and cleanup

Cleanup is **on by default**. Initial export writes the three read sets together. Subsequent alignment, modification restoration and classification proceed one condition at a time. Thus export still requires space for all unprocessed read sets, plus the active condition's alignment/classification work.

1. Validate the final BAM and create its index.
2. Delete that condition's original-subset copies, FASTQ, name-sorted alignment/restoration BAMs and sorting temporary files.
3. Run all three classifier branches and validate their published outputs.
4. Materialise published result symlinks, archive task commands/logs and delete the completed condition's Nextflow work directory.

Final indexed BAMs, classifier results, reports, command logs, workflow snapshots, the selection database and provenance remain. Cleanup requires the package ownership marker and successful stage markers; it never follows symlinks during deletion. Each cleanup writes an audit JSON file. Use `--keep-intermediates` at `init` to disable cleanup.

A failed stage retains its unfinished intermediates and records `failure.json`. Automatic recovery/resume is not implemented: inspect the logs and use a new output directory for a fresh run. Once successful tasks have been cleaned, their cache is not available for Nextflow resume. The input data remain available for regeneration. Do not run the old experiment-specific cleanup scripts against arbitrary directories.

## Outputs and scores

```text
output/
  config.json, provenance.json, commands.jsonl
  plan.json, selection.sqlite, export.json
  resources.config, workflow/
  classification_summary.csv, nanodx_comparison.tsv
  baseline/                 # also trim_only/ and faster/
    input.bam, input.bam.bai
    modification_stats.json
    classifications/
    nextflow.log, execution-report.html
    task_provenance.tar.gz
    build_cleanup.json, classification_cleanup.json
  EXPERIMENT_SUCCESS        # only after all conditions succeed
```

`classification_summary.csv` contains the leading class, native 0–1 score and feature count for each classifier and condition. RapidCNS2 uses the **initial** score, not calibrated output. Raw calibrated outputs remain retained for audit. Regenerate summaries without alignment or classification:

```bash
as-replay report /path/to/completed-output
```

Sturgeon, nanoDx and RapidCNS2 scores have different interpretations; do not compare their numeric magnitudes or treat them as diagnostic correctness probabilities. RapidCNS2 selects at most 10,000 features, so identical capped counts do not imply identical probe sets. Historical ROBIN random-forest score files may contain calibrated scores and are not automatically equivalent to the initial score used here.

## Limitations

- Recorded rejection decisions are unchanged; shorter-prefix mapping accuracy is not tested.
- Modification calls are clipped from existing calls, not re-inferred from shortened signal.
- Gaps, encounter order and pore behaviour remain tied to the observed run; clock-based mux scans are not independently rescheduled.
- Only complete reads count; partial molecules at the cutoff are excluded.
- Reserve is finite. Review `reserve_exhausted_channels` in `plan.json`; inactivity and insufficient held-out data can both cause flags. Full-run replay cannot invent later molecules.
- Summary events absent from the BAM consume time but contribute no classifier observations; see `export.json`.
- Time windows from one run are correlated observations, not independent biological replicates. No confidence intervals or reference diagnosis are inferred.
- Random-forest training does not set an explicit seed, matching the original implementation; run-to-run variation is possible.

## Development and validation

```bash
python -m unittest discover -s tests -v
```

Tests cover paired selection and time savings, preserved gaps/overheads, reverse-strand and multi-code modification restoration, invalid events, RapidCNS2 patching, configuration validation, Slurm shell quoting, initial-score summaries, and cleanup refusal/preservation. Tests use synthetic data. They do not download models, access a Slurm cluster, or execute full classifier containers. With samtools and minimap2 on PATH, enable the full three-condition orchestration test with `ASREPLAY_INTEGRATION_TOOLCHAIN=1 python -m unittest discover -s tests -v`; it uses real alignment tools and synthetic classifier stubs. GitHub Actions runs them on Python 3.10 and 3.12.

Before publishing scientific results, validate the packaged workflow with a small real-data run, pin external dependencies, verify classifier input conventions and assess repeatability. The separate SimReadUntil pilot and sample-specific reports are not part of this package.

## Licence and upstream resources

Package code is MIT licensed; see `LICENSE`. SCARLET, ROBIN, wf-human-variation, classifier code, Docker images, training data and models remain separately licensed. This repository does not grant rights to redistribute those resources and does not vendor them. Review upstream licences before sharing extracted workflows or assets.
