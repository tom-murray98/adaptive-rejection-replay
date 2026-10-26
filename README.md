# AS off-target replay

Simulate 'earlier' adaptive sampling unblocking of off-target reads from an existing ONT run, then compare methylation classification within the same sequencing time.

The package constructs **baseline**, **trim-only**, and **faster-rejection** BAMs, realigns them, restores retained modification tags, and runs SCARLET's Sturgeon, nanoDx and RapidCNS2 methylation branches. It supports Slurm submission and cleans validated intermediate files to reduce bulky file size.


## Requirements

- `samtools`, `minimap2`, Nextflow and Docker
- An existing SCARLET checkout of EPI2ME `wf-human-variation`.
- Matching BAM, sequencing summary fiel, unblock-ID list and indexed hg38 reference (same as one used in og basecalling (usually hg38_simple.fa for samples before 2025)).
  


## Install

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install .
as-replay --help
```

## Prepare the methylation workflow

Extraction reads your local upstream checkouts and writes a fresh workflow directory:

```bash
as-replay prepare-workflow \
  --scarlet /path/to/SCARLET \
  --human-variation /path/to/wf-human-variation \
  --outdir /path/to/methylation-workflow
```

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
  --threads 16

as-replay validate experiment-1h.json
```

## Submit to Slurm

```bash

# Submit; prints the job ID and saves submission metadata
as-replay submit experiment-1h.json --hours 48

# Inspect markers and recorded failures
as-replay status /path/to/new-output-1h
```

The batch script uses the active Python interpreter and inherited `PATH`, ensure those tools are available on the compute node.

## Storage and cleanup

Cleanup is **on by default**.

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

Sturgeon, nanoDx and RapidCNS2 scores have different interpretations, so don't compare their results 1 to 1.


