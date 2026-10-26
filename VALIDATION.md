# Validation record

Package validation performed locally on 22 October 2026:

- 16 tests passed with `ASREPLAY_INTEGRATION_TOOLCHAIN=1`, including the full synthetic three-condition pipeline with real samtools/minimap2 and stub classifier executables.
- Built an installable wheel and installed it into an isolated Python environment with existing system dependencies.
- Checked the installed `as-replay --version` and `--help` entrypoints.
- Extracted a workflow from the existing SCARLET checkout and Git-tracked human-variation asset into a temporary directory.
- Created and validated a configuration against existing input paths; no experiment was started.
- Generated a Slurm script using `submit --dry-run`; shell syntax validation passed. No new Slurm job was submitted.
- Parsed the resource overrides with real Nextflow and confirmed the effective CPU, memory and container settings.
- Searched the package for local account paths and sample-specific identifiers; none were found in the distribution sources.

The packaged interface has not completed a new full real-data classifier run.
Regression tests use classifier stubs for orchestration, not model predictions.
The legacy experiment scripts and completed experiment folders were not changed
or rerun during packaging. The separate SimReadUntil code is not included.
