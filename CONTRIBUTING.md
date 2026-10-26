# Contributing

Use a virtual environment and install with `python -m pip install -e .`.
Run `python -m unittest discover -s tests -v` before proposing changes.

Changes to read selection, timing, modification clipping or cleanup require
regression tests with small synthetic fixtures. Preserve original read IDs,
sequencing orientation and modification-tag meaning. Cleanup changes must prove
that final results and input files remain intact.

Do not commit patient identifiers, clinical tables, sequencing data, trained
models, local configuration or generated experiment results. Use generic paths
and fabricated sample labels in examples. Describe scientific assumptions and
external dependency requirements in pull requests.

For upstream compatibility issues, include tool versions and a minimal synthetic
reproducer. Do not upload sensitive data in issue reports.
