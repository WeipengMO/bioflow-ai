# External tools for hic_hichip

This directory is reserved for external tools that are not cleanly installed as
standard Snakemake conda environments.

FitHiChIP can be installed here for workflow-level development:

```bash
scripts/install_fithichip.sh
```

For real projects, prefer a project-local install so deployment stays
self-contained:

```bash
./scripts/install_fithichip.sh --prefix tools/FitHiChIP
```

Downloaded tool source directories such as `FitHiChIP/` are intentionally not
tracked by git.
