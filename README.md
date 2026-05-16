# BioFlowAI

**Agent-ready Snakemake workflows for bioinformatics.**

BioFlowAI is a collection of reusable Snakemake workflows for bioinformatics analysis. It is refactored from the Snakemake pipelines previously maintained inside `mowp_scripts`, with a focus on modularity, reproducibility, clarity, and future AI-agent integration.

This repository is not a single analysis pipeline. It is a workflow collection. Each workflow should be self-contained, documented, configurable, and readable by both humans and AI agents.

---

## Goals

BioFlowAI aims to provide:

- reusable Snakemake workflows for omics analysis
- clear project structure for long-term maintenance
- reproducible environments and config-driven execution
- standardized inputs, outputs, logs, and reports
- documentation that is easy for humans and AI agents to understand
- a foundation for future agent-assisted workflow execution and debugging

---

## Repository Concept

The original `mowp_scripts` repository contains a mixture of analysis scripts, visualization utilities, notebooks, and Snakemake pipelines.

BioFlowAI separates the reusable Snakemake workflows into a dedicated repository.

Project-specific scripts can stay outside this repository. Workflows that are reusable, generalizable, or worth maintaining should be moved here.

---

## Basic Layout

```text
BioFlowAI/
├── README.md
├── workflows/
│   └── <workflow_name>/
│       ├── README.md
│       ├── Snakefile
│       ├── config/
│       ├── rules/
│       ├── envs/
│       ├── scripts/
│       ├── docs/
│       ├── agent/
│       └── tests/
├── shared/
├── docs/
└── tests/
```

Each workflow should live under `workflows/<workflow_name>/` and include its own README, example config, rule files, environments, documentation, and agent metadata.

---

## Workflow Design

Each workflow should be designed around a few simple rules.

### 1. Self-contained workflow directory

A workflow should be runnable from its own directory.

It should include:

- `Snakefile`
- `config/`
- `rules/`
- `envs/`
- `scripts/`
- `docs/`
- `agent/`

### 2. Config-driven execution

Project-specific values should be placed in config files instead of being hard-coded in rules.

Example:

```yaml
genome: hg38
samples: config/samples.tsv
outdir: results
threads: 16
```

### 3. Clear inputs and outputs

Each workflow should document:

- required input files
- required metadata columns
- main output files
- report files
- log files
- temporary or regenerable files

### 4. Modular rules

Large monolithic Snakefiles should be avoided.

Rules should be split by function, for example:

```text
rules/common.smk
rules/preprocess.smk
rules/align.smk
rules/qc.smk
rules/report.smk
```

### 5. Reproducible environments

Each workflow should define Conda environments under `envs/`.

Do not rely on packages installed in the user's current shell environment.

---

## Agent-Ready Design

BioFlowAI is designed so that future AI agents can inspect, run, debug, and refactor workflows safely.

Each workflow should provide an `agent/` directory containing structured metadata.

Recommended files:

```text
agent/
├── manifest.yml
├── io.contract.yml
└── task_templates.md
```

The agent metadata should answer:

- What does this workflow do?
- What files are required before running it?
- Which config keys are required?
- Which outputs are final results?
- Where are the logs?
- Which files are user-provided and should not be modified?
- Which files are safe to regenerate?
- What is the safest dry-run command?

Example `agent/manifest.yml`:

```yaml
workflow:
  name: example_workflow
  description: Short description of the workflow.
  status: development

entrypoint:
  snakefile: Snakefile
  example_config: config/config.example.yml

commands:
  dry_run: snakemake -np
  run: snakemake --use-conda -j 16
  unlock: snakemake --unlock

inputs:
  config:
    path: config/config.yml
    required: true
  samples:
    path: config/samples.tsv
    required: true

outputs:
  results:
    path: results/
  logs:
    path: logs/
  reports:
    path: reports/

safety:
  user_owned:
    - raw data
    - config files
    - sample metadata
  safe_to_regenerate:
    - logs
    - benchmark files
    - rule-generated results
```

---

## Recommended Workflow README

Each workflow should have its own README with this minimal structure:

```markdown
# <Workflow Name>

## Overview
## Inputs
## Configuration
## Outputs
## Quick Start
## Troubleshooting
## Notes for AI Agents
```

The root README should stay general. Workflow-specific methods, parameters, and biological interpretation should be documented inside each workflow directory.

---

## Development Guidelines

Use clear and predictable names:

```text
workflows/chipseq
workflows/atacseq
workflows/rnaseq
workflows/wgbs
workflows/hic
```

Use stable internal structure:

```text
Snakefile
config/config.example.yml
config/samples.example.tsv
rules/*.smk
envs/*.yml
scripts/*.py
scripts/*.R
docs/*.md
agent/*.yml
```

Shell rules should be strict when possible:

```bash
set -euo pipefail
```

Important outputs should be checked when appropriate:

```bash
test -s output.file
```

Each rule should write logs to `logs/` whenever possible.

---

## For AI Agents

Agents should treat this repository as a workflow collection, not a loose script folder.

When working on a workflow, an agent should:

1. Identify the target workflow directory.
2. Read the workflow README.
3. Read `agent/manifest.yml` if available.
4. Compare the user config with the example config.
5. Inspect the `Snakefile` and included rule files.
6. Run or suggest `snakemake -np` before changing logic.
7. Check rule-specific logs after failures.
8. Propose minimal and reversible changes.

Agents should avoid:

- changing raw data
- silently editing user config files
- assuming genome build or sequencing mode
- mixing genome resources from different assemblies
- adding workflow-specific details to the root README

---

## License

MIT License
