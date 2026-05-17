# BioFlowAI

**Agent-ready Snakemake workflows for bioinformatics.**

BioFlowAI is a collection of reusable Snakemake workflows for bioinformatics analysis. It is refactored from the Snakemake pipelines previously maintained inside `mowp_scripts`, with a focus on modularity, reproducibility, clarity, and AI-agent integration.

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
│       └── agent/
├── shared/
└── tests/
```

Each workflow should live under `workflows/<workflow_name>/` and include its own README, example config, rule files, environments, documentation, and agent metadata.

---

## Workflow Design

Each workflow should be designed as a self-contained project template.

### 1. Self-contained workflow directory

A workflow should be runnable and deployable from its own directory.

Recommended contents:

- `README.md`: workflow overview, required inputs, expected outputs, and run examples
- `Snakefile`
- `config/`: example config files and project-specific metadata templates
- `rules/`: modular Snakemake rule files grouped by function
- `envs/`: Conda environment definitions used by the workflow
- `scripts/`: helper scripts such as deployment, validation, and data-prep utilities
- `agent/`: machine-readable metadata for AI-assisted deployment, execution, and debugging
- `run_snakemake.sh`: a stable wrapper command for dry-run and production execution

This structure keeps workflow logic reusable while allowing deployment into a separate project data directory.

---

## Teaching AI To Deploy A Workflow

If you want an AI agent to deploy a workflow into the correct data directory, the prompt needs to describe the target location, the source workflow, and the constraints clearly.

The minimum information to provide is:

- which workflow to deploy, such as `workflows/chipseq` or `workflows/rnaseq`
- the target project directory, such as `/data/project_A/chipseq_run_01`
- where the raw data already lives, or where the workflow should expect it
- which config files must be created or edited
- whether the agent is allowed to run only a dry-run or a full execution
- whether existing files in the target directory should be preserved

The AI should be guided to follow this sequence:

1. Read the workflow `README.md` and `agent/` metadata first.
2. Use the workflow's deployment script or documented deployment method.
3. Deploy into the target data directory instead of editing files under `workflows/` directly.
4. Fill in project config with paths that match the target directory.
5. Run a Snakemake dry-run before attempting a real execution.
6. Report which files were created, which files still need user input, and which command should be run next.

Recommended prompt in Chinese:

```text
请阅读 workflows/chipseq/README.md 以及 agent/ 目录中的说明，把这个 workflow 部署到 /data/project_A/chipseq_run_01。

要求：
1. 不要直接修改 workflows/chipseq 源目录中的模板内容。
2. 优先使用该 workflow 自带的部署脚本或 README 中说明的部署方式。
3. 在目标目录中准备好 config、raw_data、运行脚本和必要的链接或复制文件。
4. 根据目标目录调整配置文件中的路径。
5. 先执行 dry-run，不要直接正式运行。
6. 最后告诉我：创建了哪些文件、还需要我补充哪些输入、下一条推荐执行命令是什么。
```

Recommended prompt in English:

```text
Read workflows/chipseq/README.md and the files under agent/, then deploy this workflow into /data/project_A/chipseq_run_01.

Requirements:
1. Do not modify the source template under workflows/chipseq.
2. Prefer the workflow's own deployment script or the documented deployment procedure.
3. Prepare the target directory with config, raw_data, run scripts, and any required symlinks or copied files.
4. Update configuration paths so they match the target project directory.
5. Run a Snakemake dry-run only; do not start the full workflow.
6. At the end, report which files were created, which inputs are still missing, and the next recommended command.
```

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

## Inputs
## Configuration
## Outputs
## Workflow Logic
## Quick Start
## Troubleshooting
## Notes for AI Agents
## External Resource Links (if applicable)
```

The root README should stay general. Workflow-specific methods, parameters, and biological interpretation should be documented inside each workflow directory.

---

## Development Guidelines

Use clear and predictable names:

```text
workflows/chipseq
workflows/atacseq
workflows/rnaseq
...
```

Use stable internal structure:

```text
Snakefile
config/config.example.yml
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
