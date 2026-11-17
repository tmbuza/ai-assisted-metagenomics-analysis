# AI-Assisted Metagenomics Analysis

A CDI learning project for exploring microbial composition, Shannon diversity, and functional profiles with Python. Use AI to help explain code and draft interpretations, then check the inputs, calculations, and claims yourself.

Start with the complete synthetic practice. The hypertension case study is a preparation section: real taxonomic profiles, functional units, source attribution, and study design must be verified before analysis. Your own study uses the same input contract.

## 1. Open the project

Download and unzip this folder, then open it in VS Code. Open a terminal in `ai-assisted-metagenomics-analysis`, the folder containing this README. Run all commands below there. You do not need to clone a repository to use this download.

## 2. Create a Python environment

Use Python 3.10 or later. On macOS or Linux:

```bash
python3 --version
python3 -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
py --version
py -m venv .venv
.venv\Scripts\Activate.ps1
```

With the environment active, install the analysis packages:

```bash
pip install --upgrade pip
python -m pip install -r requirements.txt
```

## 3. Check the inputs

```bash
python scripts/python/run_workflow.py --data-dir data/synthetic --output-dir results/synthetic --check-only
```

Expected message: `Input checks passed: 12 samples; 2 groups.`

## 4. Run the analysis

```bash
python scripts/python/run_workflow.py --data-dir data/synthetic --output-dir results/synthetic
```

This creates composition, Shannon diversity, and functional-profile figures in `results/synthetic/figures`, and tables plus a run record in `results/synthetic/tables`. Bundled figures let you preview the result before running the code.

## 5. Verify one sample

```bash
python scripts/python/verify_sample.py data/synthetic results/synthetic S01
```

The verification independently recalculates the sample total, number of detected taxa, and Shannon diversity, then compares them with the saved table. A mismatch exits with an error. These synthetic abundance units are illustrative, not sequencing-read counts.

## 6. Explain and report

Follow `examples/synthetic/00-start-here.qmd` through `06-export-the-report.qmd`. Edit `reports/synthetic-report.qmd`, replacing the prompts with your interpretation and checked values.

Quarto is a separate application; it is not installed by pip. If Quarto is installed:

```bash
quarto render reports/synthetic-report.qmd --to html
quarto preview
```

To build the whole guide:

```bash
quarto render
```

The book is written to `docs/`. Python analysis runs explicitly in the terminal; Quarto pages do not execute analysis code.

## Input contract

All input files are tab-separated. Taxa and functional features are rows; sample IDs are columns. Column names and sample IDs must be unique. All three files must have exactly matching sample IDs; order may differ.

| File | First column | Other columns |
| --- | --- | --- |
| `taxon-abundance.tsv` | `taxon_id` | Nonnegative finite numeric abundances per sample |
| `functional-features.tsv` | `feature_id` | Nonnegative finite numeric abundances per sample |
| `sample-metadata.tsv` | `sample_id` | Required `group`; optional study covariates |

Each sample must have a positive total in both abundance tables. Record profiling tool, database/version, taxonomic rank, functional definitions, and units before using real data. Do not mix taxonomic ranks or incompatible units. Functional values are normalized within each sample for the exploratory plot; raw values are retained in the input. A relative functional profile does not measure pathway activity or absolute abundance.

For your own data, put the three files in `data/my-study/`, document provenance in its README, then substitute `data/my-study` and `results/my-study` in the check and run commands. Add generated figures to `_quarto.yml` resources when extending the guide. For repeated samples, account for participant IDs and timepoints before inference.

## Scope and limitations

This workflow starts from processed profiles, not raw FASTQ files. It provides descriptive analysis only: no significance tests, differential abundance, causal inference, or clinical conclusions. Shannon describes the supplied taxonomic profiles and depends on resolution, filtering, and profiling method. Functional features indicate inferred functional potential, not measured activity. Synthetic group differences are constructed teaching examples.

## Folder guide

- `examples/`: seven learning steps for each dataset.
- `data/`: synthetic inputs and preparation notes for other datasets.
- `scripts/python/`: input checking, analysis, and independent sample verification.
- `results/`: generated figures and tables, separated by dataset.
- `reports/`: learner report template.
- `assets/css/`: guide styling.

## Troubleshooting

`ModuleNotFoundError`: activate `.venv` and rerun installation. File not found: check that the terminal is in the main project folder. Input check failure: correct the reported IDs, missing values, or numeric values before running again. `quarto: command not found`: install Quarto separately; Python analysis remains available. To finish a session, run `deactivate`.

Dependencies are intentionally unpinned in this starter; each run records installed versions. Review and lock versions before using this for a reproducible research release.
