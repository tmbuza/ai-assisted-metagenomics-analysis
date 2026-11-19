# AI-Assisted Metagenomics Analysis

A free CDI learning project for exploring processed taxonomic and functional profiles with Python. Learn to check inputs, run an analysis, verify calculations and explain results with evidence and limitations.

Use AI to help understand code or draft an explanation, then check the calculations and claims yourself.

## Start here

**New to the project? Follow the complete Required setup page in the guide.** It covers software, the project folder, Python environments, package installation, terminal use and troubleshooting. The source is [`00-setup.qmd`](00-setup.qmd).

To open the guide locally after installing Quarto, run this from the main project folder:

```bash
quarto preview
```

Open the address printed in the terminal. Press **Ctrl+C** to stop preview. Quarto is a separate application; it is not installed by pip. You can complete the Python analysis without it.

The QMD links in this README point to lesson source files in the repository. Use the rendered guide for the formatted learning experience. Learner-facing guide pages explain their requirements directly rather than sending learners back to this README.

## The learning route

| Example | Purpose | Current data status |
|---|---|---|
| [Synthetic practice](examples/synthetic/00-start-here.qmd) | Learn the complete workflow, inspect Python calculations and verify S01 | Fictional data; functional labels are invented identifiers, not real pathways |
| [Hypertension example](examples/hypertension/00-start-here.qmd) | Repeat the workflow with 195 samples: 98 HTN, 56 pHTN and 41 Control | Taxonomic abundances are synthetic placeholders; functional source, definitions and units still require verification |
| [Your own study](examples/my-study/00-start-here.qmd) | Prepare your tables, verify a sample and write your own report | Learner-supplied data, with documented provenance and units |

Each example follows the same sequence: define the question, preview results, meet the inputs, validate, analyse, verify and explain, then write the report. Exact commands appear in the corresponding lesson pages.

## Quick start for terminal users

Download and extract the complete project. Open `ai-assisted-metagenomics-analysis` in your editor and open a terminal there. You should see `_quarto.yml`, `requirements.txt`, `scripts/` and `data/`. No Git clone is required for a downloaded copy.

Use Python 3.10 or later. Create an environment once, then activate it whenever you open a new terminal.

### macOS or Linux

```bash
python3 --version
python3 -m venv .venv
source .venv/bin/activate
```

### Windows PowerShell

```powershell
py --version
py -m venv .venv
.\.venv\Scripts\Activate.ps1
```

If PowerShell blocks activation, open Command Prompt in the project folder and use `.venv\Scripts\activate.bat`. The setup page explains alternatives if `py` is unavailable.

### Install, check, run and verify

With the environment active:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -c "import numpy, pandas, matplotlib; print('Analysis packages loaded successfully.')"
python scripts/python/run_workflow.py --data-dir data/synthetic --output-dir results/synthetic --check-only
python scripts/python/run_workflow.py --data-dir data/synthetic --output-dir results/synthetic
python scripts/python/verify_sample.py data/synthetic results/synthetic S01
```

For the bundled synthetic version, validation reports `Input checks passed: 12 samples; 2 groups.` The independently checked S01 reference values are total abundance **546**, **8** detected taxa and natural-log Shannon diversity **2.012363**. If inputs change, recalculate these values rather than forcing agreement.

`--check-only` does not refresh outputs. The full run creates or replaces named results in `results/synthetic/`. Verification recalculates one sample and compares it with the saved diversity table. The lesson also shows Python code for independent checks of taxon and functional proportions.

The longer lesson blocks using `python - <<'PY'` are for Bash/Zsh. In PowerShell or Command Prompt, save only the enclosed Python lines as a `.py` file and run it with `python filename.py`; see the setup page for the full explanation.

## Inputs and outputs

All inputs are TSV files with unique headers and row identifiers. Metadata has samples in rows; both abundance tables have samples in columns.

| Filename | First header | Remaining columns |
|---|---|---|
| `sample-metadata.tsv` | `sample_id` | Required `group`; optional sample annotations |
| `taxon-abundance.tsv` | `taxon_id` | Sample IDs with taxonomic abundance values |
| `functional-features.tsv` | `feature_id` | Sample IDs with functional feature values |

Sample IDs must match exactly across the three files; order may differ. Tables must have no blank cells. Abundances must be finite, numeric and nonnegative, with a positive total per sample in each abundance table. One descriptive group is supported. The current workflow requires both abundance tables; do not invent functional values to fill a missing file.

Each run writes three figures to `results/<dataset>/figures/`: composition, Shannon diversity and functional profile. It also writes diversity and relative-abundance tables plus `run-record.json` to `results/<dataset>/tables/`.

## Reports and rendering

Each example's final page includes a report template and correct figure paths. Complete the prompts with your checked values. For synthetic practice:

```bash
quarto render reports/synthetic-report.qmd --to html
```

To build the full guide:

```bash
quarto render
```

The guide build writes to `docs/`. Quarto displays existing results; it does not execute the Python analysis in these lesson pages. After changing inputs, rerun Python before rendering the updated figures.

## Scope and reproducibility

This workflow begins with processed profiles, not raw sequencing reads. It is descriptive: it does not test group differences, account for repeated samples, adjust for confounders or establish causation.

Taxonomic and functional values are normalized separately within samples. Confirm the units and feature definitions before deciding whether this is biologically appropriate. Relative profiles do not directly establish absolute abundance or functional activity. Hypertension placeholder composition and Shannon values cannot support conclusions about the source study.

Keep original inputs, preparation/source notes, scripts, outputs and the report together. Record the input version and analysis date. Dependencies are unpinned in this starter; the run record captures installed package versions but does not recreate the environment or record input checksums by itself.

## Folder guide

| Folder | Contents |
|---|---|
| `examples/` | Seven lesson pages per dataset |
| `data/` | Input tables and provenance notes |
| `scripts/python/` | Analysis and verification scripts |
| `results/` | Figures, tables and run records by dataset |
| `reports/` | QMD reports |
| `assets/css/` | Guide styling |
| `docs/` | Rendered guide |

## Troubleshooting

- **File not found:** confirm you are in the main project folder.
- **`ModuleNotFoundError`:** activate `.venv` and reinstall requirements using `python -m pip`.
- **Input check failed:** inspect the named file and follow the example's input-check page.
- **`quarto` not found:** install Quarto separately and reopen the terminal.
- **Results disagree:** check the sample ID, input version, output folder and whether the analysis was rerun.

For a detailed walkthrough, use [Required setup](00-setup.qmd). For definitions, formulas and interpretation guidance, use the [Appendix](99-appendix.qmd). To finish a Python session, run `deactivate` in the terminal.
