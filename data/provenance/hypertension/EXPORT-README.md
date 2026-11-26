# Verified LiJ_2017 inputs for the Python metagenomics guide

Use `data/hypertension/` for the main Control versus HTN worked example. The optional `data/hypertension-all-groups/` extension retains pre-hypertension.

| Dataset | Control | pHTN | HTN | Total |
|---|---:|---:|---:|---:|
| hypertension | 41 | 0 | 99 | 140 |
| hypertension-all-groups | 41 | 56 | 99 | 196 |

Both contain 440 species and 483 named MetaCyc pathways. Each sample represents a distinct subject. These are real, processed shotgun metagenomic profiles, replacing the provisional placeholder tables.

## Install and run

Copy this package's `data/` and `scripts/` folders into the existing guide, merging folders and replacing its provisional `data/hypertension/` files. Keep your guide's existing Python analysis scripts. Run from the guide root:

```bash
python scripts/python/run_workflow.py --data-dir data/hypertension --output-dir results/hypertension
python scripts/python/verify_sample.py data/hypertension results/hypertension ERR1398068
```

For the optional three-group extension, use `data/hypertension-all-groups` and `results/hypertension-all-groups` instead. Retire previous hypertension outputs before presenting newly generated outputs. This bundle replaces inputs; it does not update the guide pages or introduce differential abundance analysis.

## What the values mean

- `taxon-abundance.tsv`: species percentages from curatedMetagenomicData, approximately 100 per sample. These are NOT read counts or absolute bacterial abundances.
- `functional-features.tsv`: unstratified named pathway relative abundances from HUMAnN, on their source scale. Taxon-stratified contribution rows are excluded to avoid double counting. `UNMAPPED` and `UNINTEGRATED` are excluded from this feature table and retained in `functional-denominators.tsv`.
- The existing Python workflow divides each functional column by its named-pathway total. Its functional plots therefore show proportions within the retained named pathways. They do not show the fraction of all microbial functions or pathway activity. Pathways also overlap biologically.
- `number_reads` is source sample metadata, not the sum of either abundance table. Do not use percentage totals as sequencing depth. Do not feed either table to methods requiring integer read counts, and do not manufacture counts by multiplying percentages by read depth.
- Annotation tables map compact IDs to labels and complete original feature IDs. No abundance filtering or artificial pseudocounts were applied. The main subset includes 31 all-zero species and five all-zero pathways; filter explicitly if a downstream method requires it.

## Metadata and study interpretation

Groups come from `study_condition`, which distinguishes control, pre-hypertension and hypertension. The broader `disease` field combines pre-hypertension and hypertension and must not define this comparison. Analysis sample IDs are ENA run accessions; `source_sample_id` preserves the original study identifier.

The previous input set had 195 samples. The source contains 196, including hypertension sample ERR1398137 (H2M413924), which was missing before. The other 195 accession-based group labels agree with the replacement metadata. Labels were also checked against the original paper's supplementary Table S2.

Age, sex and BMI are unavailable in the pinned curated metadata snapshot. They were not invented or filled from sample names. This is an observational association example: group differences do not establish that microbes cause hypertension. Any unadjusted analysis needs that limitation in its report.

## Reproduce the export

The analysis remains Python. Maintainers can regenerate the TSVs using only base R, from the included source matrices and scoped metadata:

```bash
Rscript scripts/R/prepare_hypertension.R
```

To retrieve the exact original resources again (Python standard library only):

```bash
python scripts/python/fetch_hypertension_sources.py
Rscript scripts/R/prepare_hypertension.R
```

The downloader verifies pinned SHA-256 checksums and refuses mismatched existing files. To regenerate the scoped metadata from the downloaded full snapshot, remove only `data/provenance/hypertension/LiJ_2017-curated-metadata.tsv` before running the R preparation script. The full all-study metadata file is deliberately not bundled.

`data/provenance/hypertension/` records ExperimentHub resource identifiers, download URLs, a pinned metadata commit, checksums, sample checks and group summaries. Profiles are the 2021-10-14 curatedMetagenomicData releases, not a rerun of raw sequencing. Exact historical profiler software versions are not embedded in these resource matrices.

## Validation

The existing Python workflow completed for both datasets; its sample verifier passed. Independent calculations reproduced every saved species Shannon value to absolute tolerance 1e-12. Identifier matching, finite nonnegative values, source totals, group counts, distinct subjects and annotations were checked. Re-exporting from the original metadata snapshot reproduces all twelve exported TSV files byte for byte. See `data/provenance/hypertension/validation/`.

## Sources

- Li et al. (2017), *Gut microbiota dysbiosis contributes to the development of hypertension*: https://doi.org/10.1186/s40168-016-0222-x ; sequencing study PRJEB13870.
- curatedMetagenomicData: https://waldronlab.io/curatedMetagenomicData/ ; study `LiJ_2017`.
- ExperimentHub EH7181, species-containing relative abundance profiles; EH7179, pathway abundance profiles. Exact download URLs and SHA-256 hashes are in `download-manifest.json`.

Cite the original study, curatedMetagenomicData and the profiling tools when using these data; consult their published documentation and data-use terms for your intended use.
