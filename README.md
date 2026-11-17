# AI-assisted metagenomics analysis

Explore microbial composition and functional pathway profiles using Python, with AI assistance and checks you can perform yourself.

Start with a synthetic dataset, apply the workflow to a hypertension case study, then repeat it with your own data.

## One core workflow for every dataset

Learn these seven steps once, then repeat them with each worked example and your own study. The checks and analysis methods depend on your data and study design.

1. **Define the question** you want to answer.
2. **Preview the results** you will produce.
3. **Meet the inputs** and understand their structure.
4. **Check the inputs** for completeness and consistency.
5. **Run the analysis** using the supplied Python command.
6. **Verify and explain** the patterns using figures and tables.
7. **Write the report**, supporting your findings with figures and numbers and explaining what the analysis cannot establish.

## Choose your starting point

### 1. Synthetic community comparison

**Question:** How do microbial composition and functional pathway profiles differ between two hypothetical groups?

Start here to understand the input files, run the Python workflow, and verify selected results. The data contain deliberately constructed patterns for learning; they do not represent a real biological study.

### 2. Hypertension case study

**Question:** How do microbial composition and functional pathway profiles vary across blood-pressure groups?

Apply the workflow to study data. Examine differences between groups, check the supporting numbers, and consider other factors that could explain the patterns.

The example must document the source and processing of both the taxonomic and functional profiles. Group differences alone do not establish that the microbiome causes hypertension.

### 3. Your study

Bring your own taxonomic profiles, functional profiles, and sample metadata.

Define a question your data can address, prepare the required files, and follow the same seven steps. Use the worked examples as references, adapting the analysis to your study design.

## Meet the inputs

Each example uses three tab-separated files:

| File | Contents |
|---|---|
| `sample-metadata.tsv` | Sample identifiers, comparison groups, and relevant study information |
| `taxon-abundance.tsv` | Taxonomic abundance profiles for each sample |
| `functional-features.tsv` | Functional feature or pathway abundance profiles for each sample |

Each dataset's README explains its table structure, abundance units, feature identifiers, and data source.

These examples begin with prepared abundance tables. Generating taxonomic and functional profiles from raw sequencing reads is a separate upstream step.

## Understand the outputs

Use the figures and tables to explore:

- Microbial composition across samples and groups.
- Within-sample diversity calculated from taxonomic profiles.
- Functional feature or pathway profiles across samples and groups.
- Selected numerical results you can independently verify.

Interpret each result according to its measurement units and how the input profiles were generated.

## Use AI, then check its work

Ask AI to help explain code, investigate errors, or draft an interpretation. Check its suggestions against your input files, saved results, and study documentation.

Keep a record of the commands you ran, the checks you completed, and the decisions that affected your analysis.