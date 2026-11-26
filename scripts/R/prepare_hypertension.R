# Reproduce prepared TSVs from the included curatedMetagenomicData sources.
# Maintainer preparation only; the learner's analysis remains Python.
# Requires only base R. Run from this package/project root:
# Rscript scripts/R/prepare_hypertension.R
# Optional: Rscript scripts/R/prepare_hypertension.R SOURCE_DIR OUTPUT_DATA_DIR
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(0L, 2L)) stop("Provide either no arguments or SOURCE_DIR OUTPUT_DATA_DIR.")
source_dir <- if (length(args)) args[1] else "data/provenance/hypertension"
data_dir <- if (length(args)) args[2] else "data"
options(digits = 17)
fail <- function(ok, text) if (!isTRUE(ok)) stop(text, call. = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE,
                                         row.names = FALSE, na = "NA")
read_matrix <- function(type) {
  title <- paste0("2021-10-14.LiJ_2017.", type)
  en <- new.env(parent = emptyenv())
  names <- load(file.path(source_dir, paste0(title, ".rda")), envir = en)
  fail(identical(names, title), paste("Unexpected source object:", type))
  a <- en[[title]]
  fail(is.matrix(a) && is.numeric(a) && all(is.finite(a)) && all(a >= 0),
       paste("Invalid source matrix:", type))
  fail(!anyDuplicated(rownames(a)) && !anyDuplicated(colnames(a)), "Duplicate source IDs.")
  a
}
metadata_path <- file.path(source_dir, "LiJ_2017-curated-metadata.tsv")
if (!file.exists(metadata_path)) {
  en <- new.env(parent = emptyenv())
  nm <- load(file.path(source_dir, "sampleMetadata.rda"), envir = en)
  fail("sampleMetadata" %in% nm, "Missing original metadata snapshot.")
  full <- en[["sampleMetadata"]]
  write_tsv(full[full$study_name == "LiJ_2017", , drop = FALSE], metadata_path)
}
m <- read.delim(metadata_path,
                check.names = FALSE, colClasses = "character", na.strings = "NA")
fail(nrow(m) == 196L && all(m$study_name == "LiJ_2017"), "Unexpected study metadata.")
for (field in c("sample_id", "subject_id", "NCBI_accession"))
  fail(!anyNA(m[[field]]) && all(nzchar(m[[field]])) && !anyDuplicated(m[[field]]),
       paste("Invalid metadata identifier:", field))
tax <- read_matrix("relative_abundance")
fun <- read_matrix("pathway_abundance")
fail(setequal(colnames(tax), m$sample_id) && setequal(colnames(fun), m$sample_id),
     "Source profiles and metadata do not contain exactly the same samples.")

# Species profiles only: including multiple ranks would count taxonomy repeatedly.
# Preserve species groups/complexes; no NCBI-name matching or rowTree filtering.
species <- grepl("(^|\\|)s__[^|]+$", rownames(tax))
tax <- tax[species, m$sample_id, drop = FALSE]
tax_ids <- sub("^.*\\|", "", rownames(tax))
fail(sum(species) == 440L && !anyDuplicated(tax_ids), "Unexpected species identifiers.")
fail(max(abs(colSums(tax) - 100)) < 0.001,
     "Species percentages do not sum approximately to 100.")
tax_annotation <- data.frame(taxon_id = tax_ids,
    taxon_label = gsub("_", " ", sub("^s__", "", tax_ids)),
    source_feature_id = rownames(tax), stringsAsFactors = FALSE)

# HUMAnN contains totals plus taxon-stratified contributions. Keep only totals.
unstratified <- !grepl("|", rownames(fun), fixed = TRUE)
special <- rownames(fun) %in% c("UNMAPPED", "UNINTEGRATED")
named <- unstratified & !special
fail(all(c("UNMAPPED", "UNINTEGRATED") %in% rownames(fun)), "Missing special rows.")
fail(max(abs(colSums(fun[unstratified, , drop = FALSE]) - 1)) < 0.0001,
     "Unstratified HUMAnN relative-abundance totals are not approximately one.")
function_source_ids <- rownames(fun)[named]
function_ids <- sub(":.*$", "", function_source_ids)
fail(sum(named) == 483L && !anyDuplicated(function_ids), "Unexpected pathway identifiers.")
functional_annotation <- data.frame(feature_id = function_ids,
    feature_label = sub("^[^:]+:[[:space:]]*", "", function_source_ids),
    source_feature_id = function_source_ids, stringsAsFactors = FALSE)
f <- fun[named, m$sample_id, drop = FALSE]
fail(all(colSums(f) > 0), "A sample lacks positive named-pathway abundance.")

groups <- c(control = "Control", `pre-hypertension` = "pHTN", hypertension = "HTN")
group <- unname(groups[m$study_condition])
fail(!anyNA(group), "Unexpected clinical label; do not infer labels from sample names.")
fail(identical(as.integer(table(factor(group, levels = c("Control", "pHTN", "HTN")))),
               c(41L, 56L, 99L)), "Unexpected clinical group counts.")
metadata <- data.frame(sample_id = m$NCBI_accession, group = group,
    subject_id = m$subject_id, source_sample_id = m$sample_id,
    study_name = m$study_name, source_study_condition = m$study_condition,
    body_site = m$body_site, country = m$country, number_reads = m$number_reads,
    stringsAsFactors = FALSE)
fail(!anyNA(metadata) && all(nzchar(as.matrix(metadata))), "Missing selected metadata fields.")
annotation_mass <- data.frame(sample_id = metadata$sample_id,
    source_sample_id = metadata$source_sample_id, group = group,
    unmapped = as.numeric(fun["UNMAPPED", m$sample_id]),
    unintegrated = as.numeric(fun["UNINTEGRATED", m$sample_id]),
    named_pathway_total = colSums(f))

export <- function(folder, keep) {
  ix <- which(keep)
  ix <- ix[order(factor(group[ix], levels = c("Control", "pHTN", "HTN")), metadata$sample_id[ix])]
  out <- file.path(data_dir, folder)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  t <- tax[, ix, drop = FALSE]; colnames(t) <- metadata$sample_id[ix]
  ff <- f[, ix, drop = FALSE]; colnames(ff) <- metadata$sample_id[ix]
  write_tsv(data.frame(taxon_id = tax_ids, t, check.names = FALSE), file.path(out, "taxon-abundance.tsv"))
  write_tsv(data.frame(feature_id = function_ids, ff, check.names = FALSE), file.path(out, "functional-features.tsv"))
  write_tsv(metadata[ix, ], file.path(out, "sample-metadata.tsv"))
  write_tsv(tax_annotation, file.path(out, "taxon-annotations.tsv"))
  write_tsv(functional_annotation, file.path(out, "functional-annotations.tsv"))
  write_tsv(annotation_mass[ix, ], file.path(out, "functional-denominators.tsv"))
  cat(folder, ":", length(ix), "samples;", nrow(t), "species;", nrow(ff), "named pathways\n")
  # No rounding, pseudo-counts, rarefaction, zero-row filtering or renormalisation.
  # Percentages and pathway relative abundance remain in their source scales.
}
export("hypertension", group %in% c("Control", "HTN"))
export("hypertension-all-groups", rep(TRUE, nrow(m)))
