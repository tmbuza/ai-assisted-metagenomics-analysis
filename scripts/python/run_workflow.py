"""Descriptive analysis of processed metagenomic profiles."""
import argparse, csv, json, platform, hashlib, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_table(path, key):
    with path.open(newline="") as f:
        header = next(csv.reader(f, delimiter="\t"), [])
    if not header or header[0] != key or len(header) != len(set(header)):
        raise ValueError(f"{path}: expected unique headers beginning with {key}")
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if t.empty or (t == "").any().any() or t[key].duplicated().any():
        raise ValueError(f"{path}: empty table, blank values, or duplicate IDs")
    return t.set_index(key)

def load_inputs(folder):
    meta = read_table(folder / "sample-metadata.tsv", "sample_id")
    if "group" not in meta or meta["group"].str.strip().eq("").any():
        raise ValueError("Metadata requires a nonempty group column")
    tables = []
    for filename, key in [("taxon-abundance.tsv", "taxon_id"), ("functional-features.tsv", "feature_id")]:
        t = read_table(folder / filename, key).apply(pd.to_numeric, errors="raise")
        if set(t.columns) != set(meta.index):
            raise ValueError(f"{filename}: sample IDs must match metadata exactly")
        t = t.loc[:, meta.index]
        if not np.isfinite(t.to_numpy()).all() or (t < 0).any().any() or (t.sum() <= 0).any():
            raise ValueError(f"{filename}: values must be finite, nonnegative, with positive sample totals")
        tables.append(t)
    return meta, *tables

def feature_labels(folder, filename, key, label):
    path = folder / filename
    if not path.exists():
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if key not in frame or label not in frame or frame[key].duplicated().any():
        raise ValueError(f"Invalid annotations: {path}")
    return dict(zip(frame[key], frame[label]))


def plot_stacked(values, meta, labels, name, title, figures, tables):
    # Ranking uses mean sample proportion; do not select by group differences.
    top = values.mean(axis=1).nlargest(min(15, len(values))).index
    display = values.loc[top].copy()
    remainder = values.drop(top)
    if not remainder.empty:
        display.loc["Other (remaining features)"] = remainder.sum()
    display.to_csv(tables / f"{name}-display.tsv", sep="\t")
    colors = list(plt.get_cmap("tab20").colors[:len(top)])
    if not remainder.empty:
        colors.append("#d9d9d9")
    fig, ax = plt.subplots(figsize=(15, 6))
    x = np.arange(len(meta)); bottom = np.zeros(len(meta))
    for (feature, row), color in zip(display.iterrows(), colors):
        label = labels.get(feature, feature)
        ax.bar(x, row.to_numpy(), bottom=bottom, width=.95, color=color,
               label=textwrap.fill(label, 39), linewidth=0)
        bottom += row.to_numpy()
    # Stable grouping preserves metadata order within each group.
    for group, frame in meta.groupby("group", sort=False):
        positions = meta.index.get_indexer(frame.index)
        left, right = positions.min(), positions.max()
        if left:
            ax.axvline(left-.5, color="white", linewidth=3)
        ax.text((left+right)/2, 1.02, f"{group} (n={len(frame)})",
                transform=ax.get_xaxis_transform(), ha="center", fontsize=10)
    ax.set(xlim=(-.5, len(meta)-.5), ylim=(0, 1), xlabel="Sample (grouped by metadata label)",
           ylabel="Within-sample proportion", title=title)
    ax.set_title(title, pad=38)
    ax.set_xticks([])
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(figures / f"{name}.png", dpi=160); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--check-only", action="store_true")
    a = p.parse_args()
    try:
        meta, taxa, functions = load_inputs(a.data_dir)
    except (ValueError, OSError, KeyError) as e:
        p.exit(1, f"Input check failed: {e}\n")
    print(f"Input checks passed: {len(meta)} samples; {meta.group.nunique()} groups.")
    if meta.group.nunique() == 1:
        print("One group detected. Continuing with descriptive analysis.")
    if a.check_only:
        return
    order = [sample for _, frame in meta.groupby("group", sort=False) for sample in frame.index]
    meta = meta.loc[order]
    taxa, functions = taxa.loc[:, order], functions.loc[:, order]
    figures, tables = a.output_dir / "figures", a.output_dir / "tables"
    figures.mkdir(parents=True, exist_ok=True); tables.mkdir(parents=True, exist_ok=True)
    rel = taxa.div(taxa.sum(), axis=1)
    shannon = -(rel * np.log(rel.where(rel > 0, 1))).sum()
    alpha = pd.DataFrame({"total_abundance": taxa.sum(), "detected_taxa": (taxa > 0).sum(), "shannon": shannon}).join(meta)
    alpha.index.name = "sample_id"
    alpha.to_csv(tables / "alpha-diversity.tsv", sep="\t")
    rel.to_csv(tables / "taxon-relative-abundance.tsv", sep="\t")
    frel = functions.div(functions.sum(), axis=1)
    frel.to_csv(tables / "functional-relative-abundance.tsv", sep="\t")
    tax_labels = feature_labels(a.data_dir, "taxon-annotations.tsv", "taxon_id", "taxon_label")
    function_labels = feature_labels(a.data_dir, "functional-annotations.tsv", "feature_id", "feature_label")
    plot_stacked(rel, meta, tax_labels, "composition", "Top 15 taxa by mean sample proportion; remainder = Other",
                 figures, tables)
    plot_stacked(frel, meta, {}, "functional-features", "Top 15 functional features within the supplied feature set",
                 figures, tables)
    summary = alpha.groupby("group", sort=False).agg(
        samples=("shannon", "size"), median_shannon=("shannon", "median"),
        minimum_shannon=("shannon", "min"), maximum_shannon=("shannon", "max"),
        median_detected_taxa=("detected_taxa", "median"),
        minimum_detected_taxa=("detected_taxa", "min"), maximum_detected_taxa=("detected_taxa", "max"))
    summary.to_csv(tables / "group-summary.tsv", sep="\t")
    # For verified LiJ_2017 data, retain original source-scale pathway values.
    info_path = a.data_dir / "measurement-info.json"
    info = json.loads(info_path.read_text()) if info_path.exists() else None
    if info and info.get("study") == "LiJ_2017" and function_labels:
        selected = functions.mean(axis=1).nlargest(min(15, len(functions))).index
        group_medians = functions.loc[selected].T.join(meta[["group"]]).groupby("group", sort=False).median().T
        group_medians.to_csv(tables / "pathway-source-group-medians.tsv", sep="\t")
        fig, ax = plt.subplots(figsize=(13, 8))
        image = ax.imshow(group_medians.to_numpy(), aspect="auto", cmap="viridis", vmin=0)
        ax.set_yticks(range(len(selected)), [textwrap.fill(function_labels.get(i, i)+" ["+i+"]", 60) for i in selected], fontsize=8)
        ax.set_xticks(range(len(group_medians.columns)), group_medians.columns)
        ax.set_title("Median source pathway relative abundance\nTop 15 by mean source abundance across samples", pad=12)
        fig.colorbar(image, ax=ax, label="HUMAnN source relative abundance (not activity)")
        fig.tight_layout(); fig.savefig(figures / "pathway-source-heatmap.png", dpi=160); plt.close(fig)
    # Plot every sample, with a short line marking each group's median.
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = list(alpha.groupby("group", sort=False))
    colors = plt.get_cmap("tab10")

    for i, (group, frame) in enumerate(groups):
        # Horizontal spacing separates dots without changing diversity values.
        offsets = np.linspace(-0.12, 0.12, len(frame))
        x_positions = i + offsets

        ax.scatter(
            x_positions,
            frame["shannon"],
            color=colors(i % 10),
            s=65,
            zorder=3,
        )

        median = frame["shannon"].median()
        ax.hlines(
            median,
            i - 0.25,
            i + 0.25,
            color="black",
            linewidth=2,
            zorder=2,
        )

    ax.set_xticks(
        range(len(groups)),
        [f"Group {group}\n(n = {len(frame)})" for group, frame in groups],
    )
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.margins(y=0.15)
    ax.set_ylabel("Shannon diversity (natural log)")
    ax.set_title(
        "Shannon diversity by group\n"
        "Each dot is one sample; black lines show group medians",
        pad=12,
    )
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    fig.text(
        0.5, 0.02,
        "Higher values indicate greater diversity. "
        "Horizontal dot positions only separate samples.",
        ha="center",
        fontsize=9,
    )

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(figures / "shannon-diversity.png", dpi=160)
    plt.close(fig)

    tracked = [a.data_dir / n for n in ["sample-metadata.tsv", "taxon-abundance.tsv", "functional-features.tsv",
               "taxon-annotations.tsv", "functional-annotations.tsv", "measurement-info.json", "functional-denominators.tsv"]]
    checksums = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in tracked if f.exists()}
    record = {"input_sha256": checksums, "samples": len(meta), "groups": meta.group.value_counts().to_dict(),
              "taxa": len(taxa), "functional_features": len(functions), "display_top_n": 15,
              "measurement_info": info,
              "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__, "data_dir": str(a.data_dir), "analysis": "descriptive; natural-log Shannon; within-sample normalization"}
    (tables / "run-record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Results saved to {a.output_dir}")

if __name__ == "__main__":
    main()
