"""Descriptive analysis of processed metagenomic profiles."""
import argparse, csv, json, platform
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
    for values, name, title in [(rel, "composition", "Taxonomic composition"), (frel, "functional-features", "Relative functional profile")]:
        top = values.mean(axis=1).nlargest(10).index
        display = values.loc[top].copy()
        if len(values) > len(top):
            display.loc["Other (remaining features)"] = values.drop(top).sum()
        ax = display.T.plot.bar(stacked=True, figsize=(10, 5), width=.85)
        ax.set_xticks([])
        ax.set(xlabel="Sample", ylabel="Within-sample proportion", title=title)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
        plt.tight_layout(); plt.savefig(figures / f"{name}.png", dpi=160); plt.close()
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

        # for x, (sample_id, row) in zip(x_positions, frame.iterrows()):
        #     ax.annotate(
        #         sample_id,
        #         (x, row["shannon"]),
        #         xytext=(0, 8),
        #         ha="center",
        #         textcoords="offset points",
        #         fontsize=8,
        #         zorder=5,
        #         bbox=dict(
        #             facecolor="white",
        #             edgecolor="none",
        #             alpha=1,
        #             pad=1,
        #         ),
        #     )

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
    ax.set_ylabel("Shannon diversity")
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

    record = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__, "data_dir": str(a.data_dir), "analysis": "descriptive; natural-log Shannon; within-sample normalization"}
    (tables / "run-record.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Results saved to {a.output_dir}")

if __name__ == "__main__":
    main()
