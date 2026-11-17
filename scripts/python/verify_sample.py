"""Independently verify a saved sample using the Python standard library."""
import csv, math, sys
from pathlib import Path

def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: python scripts/python/verify_sample.py DATA_DIR RESULTS_DIR SAMPLE_ID")
    data, results, sample = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    with (data / "taxon-abundance.tsv").open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        if sample not in reader.fieldnames:
            raise SystemExit(f"Sample not found: {sample}")
        values = [float(row[sample]) for row in reader]
    total = sum(values)
    if total <= 0 or any(not math.isfinite(x) or x < 0 for x in values):
        raise SystemExit("Invalid abundance values")
    expected = {"total_abundance": total, "detected_taxa": sum(x > 0 for x in values), "shannon": -sum((x/total)*math.log(x/total) for x in values if x > 0)}
    with (results / "tables/alpha-diversity.tsv").open() as f:
        saved = next((row for row in csv.DictReader(f, delimiter="\t") if row["sample_id"] == sample), None)
    if saved is None:
        raise SystemExit("Sample missing from saved results")
    for key, value in expected.items():
        actual = float(saved[key])
        print(f"{key}: recalculated={value:.9f}; saved={actual:.9f}")
        if not math.isclose(value, actual, rel_tol=1e-9, abs_tol=1e-9):
            raise SystemExit(f"Verification failed: {key}")
    print("Sample verification passed.")

if __name__ == "__main__":
    main()
