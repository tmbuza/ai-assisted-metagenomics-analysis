"""Validate the same feature-row TSV contract as the main workflow."""
import argparse
from pathlib import Path
from run_workflow import load_inputs


def verify_inputs(folder):
    meta, taxa, functions = load_inputs(Path(folder))
    return meta.group.to_dict(), taxa.T.to_dict(orient="index"), functions.T.to_dict(orient="index")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    args = parser.parse_args()
    try:
        meta, _, _ = load_inputs(args.inputs)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"Input check failed: {error}\n")
    print(f"Verified {len(meta)} matched samples. Groups: {meta.group.value_counts().to_dict()}")
