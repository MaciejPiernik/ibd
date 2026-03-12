#!/usr/bin/env python3
"""Filter a per-dataset gene stats file to keep only genes present in a meta-analysis summary file."""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Filter first CSV to genes present in second CSV."
    )
    parser.add_argument("per_dataset", help="Per-dataset stats file (gene + dataset rows)")
    parser.add_argument("meta", help="Meta-analysis summary file (one row per gene)")
    parser.add_argument("-o", "--output", default="filtered_output.tsv",
                        help="Output filename (default: filtered_output.tsv)")
    args = parser.parse_args()

    df = pd.read_csv(args.per_dataset)
    meta = pd.read_csv(args.meta)

    genes_to_keep = set(meta["gene"])
    filtered = df[df["gene"].isin(genes_to_keep)]

    print(f"Genes in meta file:      {len(genes_to_keep)}")
    print(f"Rows before filtering:   {len(df)}")
    print(f"Rows after filtering:    {len(filtered)}")
    print(f"Unique genes kept:       {filtered['gene'].nunique()}")

    filtered.to_csv(args.output, index=False)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
