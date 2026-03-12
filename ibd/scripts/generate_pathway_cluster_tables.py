#!/usr/bin/env python3
"""Generate pathway cluster tables for UC (main), uninflamed and CD (supplementary).

Usage:
    python generate_pathway_cluster_tables.py \
        --uc-clusters results/uc_vs_hc/hierarchical_pathway_clusters.csv \
        --uninfl-clusters results/uninf_vs_hc/hierarchical_pathway_clusters.csv \
        --cd-clusters results/cd_vs_hc/hierarchical_pathway_clusters.csv \
        --output-dir tables/
"""
import argparse, os
import pandas as pd


def clean_name(term):
    """Remove Reactome ID from pathway name."""
    return term.split(" R-HSA")[0]


def make_cluster_table(cl_df, caption, label, top_pathway_col="Representative"):
    """Generate a LaTeX table from a hierarchical_pathway_clusters.csv."""
    rows = []
    for _, r in cl_df.iterrows():
        name = clean_name(r[top_pathway_col])
        # Truncate very long names
        if len(name) > 60:
            name = name[:57] + "..."
        size = int(r["Size"])
        nes = r["Mean_NES"]
        direction = "\\Up" if nes > 0 else "\\Down"
        rows.append(f"  {name} & {size} & {nes:+.2f} \\\\")

    return (
        "\\begin{table}[t]\n"
        "\\centering\n"
        "\\scriptsize\n"
        "\\renewcommand{\\arraystretch}{1.15}\n"
        "\\begin{tabular}{lrr}\n"
        "\\hline\n"
        "Cluster representative & Pathways & Mean NES \\\\\n"
        "\\hline\n"
        + "\n".join(rows) + "\n"
        "\\hline\n"
        "\\end{tabular}\n"
        f"\\caption{{{caption}}}\n"
        f"\\label{{{label}}}\n"
        "\\end{table}\n"
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--uc-clusters", required=True)
    p.add_argument("--uninfl-clusters", required=True)
    p.add_argument("--cd-clusters", required=True)
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    uc = pd.read_csv(args.uc_clusters, encoding="latin1")
    un = pd.read_csv(args.uninfl_clusters, encoding="latin1")
    cd = pd.read_csv(args.cd_clusters, encoding="latin1")

    # Sort by |Mean_NES| descending within direction (up first, then down)
    for df in [uc, un, cd]:
        df["_dir"] = (df["Mean_NES"] > 0).astype(int)
        df["_abs"] = df["Mean_NES"].abs()
        df.sort_values(["_dir", "_abs"], ascending=[False, False], inplace=True)
        df.drop(columns=["_dir", "_abs"], inplace=True)

    # UC main table
    tex_uc = make_cluster_table(
        uc,
        "Pathway clusters in inflamed UC. Each row is a community of co-enriched "
        "Reactome pathways identified by Louvain clustering on shared leading-edge "
        "genes. Pathways: number of member pathways; Mean NES: average normalised "
        "enrichment score across member pathways (positive = upregulated in UC).",
        "tab:pathway_clusters_uc",
    )
    with open(os.path.join(args.output_dir, "pathway_clusters_uc.tex"), "w") as f:
        f.write(tex_uc)
    print(f"UC: {len(uc)} clusters -> {args.output_dir}/pathway_clusters_uc.tex")

    # Uninflamed supplementary table
    tex_un = make_cluster_table(
        un,
        "Pathway clusters in uninflamed UC mucosa versus healthy controls.",
        "tab:pathway_clusters_uninfl",
    )
    with open(os.path.join(args.output_dir, "pathway_clusters_uninfl.tex"), "w") as f:
        f.write(tex_un)
    print(f"Uninfl: {len(un)} clusters -> {args.output_dir}/pathway_clusters_uninfl.tex")

    # CD supplementary table
    tex_cd = make_cluster_table(
        cd,
        "Pathway clusters in inflamed Crohn's disease mucosa versus healthy controls.",
        "tab:pathway_clusters_cd",
    )
    with open(os.path.join(args.output_dir, "pathway_clusters_cd.tex"), "w") as f:
        f.write(tex_cd)
    print(f"CD: {len(cd)} clusters -> {args.output_dir}/pathway_clusters_cd.tex")


if __name__ == "__main__":
    main()
