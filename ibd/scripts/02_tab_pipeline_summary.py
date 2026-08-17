#!/usr/bin/env python3
"""Generate the summary table covering every parallel meta-analysis.

The analyses form the complete set of pairwise comparisons between the three UC
sample groups (inflamed, uninflamed, control), together with the two involving
Crohn's disease.

Usage:
    python ibd/scripts/02_tab_pipeline_summary.py \
        --results results/revision_v2 \
        --output-tex results/revision_v2/pipeline_summary.tex
"""
import argparse, os
import pandas as pd

# column label, directory, sample count column for each arm
ANALYSES = [
    ("UC inf.", "uc_vs_hc", "n_uc", "n_ctrl"),
    ("UC uninf.", "uninf_vs_hc", "n_uc", "n_ctrl"),
    ("Inf. vs uninf.", "infl_vs_uninfl", "n_uc_inflamed", "n_uc_uninflamed"),
    ("CD inf.", "cd_vs_hc", "n_cd", "n_ctrl"),
    ("UC vs CD", "uc_vs_cd", "n_uc", "n_cd"),
]

ROWS = [
    ("Datasets", "datasets"),
    ("First group samples", "n_disease"),
    ("Second group samples", "n_ctrl"),
    ("Genes tested", "genes_tested"),
    ("Significant genes ($q < 0.05$)", "genes_sig"),
    ("\\quad Upregulated", "genes_up"),
    ("\\quad Downregulated", "genes_down"),
    ("Median $|g|$", "median_g"),
    ("Pathways tested", "pw_tested"),
    ("Pathways significant", "pw_sig"),
    ("Pathways after knee filter", "pw_knee"),
    ("\\quad Upregulated", "pw_knee_up"),
    ("\\quad Downregulated", "pw_knee_down"),
    ("Pathway clusters", "clusters"),
]


def load_analysis(d, disease_col, ctrl_col):
    ds = pd.read_csv(os.path.join(d, "dataset_summary.csv"))
    meta = pd.read_csv(os.path.join(d, "gene_meta_all.csv"), encoding="latin1")

    min_datasets = meta['datasets'].max() / 2
    sig = meta[
        (meta['datasets'] >= min_datasets)
        & (meta['direction_ratio'] >= 0.8)
        & (meta['q_value'] < 0.05)
    ]
    pw = pd.read_csv(os.path.join(d, "pathway_meta_all.csv"), encoding="latin1")

    pw_sig = pw[
        (pw['datasets'] >= min_datasets)
        & (pw['direction_ratio'] >= 0.8)
        & (pw['q_value'] < 0.05)
    ]

    knee_path = os.path.join(d, "pathway_meta_significant_knee_with_clusters.csv")
    pw_knee = pd.read_csv(knee_path, encoding="latin1") if os.path.exists(knee_path) else None

    cl_path = os.path.join(d, "hierarchical_pathway_clusters.csv")
    clusters = pd.read_csv(cl_path, encoding="latin1") if os.path.exists(cl_path) else None

    return {
        "datasets": len(ds),
        "n_disease": int(ds[disease_col].sum()),
        "n_ctrl": int(ds[ctrl_col].sum()),
        "genes_tested": len(meta),
        "genes_sig": len(sig),
        "genes_up": int((sig["direction"] == "up").sum()),
        "genes_down": int((sig["direction"] == "down").sum()),
        "median_g": sig["hedges_g"].abs().median(),
        "pw_tested": len(pw),
        "pw_sig": len(pw_sig),
        "pw_knee": len(pw_knee) if pw_knee is not None else "---",
        "pw_knee_up": int((pw_knee["direction"] == "up").sum()) if pw_knee is not None else "---",
        "pw_knee_down": int((pw_knee["direction"] == "down").sum()) if pw_knee is not None else "---",
        "clusters": len(clusters) if clusters is not None else "---",
    }


def fmt(val):
    if isinstance(val, float):
        return f"{val:.2f}"
    if isinstance(val, int):
        return f"{val:,}".replace(",", "{,}")
    return str(val)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", required=True)
    p.add_argument("--output-tex", required=True)
    args = p.parse_args()

    loaded = [(label, load_analysis(os.path.join(args.results, d), disease, ctrl))
              for label, d, disease, ctrl in ANALYSES]

    lines = [
        f"  {label} & " + " & ".join(fmt(a[key]) for _, a in loaded) + " \\\\"
        for label, key in ROWS
    ]

    table = (
        "\\begin{table}[t]\n"
        "\\centering\n"
        "\\scriptsize\n"
        "\\renewcommand{\\arraystretch}{1.15}\n"
        "\\begin{tabular}{l" + "r" * len(loaded) + "}\n"
        "\\hline\n"
        " & " + " & ".join(label for label, _ in loaded) + " \\\\\n"
        "\\hline\n"
        + "\n".join(lines) + "\n"
        "\\hline\n"
        "\\end{tabular}\n"
        "\\caption{Summary of the parallel meta-analyses. The first three columns are the "
        "pairwise comparisons between inflamed UC, uninflamed UC, and control mucosa; the last "
        "two involve Crohn's disease. Upregulated and downregulated refer to the direction in "
        "the first group relative to the second. Significant genes: random-effects "
        "meta-analysis, $q < 0.05$, direction consistency $\\geq 80\\%$, presence in at least "
        "half the contributing datasets. Pathway filtering: knee-point on ranked "
        "$|\\mathrm{NES}|$; clusters: Louvain community detection on shared leading-edge "
        "genes.}\n"
        "\\label{tab:pipeline_summary}\n"
        "\\end{table}\n"
    )

    with open(args.output_tex, "w") as f:
        f.write(table)

    print(f"Written: {args.output_tex}")
    for label, a in loaded:
        print(f"  {label:15s} {a['datasets']:2d} datasets, {a['genes_sig']:5d} genes, "
              f"{a['pw_knee']} knee pathways, {a['clusters']} clusters")


if __name__ == "__main__":
    main()
