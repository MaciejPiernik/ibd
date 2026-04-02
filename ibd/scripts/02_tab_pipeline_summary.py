#!/usr/bin/env python3
"""Generate a three-way pipeline summary table (UC inflamed, UC uninflamed, CD).

Usage:
    python ibd/scripts/02_tab_pipeline_summary.py \
        --uc-dir results/paper/uc_vs_hc \
        --uninfl-dir results/paper/uninf_vs_hc \
        --cd-dir results/paper/cd_vs_hc \
        --uccd-dir results/paper/uc_vs_cd \
        --output-tex results/paper/pipeline_summary.tex
"""
import argparse, os
import pandas as pd


def load_analysis(d, disease_col="n_uc"):
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

    n_disease = ds[disease_col].sum()
    n_ctrl = ds["n_ctrl"].sum() if "n_ctrl" in ds.columns else ds["n_cd"].sum() # for direct UC vs CD, "n_cd" column is used instead of "n_ctrl"

    return {
        "datasets": len(ds),
        "n_disease": int(n_disease),
        "n_ctrl": int(n_ctrl),
        "genes_tested": len(meta),
        "genes_sig": len(sig),
        "genes_up": int((sig["direction"] == "up").sum()),
        "genes_down": int((sig["direction"] == "down").sum()),
        "median_g": sig["hedges_g"].abs().median(),
        "pw_tested": len(pd.read_csv(os.path.join(d, "pathway_meta_all.csv"), encoding="latin1")),
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
    p.add_argument("--uc-dir", required=True)
    p.add_argument("--uninfl-dir", required=True)
    p.add_argument("--cd-dir", required=True)
    p.add_argument("--uccd-dir", required=True)
    p.add_argument("--output-tex", required=True)
    args = p.parse_args()

    uc = load_analysis(args.uc_dir, "n_uc")
    un = load_analysis(args.uninfl_dir, "n_uc")
    cd = load_analysis(args.cd_dir, "n_cd")
    uccd = load_analysis(args.uccd_dir, "n_uc")

    rows = [
        ("Datasets", "datasets"),
        ("Disease samples", "n_disease"),
        ("Control samples", "n_ctrl"),
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

    lines = []
    for label, key in rows:
        v_uc = fmt(uc[key])
        v_un = fmt(un[key])
        v_cd = fmt(cd[key])
        v_uccd = fmt(uccd[key])
        lines.append(f"  {label} & {v_uc} & {v_un} & {v_cd} & {v_uccd} \\\\")

    table = (
        "\\begin{table}[t]\n"
        "\\centering\n"
        "\\scriptsize\n"
        "\\renewcommand{\\arraystretch}{1.15}\n"
        "\\begin{tabular}{lrrrr}\n"
        "\\hline\n"
        " & UC inf. & UC uninf. & CD inf. & UC vs CD \\\\\n"
        "\\hline\n"
        + "\n".join(lines) + "\n"
        "\\hline\n"
        "\\end{tabular}\n"
        "\\caption{Summary of the four parallel meta-analyses. "
        "Summary of the four parallel meta-analyses. For the three disease vs control comparisons, upregulated/downregulated refers to the direction in disease relative to healthy controls. For the direct UC\,vs\,CD comparison, upregulated means higher in UC and downregulated means lower in UC relative to CD. "
        "Significant genes: random-effects meta-analysis, $q < 0.05$, "
        "direction consistency $\\geq 80\\%$. "
        "Pathway filtering: knee-point on ranked $|\\mathrm{NES}|$; "
        "clusters: Louvain community detection on shared leading-edge genes.}\n"
        "\\label{tab:pipeline_summary}\n"
        "\\end{table}\n"
    )

    with open(args.output_tex, "w") as f:
        f.write(table)

    print(f"Written: {args.output_tex}")
    for label, d in [("UC", uc), ("Uninfl", un), ("CD", cd), ("UC vs CD", uccd)]:
        print(f"  {label}: {d['datasets']} datasets, {d['genes_sig']} genes, "
              f"{d['pw_knee']} knee pathways, {d['clusters']} clusters")


if __name__ == "__main__":
    main()
