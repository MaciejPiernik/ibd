"""
Cross-comparison gene table for the Discussion section.

Merges gene-level meta-analysis results across four analyses and produces a
LaTeX table showing key genes grouped by functional category, with effect
sizes, significance, and the activity-layer classification from S_gene_layers.py.

Usage:
    python ibd/scripts/09_cross_comparison_table.py \
        --uc-dir  uc_vs_hc \
        --uninf-dir uninf_vs_hc \
        --cd-dir  cd_vs_hc \
        --uccd-dir uc_vs_cd \
        --output  cross_comparison_table.tex
"""

import argparse
import os
import numpy as np
import pandas as pd


# ── gene list organised by Discussion section ─────────────────────────────
GENE_GROUPS = {
    "Therapeutic targets --- anti-TNF / integrin / S1P": [
        "TNF", "TNFRSF1B", "ITGA4", "ITGB7", "MADCAM1",
        "S1PR1", "S1PR4",
    ],
    "Therapeutic targets --- IL-23 / JAK-STAT": [
        "IL23A", "IL12A", "IL12B", "JAK1", "JAK2", "JAK3", "STAT1",
    ],
    "Therapeutic targets --- inflammasome / complement / CXCL": [
        "CASP1", "IL1B", "CASP4", "CASP5", "GSDMD", "IL18",
        "CFB", "C2", "C3", "C4A", "CD55",
        "CXCL1", "CXCL2", "CXCL3", "CXCL5",
    ],
    "Metabolic axis --- mitochondrial regulators": [
        "PPARGC1A", "PPARGC1B", "ESRRA",
    ],
    "Metabolic axis --- energy metabolism genes": [
        "OPA1", "PDP2", "IDH3A", "ACSF2", "SLC22A5",
    ],
    "Ferroptosis / ROS / metallothionein": [
        "DUOX2", "DUOXA2", "ACSL4", "GPX4", "SLC40A1",
        "MT1E", "MT1F", "MT1G", "MT1M", "MT1X", "MT2A",
    ],
    "Barrier integrity": [
        "CLDN8", "AQP8", "OCLN", "SLC26A2",
    ],
    "Tryptophan / kynurenine pathway": [
        "IDO1", "KYNU", "HAAO", "QPRT",
    ],
    "Glucuronidation (selected UGTs)": [
        "UGT1A1", "UGT1A6", "UGT2A3", "UGT2B17",
    ],
}


VALIDATED = r"$^{\dagger}$"


def _fmt_g(val: float, q: float, replicated: bool = False) -> str:
    """Effect size with significance stars, daggered where RNA-seq agrees."""
    if np.isnan(val):
        return "---"
    stars = ""
    if q < 0.001:
        stars = "***"
    elif q < 0.01:
        stars = "**"
    elif q < 0.05:
        stars = "*"
    sign = "+" if val > 0 else ""
    return f"{sign}{val:.2f}{stars}{VALIDATED if replicated else ''}"


TICK = r"$\bullet$"


def _layers(gene: str, layers):
    """The two properties from S_gene_layers.py; a gene may carry both.

    Returned as a pair of cells so that the overlap reads as a subset of each
    layer rather than as a category of its own.
    """
    if gene not in layers.index:
        return "---", "---"

    row = layers.loc[gene]

    return (TICK if row["constitutive"] else "",
            TICK if row["inflammation_dependent"] else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master", default="results/revision_v2/gene_master.csv",
                    help="gene_master.csv from S_master_tables.py")
    ap.add_argument("--output", default="results/revision_v2/cross_comparison_table.tex")
    args = ap.parse_args()

    # direction and support are text columns with gaps, so read without chunk inference
    master = pd.read_csv(args.master, low_memory=False).set_index("gene")

    # Only the three UC contrasts have an RNA-seq counterpart to validate against.
    columns = [("uc", "infl", "infl"), ("uninf", "uninf", "uninf"),
               ("direct", "direct", "direct"), ("cd", "cd", None), ("uccd", "uccd", None)]

    rows = []
    for group_name, genes in GENE_GROUPS.items():
        rows.append({"_group": group_name})
        for gene in genes:
            r = {"gene": gene}
            row = master.loc[gene] if gene in master.index else None
            for label, suffix, validated in columns:
                if row is None:
                    r[f"g_{label}"], r[f"q_{label}"] = np.nan, np.nan
                    r[f"v_{label}"] = False
                    continue
                r[f"g_{label}"] = row[f"hedges_g_{suffix}"]
                r[f"q_{label}"] = row[f"q_value_{suffix}"]
                r[f"v_{label}"] = bool(row[f"replicated_{validated}"]) if validated else False

            r["constit"], r["infldep"] = _layers(gene, master)
            rows.append(r)

    # ── build LaTeX ───────────────────────────────────────────────────────
    lines = []
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    lines.append(r"\scriptsize")
    lines.append(r"\renewcommand{\arraystretch}{1.10}")
    lines.append(r"\begin{tabular}{l rrrrr cc}")
    lines.append(r"\toprule")
    lines.append(r"Gene & $g_\mathrm{UC\,infl}$ & $g_\mathrm{UC\,uninfl}$"
                 r" & $g_\mathrm{infl\,vs\,uninfl}$"
                 r" & $g_\mathrm{CD}$ & $g_\mathrm{UC\,vs\,CD}$"
                 r" & Constit & Infl-dep \\")
    lines.append(r"\midrule")

    for r in rows:
        if "_group" in r:
            lines.append(r"\midrule")
            lines.append(r"\multicolumn{8}{l}{\textit{"
                         + r["_group"] + r"}} \\")
            continue

        gene_tex = r["gene"].replace("_", r"\_")
        g_uc = _fmt_g(r["g_uc"], r["q_uc"], r["v_uc"])
        g_un = _fmt_g(r["g_uninf"], r["q_uninf"], r["v_uninf"])
        g_direct = _fmt_g(r["g_direct"], r["q_direct"], r["v_direct"])
        g_cd = _fmt_g(r["g_cd"], r["q_cd"])
        g_uccd = _fmt_g(r["g_uccd"], r["q_uccd"])
        constit, infldep = r["constit"], r["infldep"]
        lines.append(f"  {gene_tex} & {g_uc} & {g_un} & {g_direct} & {g_cd}"
                      f" & {g_uccd} & {constit} & {infldep} \\\\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(
        r"\caption{Cross-comparison gene-level effect sizes for key genes"
        r" discussed in the text. Each column shows the meta-analytic"
        r" Hedges'~$g$ from one of the five analyses; significance is indicated"
        r" by stars ($^{*}q < 0.05$, $^{**}q < 0.01$, $^{***}q < 0.001$)."
        r" The last two columns mark the two activity layers, which are independent"
        r" properties rather than exclusive categories, so a gene may carry both or"
        r" neither: Constit = altered already in uninflamed mucosa; Infl-dep ="
        r" measurably changed by inflammation."
        r" A dagger marks an effect the independent RNA-seq cohorts confirm --- significant"
        r" there in the same direction; only the three UC contrasts have an RNA-seq"
        r" counterpart, so the CD and UC\,vs\,CD columns carry no dagger."
        r" Stars indicate $q$-values only; the layers additionally require direction"
        r" consistency $\geq 80\%$ and presence in at least half the contributing"
        r" cohorts, so a starred effect may still be unclassified."
        r" The UC\,vs\,CD column shows"
        r" the direct comparison (positive = higher in UC).}"
    )
    lines.append(r"\label{tab:cross_comparison}")
    lines.append(r"\end{table*}")

    tex = "\n".join(lines)
    with open(args.output, "w") as f:
        f.write(tex)
    print(f"Saved {args.output}")

    # Also write a CSV version for convenience
    csv_out = args.output.replace(".tex", ".csv")
    data_rows = [r for r in rows if "_group" not in r]
    pd.DataFrame(data_rows).to_csv(csv_out, index=False)
    print(f"Saved {csv_out}")


if __name__ == "__main__":
    main()
