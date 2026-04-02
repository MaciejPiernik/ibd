"""
Cross-comparison gene table for the Discussion section.

Merges gene-level meta-analysis results across four analyses and produces a
LaTeX table showing key genes grouped by functional category, with effect
sizes, significance, and constitutive/inflammation-dependent classification.

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
    "Therapeutic targets — anti-TNF / integrin / S1P": [
        "TNF", "TNFRSF1B", "ITGA4", "ITGB7", "MADCAM1",
        "S1PR1", "S1PR2", "S1PR3", "S1PR4",
    ],
    "Therapeutic targets — IL-23 / JAK-STAT": [
        "IL23A", "IL12A", "IL12B", "JAK1", "JAK2", "JAK3", "STAT1", "STAT3",
    ],
    "Therapeutic targets — inflammasome / complement / CXCL": [
        "CASP1", "IL1B", "CASP4", "CASP5", "GSDMD", "IL18",
        "CFB", "C2", "C3", "C4A", "CD55",
        "CXCL1", "CXCL2", "CXCL3", "CXCL5",
    ],
    "Metabolic axis — mitochondrial regulators": [
        "PPARGC1A", "PPARGC1B", "ESRRA",
    ],
    "Metabolic axis — energy metabolism genes": [
        "OPA1", "PDP2", "IDH3A", "ACSF2", "SLC22A5",
    ],
    "Ferroptosis / ROS / metallothionein": [
        "DUOX2", "DUOXA2", "ACSL4", "GPX4", "FTH1", "SLC40A1",
        "MT1E", "MT1F", "MT1G", "MT1M", "MT1X", "MT2A",
    ],
    "Barrier integrity": [
        "CLDN8", "AQP8", "OCLN", "SLC26A2",
    ],
    "Tryptophan / kynurenine pathway": [
        "IDO1", "KYNU", "HAAO", "QPRT",
    ],
    "Glucuronidation (selected UGTs)": [
        "UGT1A1", "UGT1A6", "UGT2A3", "UGT2B7", "UGT2B15", "UGT2B17",
    ],
}


def _load_gene_all(directory: str) -> pd.DataFrame:
    """Load gene_meta_all.csv, handling the space-in-filename issue."""
    for candidate in ["gene_meta_all.csv", " gene_meta_all.csv"]:
        path = os.path.join(directory, candidate)
        if os.path.exists(path):
            return pd.read_csv(path)
    raise FileNotFoundError(f"gene_meta_all.csv not found in {directory}")


def _fmt_g(val: float, q: float) -> str:
    """Format effect size with significance stars."""
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
    return f"{sign}{val:.2f}{stars}"


def _classify(gene: str, uc_row, uninf_row) -> str:
    """Classify as constitutive, inflammation-dependent, or amplified."""
    if uc_row is None or uninf_row is None:
        return "---"

    uc_q = uc_row["q_value"]
    un_q = uninf_row["q_value"]
    uc_g = uc_row["hedges_g"]
    un_g = uninf_row["hedges_g"]

    uc_sig = uc_q < 0.05
    un_sig = un_q < 0.05

    if uc_sig and not un_sig:
        return "Infl-dep"
    elif uc_sig and un_sig:
        # Check ratio
        if abs(uc_g) < 0.1:
            return "Both sig"
        ratio = un_g / uc_g
        if 0.70 <= ratio <= 1.30:
            return "Constit"
        else:
            return "Amplified"
    elif not uc_sig and un_sig:
        return "Uninfl-only"
    else:
        return "NS"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uc-dir", default="results/paper/uc_vs_hc")
    ap.add_argument("--uninf-dir", default="results/paper/uninf_vs_hc")
    ap.add_argument("--cd-dir", default="results/paper/cd_vs_hc")
    ap.add_argument("--uccd-dir", default="results/paper/uc_vs_cd")
    ap.add_argument("--output", default="results/paper/cross_comparison_table.tex")
    args = ap.parse_args()

    uc = _load_gene_all(args.uc_dir).set_index("gene")
    uninf = _load_gene_all(args.uninf_dir).set_index("gene")
    cd = _load_gene_all(args.cd_dir).set_index("gene")
    uccd = _load_gene_all(args.uccd_dir).set_index("gene")

    rows = []
    for group_name, genes in GENE_GROUPS.items():
        rows.append({"_group": group_name})
        for gene in genes:
            r = {"gene": gene}
            for label, df in [("uc", uc), ("uninf", uninf),
                              ("cd", cd), ("uccd", uccd)]:
                if gene in df.index:
                    row = df.loc[gene]
                    r[f"g_{label}"] = row["hedges_g"]
                    r[f"q_{label}"] = row["q_value"]
                else:
                    r[f"g_{label}"] = np.nan
                    r[f"q_{label}"] = np.nan

            uc_row = uc.loc[gene] if gene in uc.index else None
            un_row = uninf.loc[gene] if gene in uninf.index else None
            r["class"] = _classify(gene, uc_row, un_row)
            rows.append(r)

    # ── build LaTeX ───────────────────────────────────────────────────────
    lines = []
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    lines.append(r"\scriptsize")
    lines.append(r"\renewcommand{\arraystretch}{1.10}")
    lines.append(r"\begin{tabular}{l rrrr l}")
    lines.append(r"\toprule")
    lines.append(r"Gene & $g_\mathrm{UC\,infl}$ & $g_\mathrm{UC\,uninfl}$"
                 r" & $g_\mathrm{CD}$ & $g_\mathrm{UC\,vs\,CD}$"
                 r" & Classification \\")
    lines.append(r"\midrule")

    for r in rows:
        if "_group" in r:
            lines.append(r"\midrule")
            lines.append(r"\multicolumn{6}{l}{\textit{"
                         + r["_group"] + r"}} \\")
            continue

        gene_tex = r["gene"].replace("_", r"\_")
        g_uc = _fmt_g(r["g_uc"], r["q_uc"])
        g_un = _fmt_g(r["g_uninf"], r["q_uninf"])
        g_cd = _fmt_g(r["g_cd"], r["q_cd"])
        g_uccd = _fmt_g(r["g_uccd"], r["q_uccd"])
        cl = r["class"]
        lines.append(f"  {gene_tex} & {g_uc} & {g_un} & {g_cd}"
                      f" & {g_uccd} & {cl} \\\\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(
        r"\caption{Cross-comparison gene-level effect sizes for key genes"
        r" discussed in the text. Each column shows the meta-analytic"
        r" Hedges'~$g$ from one of the four analyses; significance is"
        r" indicated by stars ($^{*}q < 0.05$, $^{**}q < 0.01$,"
        r" $^{***}q < 0.001$). Classification: Constit = constitutive"
        r" (significant in both inflamed and uninflamed analyses, ratio"
        r" $0.70$--$1.30$); Amplified = inflammation-amplified (significant"
        r" in both, ratio outside $0.70$--$1.30$); Infl-dep ="
        r" inflammation-dependent (significant only in inflamed);"
        r" NS = not significant in either. The UC\,vs\,CD column shows"
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
