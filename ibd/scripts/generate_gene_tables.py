#!/usr/bin/env python3
"""Generate condensed top-gene tables (Gene, Hedges' g, Pathway).

Usage:
    python generate_gene_tables.py \
        --gene-sig results/gene_meta_significant.csv \
        --pathway-sig results/pathway_meta_significant.csv \
        --output-tex tables/gene_tables.tex \
        --top-n 20
"""
import argparse, re
import pandas as pd


def build_gene_to_pathway(pw_knee):
    gene_best = {}
    for _, row in pw_knee.iterrows():
        pathway = row["term"].split(" R-HSA")[0]
        abs_nes = abs(row["combined_nes"])
        for col in ["lead_genes_core", "lead_genes_all"]:
            if col not in row or pd.isna(row[col]):
                continue
            for gene in str(row[col]).split(";"):
                gene = gene.strip()
                if not gene:
                    continue
                if gene not in gene_best or abs_nes > gene_best[gene][0]:
                    gene_best[gene] = (abs_nes, pathway)
    return {g: v[1] for g, v in gene_best.items()}


def is_noncoding(gene):
    return bool(re.match(r"^(MIR\d|LOC\d|LINC\d|SNORD|SCARNA|SNORA)", gene))


def make_table(df, direction, gene_to_pw, top_n):
    sub = df[df["direction"] == direction].copy()
    sub = sub[~sub["gene"].apply(is_noncoding)]
    sub["abs_g"] = sub["hedges_g"].abs()
    sub = sub.nlargest(top_n, "abs_g")

    rows = []
    for _, r in sub.iterrows():
        pw = gene_to_pw.get(r["gene"], "---")
        rows.append(f"  {r['gene']} & {r['hedges_g']:+.2f} & {pw} \\\\")

    dir_label = "upregulated" if direction == "up" else "downregulated"
    label = "tab:top_up" if direction == "up" else "tab:top_down"

    return (
        "\\begin{table}[t]\n"
        "\\centering\n"
        "\\scriptsize\n"
        "\\renewcommand{\\arraystretch}{1.15}\n"
        "\\begin{tabular}{lrl}\n"
        "\\hline\n"
        "Gene & $g$ & Pathway \\\\\n"
        "\\hline\n"
        + "\n".join(rows) + "\n"
        "\\hline\n"
        "\\end{tabular}\n"
        f"\\caption{{Top {top_n} {dir_label} protein-coding genes by Hedges' $g$. "
        "Pathway: highest-scoring enriched Reactome pathway containing the gene "
        "in its leading edge; ``---'': not in any enriched pathway.}}\n"
        f"\\label{{{label}}}\n"
        "\\end{table}\n"
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gene-sig", required=True)
    p.add_argument("--pathway-sig", required=True,
                   help="pathway_meta_significant.csv (all significant pathways for broader mapping)")
    p.add_argument("--output-tex", required=True)
    p.add_argument("--top-n", type=int, default=20)
    args = p.parse_args()

    sig = pd.read_csv(args.gene_sig, encoding="latin1")
    pw = pd.read_csv(args.pathway_sig, encoding="latin1")
    g2p = build_gene_to_pathway(pw)

    print(f"Significant genes: {len(sig)} "
          f"({(sig['direction']=='up').sum()} up, "
          f"{(sig['direction']=='down').sum()} down)")
    print(f"Mapped to pathways: {sum(1 for g in sig['gene'] if g in g2p)}")

    with open(args.output_tex, "w") as f:
        f.write(make_table(sig, "up", g2p, args.top_n))
        f.write("\n")
        f.write(make_table(sig, "down", g2p, args.top_n))

    print(f"Written: {args.output_tex}")


if __name__ == "__main__":
    main()
