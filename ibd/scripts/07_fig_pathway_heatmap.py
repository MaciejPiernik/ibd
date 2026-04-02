#!/usr/bin/env python3
"""Generate a three-way pathway comparison heatmap (butterfly layout).

Left panel: upregulated pathways (labels on left, cells on right).
Right panel: downregulated pathways (cells on left, labels on right).
Heatmap cells meet in the middle. Colorbar at bottom.

Cells are visually coded in three tiers:
  ■  Full colour          – significant (q_value < 0.05)
  ■  Desaturated + hatch  – not significant (q_value ≥ 0.05)
  ■  Light grey            – pathway not tested / missing

Usage:
    python ibd/scripts/07_fig_pathway_heatmap.py \
        --uc-clusters results/paper/uc_vs_hc/hierarchical_pathway_clusters.csv \
        --uninfl-clusters results/paper/uninf_vs_hc/hierarchical_pathway_clusters.csv \
        --cd-clusters results/paper/cd_vs_hc/hierarchical_pathway_clusters.csv \
        --uc-pw-all results/paper/uc_vs_hc/pathway_meta_all.csv \
        --uninfl-pw-all results/paper/uninf_vs_hc/pathway_meta_all.csv \
        --cd-pw-all results/paper/cd_vs_hc/pathway_meta_all.csv \
        --output-pdf figures/pathway_heatmap.pdf \
        --sig-threshold 0.05
"""
import argparse, os, colorsys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def clean_name(term):
    return term.split(" R-HSA")[0]


def desaturate(rgba, factor=0.30):
    """Blend an RGBA colour toward grey, keeping `factor` of the original
    saturation.  Returns an (r, g, b, a) tuple."""
    r, g, b, a = rgba
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    # push lightness toward 0.85 (pale) and crush saturation
    l_new = l + (0.85 - l) * (1 - factor)
    s_new = s * factor
    r2, g2, b2 = colorsys.hls_to_rgb(h, l_new, s_new)
    return (r2, g2, b2, a)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--uc-clusters", required=True)
    p.add_argument("--uninfl-clusters", required=True)
    p.add_argument("--cd-clusters", required=True)
    p.add_argument("--uc-pw-all", required=True)
    p.add_argument("--uninfl-pw-all", required=True)
    p.add_argument("--cd-pw-all", required=True)
    p.add_argument("--output-pdf", required=True)
    p.add_argument("--sig-threshold", type=float, default=0.05,
                   help="q-value threshold for significance (default 0.05)")
    args = p.parse_args()

    SIG_THR = args.sig_threshold
    os.makedirs(os.path.dirname(args.output_pdf) or ".", exist_ok=True)

    # ── Load cluster representatives ──
    uc_cl = pd.read_csv(args.uc_clusters, encoding="latin1")
    un_cl = pd.read_csv(args.uninfl_clusters, encoding="latin1")
    cd_cl = pd.read_csv(args.cd_clusters, encoding="latin1")

    all_reps = set()
    for cl in [uc_cl, un_cl, cd_cl]:
        all_reps.update(cl["Representative"].tolist())

    # ── NES + significance lookups ──
    uc_pw = pd.read_csv(args.uc_pw_all, encoding="latin1")
    un_pw = pd.read_csv(args.uninfl_pw_all, encoding="latin1")
    cd_pw = pd.read_csv(args.cd_pw_all, encoding="latin1")

    uc_nes = dict(zip(uc_pw["term"], uc_pw["combined_nes"]))
    un_nes = dict(zip(un_pw["term"], un_pw["combined_nes"]))
    cd_nes = dict(zip(cd_pw["term"], cd_pw["combined_nes"]))

    uc_qval = dict(zip(uc_pw["term"], uc_pw["q_value"]))
    un_qval = dict(zip(un_pw["term"], un_pw["q_value"]))
    cd_qval = dict(zip(cd_pw["term"], cd_pw["q_value"]))

    uc_ds = dict(zip(uc_pw["term"], uc_pw["datasets"]))
    un_ds = dict(zip(un_pw["term"], un_pw["datasets"]))
    cd_ds = dict(zip(cd_pw["term"], cd_pw["datasets"]))

    uc_dr = dict(zip(uc_pw["term"], uc_pw["direction_ratio"]))
    un_dr = dict(zip(un_pw["term"], un_pw["direction_ratio"]))
    cd_dr = dict(zip(cd_pw["term"], cd_pw["direction_ratio"]))

    uc_ds_max = uc_pw["datasets"].max()
    un_ds_max = un_pw["datasets"].max()
    cd_ds_max = cd_pw["datasets"].max()

    # ── Build data ──
    data = []
    for term in all_reps:
        name = clean_name(term)
        if len(name) > 55:
            name = name[:52] + "..."
        uc_val = uc_nes.get(term, np.nan)
        un_val = un_nes.get(term, np.nan)
        cd_val = cd_nes.get(term, np.nan)
        sort_nes = uc_val if not np.isnan(uc_val) else (
            cd_val if not np.isnan(cd_val) else un_val)

        # Significance flags (True = significant)
        # Requires: q_value < threshold AND datasets >= max/2 AND direction_ratio >= 0.8
        def _is_sig(qval_dict, ds_dict, dr_dict, ds_max, term):
            q = qval_dict.get(term, np.nan)
            d = ds_dict.get(term, np.nan)
            r = dr_dict.get(term, np.nan)
            if np.isnan(q):
                return np.nan
            return (q < SIG_THR) and (d >= ds_max / 2) and (r >= 0.8)

        uc_sig = _is_sig(uc_qval, uc_ds, uc_dr, uc_ds_max, term)
        un_sig = _is_sig(un_qval, un_ds, un_dr, un_ds_max, term)
        cd_sig = _is_sig(cd_qval, cd_ds, cd_dr, cd_ds_max, term)

        data.append({"name": name,
                      "uc": uc_val, "un": un_val, "cd": cd_val,
                      "uc_sig": uc_sig, "un_sig": un_sig, "cd_sig": cd_sig,
                      "sort_nes": sort_nes})

    df = pd.DataFrame(data)
    df_up = df[df["sort_nes"] > 0].sort_values(
        "sort_nes", ascending=False).reset_index(drop=True)
    df_down = df[df["sort_nes"] <= 0].sort_values(
        "sort_nes", ascending=True).reset_index(drop=True)

    mat_up   = df_up[["uc", "un", "cd"]].values
    sig_up   = df_up[["uc_sig", "un_sig", "cd_sig"]].values
    mat_down = df_down[["uc", "un", "cd"]].values
    sig_down = df_down[["uc_sig", "un_sig", "cd_sig"]].values

    n_up, n_down = len(df_up), len(df_down)
    n_max = max(n_up, n_down)

    # Pad shorter panel so both have same height
    def pad_matrices(mat, sig, names_list, n_actual, n_target):
        if n_actual < n_target:
            pad_n = n_target - n_actual
            mat = np.vstack([mat, np.full((pad_n, 3), np.nan)])
            sig = np.vstack([sig, np.full((pad_n, 3), np.nan)])
            names_list = names_list + [""] * pad_n
        return mat, sig, names_list

    mat_up, sig_up, names_up = pad_matrices(
        mat_up, sig_up, df_up["name"].tolist(), n_up, n_max)
    mat_down, sig_down, names_down = pad_matrices(
        mat_down, sig_down, df_down["name"].tolist(), n_down, n_max)

    # ── Global colour scale ──
    all_vals = np.concatenate([mat_up.ravel(), mat_down.ravel()])
    vmax = np.nanmax(np.abs(all_vals[np.isfinite(all_vals)]))
    cmap = plt.cm.RdBu_r
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    col_labels = ["UC\ninfl.", "UC\nuninfl.", "CD\ninfl."]

    # ── Layout ──
    cell_w = 1.0
    cell_h = 1.0
    gap = 0.4

    x_up   = [0, cell_w, 2 * cell_w]
    x_down = [3 * cell_w + gap, 4 * cell_w + gap, 5 * cell_w + gap]

    fig_w = 13.5
    fig_h = max(5.5, 0.36 * n_max + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ── Cell drawing ──
    def draw_cell(x, y, val, significant):
        """Draw one heatmap cell.

        significant can be True, False, or nan (= pathway missing).
        """
        # --- missing / not tested ---
        if np.isnan(val):
            rect = plt.Rectangle((x, y), cell_w, cell_h,
                                  facecolor="#ededed", edgecolor="white",
                                  linewidth=0.8)
            ax.add_patch(rect)
            return

        base_color = cmap(norm(val))

        if significant is True or significant == 1.0:
            # Significant: full colour, no overlay
            face = base_color
            hatch = None
        else:
            # Not significant: desaturated + diagonal hatch
            face = desaturate(base_color, factor=0.1)
            hatch = "////"

        rect = plt.Rectangle((x, y), cell_w, cell_h,
                              facecolor=face, edgecolor="white",
                              linewidth=0.8, hatch=hatch)
        if hatch:
            rect.set_edgecolor("white")       # keep border white
            # matplotlib draws hatches in the edge colour; override via a
            # second invisible-fill patch that carries the hatch colour
            ax.add_patch(rect)
            hatch_rect = plt.Rectangle(
                (x, y), cell_w, cell_h,
                facecolor="none", edgecolor="#888888",
                linewidth=0, hatch=hatch, alpha=0.45)
            ax.add_patch(hatch_rect)
        else:
            ax.add_patch(rect)

        # NES text
        text = f"{val:+.1f}"
        brightness = (0.299 * face[0] + 0.587 * face[1] + 0.114 * face[2])
        txt_color = "white" if brightness < 0.55 else "#222222"
        ax.text(x + cell_w / 2, y + cell_h / 2, text,
                ha="center", va="center", fontsize=6.5,
                color=txt_color, fontweight="medium")

    # Draw upregulated panel (left)
    for i in range(n_max):
        y = (n_max - 1 - i) * cell_h
        for j in range(3):
            if i < n_up:
                draw_cell(x_up[j], y, mat_up[i, j], sig_up[i, j])
            else:
                rect = plt.Rectangle((x_up[j], y), cell_w, cell_h,
                                     facecolor="#fafafa", edgecolor="#fafafa",
                                     linewidth=0)
                ax.add_patch(rect)

    # Draw downregulated panel (right)
    for i in range(n_max):
        y = (n_max - 1 - i) * cell_h
        for j in range(3):
            if i < n_down:
                draw_cell(x_down[j], y, mat_down[i, j], sig_down[i, j])
            else:
                rect = plt.Rectangle((x_down[j], y), cell_w, cell_h,
                                     facecolor="#fafafa", edgecolor="#fafafa",
                                     linewidth=0)
                ax.add_patch(rect)

    # ── Row labels ──
    for i, name in enumerate(names_up):
        if not name:
            continue
        y = (n_max - 1 - i) * cell_h + cell_h / 2
        ax.text(x_up[0] - 0.15, y, name,
                ha="right", va="center", fontsize=6.5)

    for i, name in enumerate(names_down):
        if not name:
            continue
        y = (n_max - 1 - i) * cell_h + cell_h / 2
        ax.text(x_down[2] + cell_w + 0.15, y, name,
                ha="left", va="center", fontsize=6.5)

    # ── Column headers ──
    header_y = n_max * cell_h + 0.15
    for j, label in enumerate(col_labels):
        ax.text(x_up[j] + cell_w / 2, header_y, label,
                ha="center", va="bottom", fontsize=7, fontweight="bold")
        ax.text(x_down[j] + cell_w / 2, header_y, label,
                ha="center", va="bottom", fontsize=7, fontweight="bold")

    # ── Panel titles ──
    title_y = n_max * cell_h + 1.5
    ax.text((x_up[0] + x_up[2] + cell_w) / 2, title_y,
            "Upregulated", ha="center", va="bottom",
            fontsize=9, fontweight="bold", color="#b22222")
    ax.text((x_down[0] + x_down[2] + cell_w) / 2, title_y,
            "Downregulated", ha="center", va="bottom",
            fontsize=9, fontweight="bold", color="#1a5276")

    # ── Continuation dots ──
    if n_up < n_max:
        ell_y_up = (n_max - n_up) * cell_h - 0.3
    else:
        ell_y_up = -0.5
    ax.text((x_up[0] + x_up[2] + cell_w) / 2, ell_y_up,
            "···", ha="center", va="center", fontsize=12, color="#888888")

    if n_down < n_max:
        ell_y_down = (n_max - n_down) * cell_h - 0.3
    else:
        ell_y_down = -0.5
    ax.text((x_down[0] + x_down[2] + cell_w) / 2, ell_y_down,
            "···", ha="center", va="center", fontsize=12, color="#888888")

    # ── Colourbar ──
    cbar_ax = fig.add_axes([0.4123, 0.1, 0.2, 0.018])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Normalized Enrichment Score (NES) of significant pathways", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)

    # ── Axis limits ──
    margin_left = 0.05
    margin_right = 0.05
    ax.set_xlim(x_up[0] - margin_left, x_down[2] + cell_w + margin_right)
    ax.set_ylim(-1.0, n_max * cell_h + 2.2)
    ax.set_aspect("equal")
    ax.axis("off")

    plt.savefig(args.output_pdf, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Written: {args.output_pdf} "
          f"({n_up} up + {n_down} down = {n_up + n_down} pathways)")


if __name__ == "__main__":
    main()