#!/usr/bin/env python3
"""
Generate a publication-ready SVG methodology flowchart for the UC meta-analysis paper.
Outputs: methods_figure.svg
"""

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import SubElement

# ── Layout constants ──
W = 1000
CX = W / 2
BOX_H = 50
BOX_H2 = 68
BOX_RX = 6
GAP_V = 36
ARROW_GAP = 6
BW = 280
FONT = "Arial, Helvetica, sans-serif"
FS = 16
FS_SM = 12
FS_XS = 11
FS_TRACK = 18
FS_HEAD = 16

# ── Colors ──
C_DATA  = "#E3E3E3"
C_PROC  = "#CFE4CF"
C_COMP  = "#D5CCE8"
C_DE    = "#FFD9B3"
C_ANAL  = "#B3D4F0"
C_FILT  = "#FFF3B3"
C_OUT   = "#C8E6C8"
C_VALID = "#F5F5F5"
C_VBDR  = "#AAAAAA"
C_PANEL = "#EDF3FA"
C_PSTK  = "#C0D4E8"

S_MAIN  = "#555555"
S_LIGHT = "#888888"
S_ARR   = "#444444"
T1      = "#1a1a1a"
T2      = "#444444"
T_V     = "#555555"


def s(v):
    return str(round(v, 1))


def make_svg():
    svg = ET.Element("svg", xmlns="http://www.w3.org/2000/svg")
    svg.set("version", "1.1")

    defs = SubElement(svg, "defs")
    for mid, color, sz in [("arr", S_ARR, 7), ("arr_v", C_VBDR, 6)]:
        m = SubElement(defs, "marker", id=mid, viewBox="0 0 10 10",
                       refX="9", refY="5", markerWidth=str(sz),
                       markerHeight=str(sz), orient="auto-start-reverse")
        SubElement(m, "path", d="M1 1.5L8 5L1 8.5", fill="none",
                   stroke=color, **{"stroke-width": "1.5",
                   "stroke-linecap": "round", "stroke-linejoin": "round"})

    g = SubElement(svg, "g")

    def rect(x, y, w, h, fill, stroke=S_MAIN, sw="0.8", rx=BOX_RX, dash=False):
        a = {"stroke-width": sw}
        if dash:
            a["stroke-dasharray"] = "4 2"
        SubElement(g, "rect", x=s(x), y=s(y), width=s(w), height=s(h),
                   rx=s(rx), fill=fill, stroke=stroke, **a)

    def txt(x, y, text, anchor="middle", size=FS, weight="400",
            color=T1, style="normal", spacing="0"):
        t = SubElement(g, "text", x=s(x), y=s(y), fill=color, **{
            "text-anchor": anchor, "dominant-baseline": "central",
            "font-family": FONT, "font-size": s(size),
            "font-weight": weight, "font-style": style,
            "letter-spacing": spacing})
        t.text = text

    def box(x, y, w, h, fill, label, sub=None, bold=True):
        rect(x, y, w, h, fill)
        if sub:
            txt(x + w/2, y + h/2 - 9, label, weight="600" if bold else "400")
            txt(x + w/2, y + h/2 + 10, sub, size=FS_SM, color=T2)
        else:
            txt(x + w/2, y + h/2, label, weight="600" if bold else "400")

    def arrv(x, y1, y2):
        SubElement(g, "line", x1=s(x), y1=s(y1), x2=s(x), y2=s(y2),
                   stroke=S_ARR, **{"stroke-width": "1.2",
                   "marker-end": "url(#arr)"})

    def linev(x, y1, y2):
        SubElement(g, "line", x1=s(x), y1=s(y1), x2=s(x), y2=s(y2),
                   stroke=S_ARR, **{"stroke-width": "0.8"})

    def lineh(x1, x2, y):
        SubElement(g, "line", x1=s(x1), y1=s(y), x2=s(x2), y2=s(y),
                   stroke=S_ARR, **{"stroke-width": "0.8"})

    def vbox(x, y, w, h, lines):
        rect(x, y, w, h, C_VALID, stroke=C_VBDR, sw="0.7", rx=4, dash=True)
        sp = 14
        sy = y + h/2 - (len(lines) - 1) * sp / 2
        for i, ln in enumerate(lines):
            txt(x + w/2, sy + i * sp, ln, size=FS_XS, color=T_V, style="italic")

    def vconn(x1, y1, x2, y2):
        SubElement(g, "line", x1=s(x1), y1=s(y1), x2=s(x2), y2=s(y2),
                   stroke=C_VBDR, **{"stroke-width": "0.8",
                   "stroke-dasharray": "3 3", "marker-end": "url(#arr_v)"})

    y = 20

    # ── Data acquisition ──
    aw = 420
    box(CX - aw/2, y, aw, BOX_H2, C_DATA,
        "GEO database search",
        "33 candidate datasets → 14 retained (1,470 samples, 9 platforms)")
    y += BOX_H2
    arrv(CX, y + ARROW_GAP, y + GAP_V - ARROW_GAP)
    y += GAP_V

    # ── Data processing ──
    pw = 270
    pg = 40
    pxl = CX - pg/2 - pw
    pxr = CX + pg/2
    box(pxl, y, pw, BOX_H2, C_PROC,
        "Probe-to-gene mapping",
        "Multi-probe collapse by mean expression")
    box(pxr, y, pw, BOX_H2, C_PROC,
        "Metadata harmonization",
        "14 dataset-specific parsers")
    y += BOX_H2
    ml = pxl + pw/2
    mr = pxr + pw/2
    my = y + 12
    linev(ml, y + 3, my)
    linev(mr, y + 3, my)
    lineh(ml, mr, my)
    arrv(CX, my, my + 18)
    y = my + 20

    # ── Four comparisons ──
    txt(CX, y + 8, "Sample selection: 4 parallel comparisons",
        size=FS_HEAD, weight="700")
    y += 26

    cw = 185
    cg = 10
    cx0 = CX - (4 * cw + 3 * cg) / 2
    ch = 60
    comps = [
        ("UC inflamed vs ctrl",   "12 datasets; 460 UC, 236 ctrl"),
        ("UC uninflamed vs ctrl", "5 datasets; 132 UC, 113 ctrl"),
        ("CD inflamed vs ctrl",   "8 datasets; 123 CD, 134 ctrl"),
        ("UC vs CD direct",       "6 datasets; 197 UC, 108 CD"),
    ]
    centers = []
    for i, (lb, sub) in enumerate(comps):
        bx = cx0 + i * (cw + cg)
        box(bx, y, cw, ch, C_COMP, lb, sub)
        centers.append(bx + cw / 2)
    y += ch

    by = y + 8
    for xc in centers:
        linev(xc, y + 2, by)
    lineh(centers[0], centers[-1], by)
    ay = by + 14
    txt(CX, ay, "Same analytical pipeline applied to each",
        size=FS_SM, style="italic", color=T2)
    arrv(CX, ay + 8, ay + GAP_V - 4)
    y = ay + GAP_V - 2

    # ── Per-dataset differential expression ──
    txt(CX, y + 8, "Per-dataset differential expression",
        size=FS_HEAD, weight="700")
    y += 26

    de_w = 280
    de_gap = 60
    de_xl = CX - de_gap/2 - de_w
    de_xr = CX + de_gap/2

    GT = de_xl + de_w / 2
    PT = de_xr + de_w / 2
    gx = GT - BW / 2
    px = PT - BW / 2

    box(de_xl, y, de_w, BOX_H2, C_DE,
        "Hedges' g effect sizes",
        "Standardized mean difference + sampling variance")
    box(de_xr, y, de_w, BOX_H2, C_DE,
        "limma moderated t-statistics",
        "Empirical Bayes; BH-corrected q-values")
    de_bot = y + BOX_H2

    # ── Gene track + Pathway track ──
    PP_X = 14
    PP_TOP = 34
    PP_BOT = 14

    r1 = de_bot + GAP_V + PP_TOP
    r2 = r1 + BOX_H2 + GAP_V
    r3 = r2 + BOX_H2 + GAP_V
    fb = r3 + BOX_H2 + GAP_V
    fb_h = 42
    r4 = fb + fb_h + GAP_V
    r5 = r4 + BOX_H2 + GAP_V
    r6 = r5 + BOX_H2 + GAP_V

    ptop = r1 - PP_TOP
    pbot = r6 + BOX_H + PP_BOT

    for tx in [gx, px]:
        rect(tx - PP_X, ptop, BW + 2 * PP_X, pbot - ptop,
             C_PANEL, stroke=C_PSTK, sw="0.7", rx=10)

    txt(GT, ptop + 18, "GENE TRACK",
        size=FS_TRACK, weight="700", color="#5A85AD", spacing="1.5")
    txt(PT, ptop + 18, "PATHWAY TRACK",
        size=FS_TRACK, weight="700", color="#5A85AD", spacing="1.5")

    arrv(GT, de_bot + ARROW_GAP, ptop - ARROW_GAP)
    arrv(PT, de_bot + ARROW_GAP, ptop - ARROW_GAP)

    # Row 1: Gene-level REML / GSEA
    box(gx, r1, BW, BOX_H2, C_ANAL,
        "Gene-level REML meta-analysis",
        "τ² by REML; inverse-variance weighting; Wald test")
    box(px, r1, BW, BOX_H2, C_ANAL,
        "GSEA (Reactome 2022, preranked)",
        "1,000 permutations; NES + leading-edge genes")

    # RRA validation (left side → Gene-level REML)
    rw, rh = 160, 56
    rx = gx - PP_X - rw - 6
    ry = r1 + 6
    vbox(rx, ry, rw, rh,
         ["Robust Rank Aggregation",
          "Ranking by |t| from limma",
          "(ρ = 0.89 with |g| rankings)"])
    vconn(rx + rw, ry + rh/2, gx - 2, ry + rh/2)

    # Row 2: Heterogeneity / Pathway-level REML
    arrv(GT, r1 + BOX_H2 + ARROW_GAP, r2 - ARROW_GAP)
    arrv(PT, r1 + BOX_H2 + ARROW_GAP, r2 - ARROW_GAP)

    box(gx, r2, BW, BOX_H2, C_ANAL,
        "Heterogeneity assessment",
        "I² quantification; |g|–I² correlation (ρ = 0.74)")
    box(px, r2, BW, BOX_H2, C_ANAL,
        "Pathway-level REML meta-analysis",
        "NES variance ≈ (NES/z)²")

    # Stouffer validation (right side → Pathway-level REML)
    sw_, sh = 160, 52
    sx = px + BW + PP_X + 6
    sy = r2 + 8
    vbox(sx, sy, sw_, sh,
         ["Stouffer's weighted z-method",
          "(91.3% concordance, r = 0.89)"])
    vconn(sx, sy + sh/2, px + BW + 2, sy + sh/2)

    # Row 3: pathway — LE aggregation; gene track continues
    arrv(PT, r2 + BOX_H2 + ARROW_GAP, r3 - ARROW_GAP)
    linev(GT, r2 + BOX_H2 + ARROW_GAP, fb - ARROW_GAP)

    box(px, r3, BW, BOX_H2, C_ANAL,
        "Leading-edge gene aggregation",
        "Core LE genes: present in ≥ 50% of datasets")
    linev(PT, r3 + BOX_H2 + ARROW_GAP, fb - ARROW_GAP)

    # Shared significance filter bar
    fb_x = gx - PP_X
    fb_w = px + BW + PP_X - fb_x
    rect(fb_x, fb, fb_w, fb_h, "#E8E4F0", stroke=S_LIGHT, sw="0.6", rx=4)
    txt(fb_x + fb_w/2, fb + fb_h/2,
        "Significance filters: q < 0.05, direction consistency ≥ 80%, ≥ 50% datasets",
        size=FS_SM, weight="500")

    arrv(PT, fb + fb_h + ARROW_GAP, r4 - ARROW_GAP)
    arrv(GT, fb + fb_h + ARROW_GAP, r6 - ARROW_GAP)

    # Row 4: pathway — knee-point
    box(px, r4, BW, BOX_H2, C_FILT,
        "Knee-point filtering",
        "Data-driven cutoff on sorted |NES|")

    # Row 5: pathway — Louvain
    arrv(PT, r4 + BOX_H2 + ARROW_GAP, r5 - ARROW_GAP)
    box(px, r5, BW, BOX_H2, C_FILT,
        "Louvain community detection",
        "Overlap coefficient network (threshold ≥ 0.33)")

    # Row 6: both result boxes
    arrv(PT, r5 + BOX_H2 + ARROW_GAP, r6 - ARROW_GAP)
    box(gx, r6, BW, BOX_H, C_OUT,
        "Significant genes per comparison", bold=True)
    box(px, r6, BW, BOX_H, C_OUT,
        "Pathway clusters per comparison", bold=True)

    # ── Cross-comparison synthesis ──
    merge_y = r6 + BOX_H + 28
    linev(GT, r6 + BOX_H + ARROW_GAP, merge_y)
    linev(PT, r6 + BOX_H + ARROW_GAP, merge_y)
    lineh(GT, PT, merge_y)
    arrv((GT + PT) / 2, merge_y, merge_y + 22)

    out_y = merge_y + 24
    out_x = gx - 10
    out_w = px + BW - gx + 20
    box(out_x, out_y, out_w, BOX_H2, C_OUT,
        "Cross-comparison synthesis",
        "Constitutive vs inflammatory (Wald test on Δg, ratio filter 0.70–1.30); UC vs CD pathway map")

    # UC–CD concordance validation (right side, vertically centered with synthesis box)
    vh = 56
    vy = out_y + (BOX_H2 - vh) / 2
    vbox(sx, vy, sw_, vh,
         ["Indirect vs direct UC–CD",
          "(Pearson r = 0.74)"])
    vconn(sx, vy + vh/2, out_x + out_w + 2, vy + vh/2)

    final_y = max(out_y + BOX_H2, vy + vh) + 24
    svg.set("width", str(W))
    svg.set("height", s(final_y))
    svg.set("viewBox", f"0 0 {W} {final_y}")
    return svg


def main():
    svg = make_svg()
    tree = ET.ElementTree(svg)
    ET.indent(tree, space="  ")
    out = "figures/methods_figure.svg"
    tree.write(out, encoding="unicode", xml_declaration=True)
    print(f"Written to {out}")
    print(f"Dimensions: {svg.get('width')} x {svg.get('height')}")


if __name__ == "__main__":
    main()
