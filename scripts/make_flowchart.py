"""Figure 1: study design as an SVG schematic, exported to PNG (600 dpi) and PDF.

Each row: data (left), a pictogram of the comparison (middle) and a small plot of the readout (right).
Marks in the readout plots use the estimates in results/.
Canvas units: 10 per mm, 180 mm wide, so 8 pt is about 28 units.
Usage: python3 scripts/make_flowchart.py   (needs rsvg-convert)

Author: Thai Huu Loc Tran
"""
import json, os, subprocess
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R, F = os.path.join(ROOT, "results/"), os.path.join(ROOT, "figures/")
os.makedirs(F, exist_ok=True)


def get(name, s, floor=0.25, kind="scaled"):
    return next(r for r in json.load(open(R + name + ".json"))["runs"] if r["set"] == s and r["umi_floor"] == floor and r["matrix"] == kind)


FIVE = ["k562_two_screens", "within_k562ess", "within_rpe1", "within_hepg2", "within_jurkat", "within_hct116", "within_hek293t"]
SHORT = ["K562 x2", "K562", "RPE1", "HepG2", "Jurkat", "HCT", "HEK"]
k562 = get("k562_two_screens", "active_pairs")
share = [get(f, "active_pairs")["f_asym"] for f in FIVE]
cpx = [get(f, "complex_active")["f_asym"] for f in FIVE]
tf = []
for f in FIVE:
    t = [get(f, "tf_direction_select_screen1_read_screen2"), get(f, "tf_direction_select_screen2_read_screen1")]
    tf.append(sum(r["concordant_active"] for r in t) / sum(r["n_active"] for r in t))
gl = [float(pd.read_csv(R + f + ".gene_level.summary.tsv", sep="\t", nrows=1).share_of_V_asym_gene_level.iloc[0]) for f in FIVE]
lr = pd.read_csv(R + "linear_response.tsv", sep="\t")
fit, ceil = lr[lr.predictor == "fitted_scale"].r_prediction.tolist(), lr[lr.predictor == "fitted_scale"].r_ceiling.tolist()

W = 1800
INK, GREY, LINE = "#111111", "#8C96A8", "#DCE1EA"
BLUE, RED, GREEN, PURPLE, ORANGE = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A"
FONT = "Arial, Helvetica, sans-serif"
svg = []
add = svg.append


def text(x, y, s, size=27, weight="normal", anchor="start", fill=INK):
    add(f"<text x='{x:.1f}' y='{y:.1f}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' text-anchor='{anchor}' fill='{fill}'>{s}</text>")


def line(x1, y1, x2, y2, col=INK, w=3, dash=None, head=False):
    d = f" stroke-dasharray='{dash}'" if dash else ""
    h = f" marker-end='url(#ah{col[1:]})'" if head else ""
    add(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='{col}' stroke-width='{w}'{d}{h}/>")


def gene(x, y, label, col, r=30):
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='{col}' fill-opacity='0.2' stroke='{col}' stroke-width='3'/>")
    text(x, y + 9, label, 26, "bold", "middle", col)


def cell(x, y, col, r=11):
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='{col}' fill-opacity='0.25' stroke='{col}' stroke-width='3'/><circle cx='{x:.1f}' cy='{y:.1f}' r='{r * 0.35:.1f}' fill='{col}'/>")


def matrix(x, y, n, s, fn):
    """n x n grid of squares with side s; fn(i, j) gives (colour, opacity) or None."""
    for i in range(n):
        for j in range(n):
            c = fn(i, j)
            if c:
                add(f"<rect x='{x + j * s:.1f}' y='{y + i * s:.1f}' width='{s - 3}' height='{s - 3}' fill='{c[0]}' fill-opacity='{c[1]}'/>")
            else:
                add(f"<rect x='{x + j * s:.1f}' y='{y + i * s:.1f}' width='{s - 3}' height='{s - 3}' fill='white' stroke='{LINE}' stroke-width='2'/>")


def bars(x0, y0, w, h, vals, col, labels, ymax, fmt="{:.2f}", ref=None, ylab=None, opac=None, lift=None):
    line(x0, y0, x0 + w, y0, INK, 2.5); line(x0, y0, x0, y0 - h, INK, 2.5)
    step = w / len(vals)
    if ref is not None:
        line(x0, y0 - ref / ymax * h, x0 + w, y0 - ref / ymax * h, GREY, 2, "6 5")
    for i, v in enumerate(vals):
        bh = v / ymax * h
        add(f"<rect x='{x0 + i * step + step * 0.18:.1f}' y='{y0 - bh:.1f}' width='{step * 0.64:.1f}' height='{bh:.1f}' fill='{col}' fill-opacity='{opac[i] if opac else 1}'/>")
        text(x0 + (i + 0.5) * step, y0 - bh - 8 - (lift[i] if lift else 0), fmt.format(v), 21, anchor="middle")
        text(x0 + (i + 0.5) * step, y0 + 28, labels[i], 20 if len(vals) > 5 else 21, anchor="middle")
    if ylab:
        add(f"<text transform='translate({x0 - 14},{y0 - h / 2}) rotate(-90)' font-family='{FONT}' font-size='23' text-anchor='middle' fill='{INK}'>{ylab}</text>")


add("<svg xmlns='http://www.w3.org/2000/svg' width='180mm' height='@HMM@mm' viewBox='0 0 1800 @H@'>")
add("<defs>" + "".join(f"<marker id='ah{c[1:]}' viewBox='0 0 10 10' refX='8' refY='5' markerWidth='5' markerHeight='5' orient='auto'><path d='M0,0 L10,5 L0,10 Z' fill='{c}'/></marker>" for c in (GREY, BLUE, RED, GREEN, PURPLE, ORANGE, INK)) + "</defs>")
add("<rect width='1800' height='@H@' fill='white'/>")
C1, C2, C3 = 40, 520, 1190
text(C1, 58, "Data", 34, "bold"); text(C2, 58, "Comparison", 34, "bold"); text(C3, 58, "Readout", 34, "bold")
RH, Y0 = 270, 90


def row(k, col, title, lines_, fig):
    y = Y0 + k * RH
    add(f"<rect x='{C1 - 20}' y='{y + 12}' width='8' height='{RH - 34}' fill='{col}'/>")
    text(C1, y + 44, title, 29, "bold", fill=col)
    for i, s in enumerate(lines_):
        text(C1, y + 84 + i * 34, s, 25)
    text(1760, y + 44, fig, 27, "bold", "end", col)
    if k:
        line(C1, y, W - 40, y, LINE, 2)
    return y


# ------------------------------------------------ row 1: two K562 screens, symmetric and skew-symmetric parts
y = row(0, BLUE, "K562, two screens", ["Replogle, genome-wide", "Replogle, essential genes", f"{k562['n_genes']:,} genes knocked down", "and measured in both"], "Fig. 2")
mx, my, s = C2 + 10, y + 70, 34
sym = {(0, 1): 0.9, (1, 0): 0.9, (2, 3): 0.5, (3, 2): 0.5, (0, 3): 0.3, (3, 0): 0.3}
skw = {(0, 2): (BLUE, 0.9), (2, 0): (RED, 0.9), (1, 3): (BLUE, 0.5), (3, 1): (RED, 0.5)}
full = {**{k: (GREEN, v) for k, v in sym.items()}, (0, 2): (BLUE, 0.9), (1, 3): (BLUE, 0.5)}
matrix(mx, my, 4, s, lambda i, j: None if i == j else full.get((i, j)))
text(mx + 2 * s, my + 4 * s + 30, "response", 23, anchor="middle"); text(mx + 2 * s, my + 4 * s + 56, "matrix", 23, anchor="middle")
text(mx + 4 * s + 22, my + 2 * s + 8, "=", 34, anchor="middle")
matrix(mx + 190, my, 4, s, lambda i, j: None if i == j else ((GREEN, sym[(i, j)]) if (i, j) in sym else None))
text(mx + 190 + 2 * s, my + 4 * s + 30, "symmetric", 23, anchor="middle", fill=GREEN); text(mx + 190 + 2 * s, my + 4 * s + 56, "part", 23, anchor="middle", fill=GREEN)
text(mx + 190 + 4 * s + 22, my + 2 * s + 8, "+", 34, anchor="middle")
matrix(mx + 380, my, 4, s, lambda i, j: None if i == j else skw.get((i, j)))
text(mx + 380 + 2 * s, my + 4 * s + 30, "skew-symmetric", 23, anchor="middle", fill=BLUE); text(mx + 380 + 2 * s, my + 4 * s + 56, "part", 23, anchor="middle", fill=BLUE)
bars(C3 + 60, y + 215, 200, 130, [k562["r_S"], k562["r_A"]], GREEN, ["symmetric", "skew"], 1.0, ylab="correlation", opac=[1, 1])
add(f"<rect x='{C3 + 60 + 100 + 18:.1f}' y='{y + 215 - k562['r_A'] * 130:.1f}' width='64' height='{k562['r_A'] * 130:.1f}' fill='{BLUE}'/>")
text(C3 + 300, y + 120, "between the", 23); text(C3 + 300, y + 148, "two screens", 23)
text(C3 + 300, y + 196, f"skew share {k562['f_asym']:.2f}", 23, "bold", fill=BLUE)

# ------------------------------------------------ row 2: six cell lines
y = row(1, GREEN, "Six cell lines", ["Replogle: K562, RPE1", "Nadig: HepG2, Jurkat", "Huang: HCT116, HEK293T", "(genome-wide)"], "Fig. 3")
for q, (cx, lab) in enumerate(((C2 + 60, "quarter 1"), (C2 + 190, "quarter 2"), (C2 + 370, "quarter 3"), (C2 + 500, "quarter 4"))):
    for k in range(6):
        cell(cx - 28 + (k % 3) * 28, y + 95 + (k // 3) * 28, GREEN)
    text(cx, y + 165, lab, 22, anchor="middle")
line(C2 + 280, y + 70, C2 + 280, y + 175, GREY, 2.5, "6 5")
text(C2 + 125, y + 205, "screen 1", 23, "bold", "middle"); text(C2 + 435, y + 205, "screen 2", 23, "bold", "middle")
text(C2 + 280, y + 240, "cells of one knockdown, and pairs of cell lines", 22, anchor="middle")
bars(C3 + 60, y + 205, 540, 120, share, BLUE, SHORT, 0.6, ref=0.5, ylab="skew share")

# ------------------------------------------------ row 3: classes of gene pairs
y = row(2, ORANGE, "Classes of gene pairs", ["CORUM: protein complexes", "DoRothEA, TRRUST, CollecTRI:", "transcription factor and target"], "Fig. 4")
gene(C2 + 50, y + 110, "A", ORANGE); gene(C2 + 210, y + 110, "B", ORANGE)
line(C2 + 86, y + 98, C2 + 172, y + 98, ORANGE, 4, head=True); line(C2 + 174, y + 122, C2 + 88, y + 122, ORANGE, 4, head=True)
text(C2 + 130, y + 185, "subunits of", 23, anchor="middle"); text(C2 + 130, y + 213, "one complex", 23, anchor="middle")
gene(C2 + 360, y + 110, "TF", PURPLE, 34); gene(C2 + 530, y + 110, "T", PURPLE)
line(C2 + 400, y + 110, C2 + 492, y + 110, PURPLE, 4, head=True)
text(C2 + 445, y + 185, "transcription factor", 23, anchor="middle"); text(C2 + 445, y + 213, "and target", 23, anchor="middle")
bars(C3 + 60, y + 205, 210, 120, [sum(share) / len(share), sum(cpx) / len(cpx)], ORANGE, ["other", "complex"], 0.6, ref=0.5, ylab="skew share", opac=[0.4, 1])
bars(C3 + 370, y + 215, 210, 100, [100 * tf[3], 100 * tf[4]], PURPLE, ["HepG2", "Jurkat"], 100, fmt="{:.0f}%", ref=50)
text(C3 + 475, y + 84, "larger effect TF to target", 20, anchor="middle")

# ------------------------------------------------ row 4: gene level
y = row(3, PURPLE, "Gene-level terms", ["Effect a gene sends as a", "knockdown and receives", "as a readout"], "Fig. 5")
gene(C2 + 280, y + 130, "G", PURPLE, 36)
for dy in (-50, 0, 50):
    line(C2 + 322, y + 130 + dy * 0.4, C2 + 470, y + 130 + dy, BLUE, 4, head=True)
    line(C2 + 90, y + 130 + dy, C2 + 236, y + 130 + dy * 0.4, RED, 4, head=True)
text(C2 + 480, y + 100, "sent", 24, fill=BLUE); text(C2 + 18, y + 100, "received", 24, fill=RED, anchor="start")
text(C2 + 280, y + 225, "removed before the pair-specific share", 22, anchor="middle")
bars(C3 + 60, y + 205, 540, 120, [100 * g for g in gl], PURPLE, SHORT, 60, fmt="{:.0f}%", ylab="gene-level")

# ------------------------------------------------ row 5: prediction from co-expression
y = row(4, RED, "Prediction of the response", ["Co-expression in control", "cells, or a model trained", "on other knockdowns"], "Figs. 6, 7")
for k in range(8):
    cell(C2 + 30 + (k % 4) * 30, y + 100 + (k // 4) * 30, GREY)
text(C2 + 75, y + 175, "control cells", 23, anchor="middle")
line(C2 + 160, y + 115, C2 + 230, y + 115, GREY, 4, head=True)
matrix(C2 + 250, y + 62, 4, 28, lambda i, j: None if i == j else ((RED, sym[(i, j)]) if (i, j) in sym else None))
text(C2 + 306, y + 200, "predicted", 23, anchor="middle", fill=RED); text(C2 + 306, y + 228, "response", 23, anchor="middle", fill=RED)
text(C2 + 395, y + 124, "vs", 26, anchor="middle")
matrix(C2 + 430, y + 62, 4, 28, lambda i, j: None if i == j else full.get((i, j)))
text(C2 + 486, y + 200, "second", 23, anchor="middle"); text(C2 + 486, y + 228, "screen", 23, anchor="middle")
bx, by, bw, bh = C3 + 60, y + 205, 540, 120
bars(bx, by, bw, bh, fit, RED, SHORT, 0.8, ylab="correlation", lift=[max(0, (c - f) / 0.8 * bh + 6) if 0 <= (c - f) / 0.8 * bh < 30 or c < f else 0 for f, c in zip(fit, ceil)])
for i, c in enumerate(ceil):
    x0 = bx + i * bw / len(ceil)
    line(x0 + 8, by - c / 0.8 * bh, x0 + bw / len(ceil) - 8, by - c / 0.8 * bh, INK, 4)
text(bx + bw, by + 58, "black line: second screen", 21, anchor="end")

H = Y0 + 5 * RH + 10
out = "\n".join(svg).replace("@H@", str(H)).replace("@HMM@", f"{H / 10:.0f}") + "</svg>"
open(F + "fig1_study_design.svg", "w").write(out)
subprocess.run(["rsvg-convert", "-d", "600", "-p", "600", F + "fig1_study_design.svg", "-o", F + "fig1_study_design.png"], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", F + "fig1_study_design.svg", "-o", F + "fig1_study_design.pdf"], check=True)
print("figure 1 written")
