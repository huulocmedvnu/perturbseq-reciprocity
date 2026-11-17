"""Figures 2 to 7 from the tables in results/, with the source data of every panel.

Usage: python3 scripts/plot_figures.py
Writes figures/fig*.pdf and .png and results/figure_source_data/*.csv.
Style: plotly_white, Plotly colours, Arial, black text, no gridlines. Figures are at most 820 px wide with
text of at least 12 px, which gives 7 pt at a printed width of 170 mm.

Author: Thai Huu Loc Tran
"""
import json, os, glob, math, numpy as np, pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
pio.kaleido.scope.mathjax = None
from plotly.subplots import make_subplots

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R, F, SD = os.path.join(ROOT, "results/"), os.path.join(ROOT, "figures/"), os.path.join(ROOT, "results/figure_source_data/")
os.makedirs(F, exist_ok=True); os.makedirs(SD, exist_ok=True)
C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692"]; K = "#111111"
LINES = ["k562ess", "rpe1", "hepg2", "jurkat", "hct116", "hek293t"]
NAME = {"k562ess": "K562", "rpe1": "RPE1", "hepg2": "HepG2", "jurkat": "Jurkat", "hct116": "HCT116", "hek293t": "HEK293T"}


def runs(name): return json.load(open(R + name + ".json"))["runs"]
def get(name, s, floor=0.25, kind="scaled"): return next(r for r in runs(name) if r["set"] == s and r["umi_floor"] == floor and r["matrix"] == kind)
def source(panel, df): pd.DataFrame(df).to_csv(SD + panel + ".csv", index=False)
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def style(fig, w, h):
    fig.update_layout(template="plotly_white", width=w, height=h, font=dict(family="Arial", size=13, color=K), margin=dict(l=70, r=20, t=40, b=60), legend=dict(font=dict(color=K, size=12)))
    fig.update_xaxes(showgrid=False, showline=True, linecolor=K, ticks="outside", tickfont=dict(color=K, size=12), title_font=dict(color=K), zeroline=False)
    fig.update_yaxes(showgrid=False, showline=True, linecolor=K, ticks="outside", tickfont=dict(color=K, size=12), title_font=dict(color=K), zeroline=False)
    for a in fig.layout.annotations:
        a.font = dict(family="Arial", size=14 if a.text[:2] in ("a ", "b ", "c ", "d ") else 12, color=K)
    return fig
def save(fig, name): fig.write_image(F + name + ".png", scale=4); fig.write_image(F + name + ".pdf")
def err(lo, hi, v): return dict(type="data", symmetric=False, array=[h - x for h, x in zip(hi, v)], arrayminus=[x - l for l, x in zip(lo, v)], color=K, thickness=1, width=3)


# ---------------------------------------------------------------- Figure 2: two independent K562 screens
d = pd.read_csv(R + "k562_two_screens.density.tsv", sep="\t")
d = d[(d.umi_floor == 0.25) & (d.matrix == "scaled") & (d.set == "active_pairs")]
cal = pd.read_csv(R + "calibration.tsv", sep="\t")
fig = make_subplots(rows=2, cols=2, subplot_titles=["a  Symmetric part", "b  Skew-symmetric part", "c  Simulated screens", "d  Skew-symmetric share, K562"], horizontal_spacing=0.16, vertical_spacing=0.17)
for col, part in ((1, "symmetric"), (2, "skew")):
    p = d[d.part == part].copy()
    lim = p.x_hi.max()
    p = p[(p.x_lo > -lim + 1e-9) & (p.x_hi < lim - 1e-9) & (p.y_lo > p.y_lo.min() + 1e-9) & (p.y_hi < p.y_hi.max() - 1e-9)]
    xs = np.sort(p.x_lo.unique()); w = xs[1] - xs[0]
    z = np.full((60, 60), np.nan)
    for r in p.itertuples(): z[int(round((r.y_lo + lim) / w)), int(round((r.x_lo + lim) / w))] = math.log10(r.count)
    ax = -lim + w * (np.arange(60) + 0.5)
    fig.add_trace(go.Heatmap(x=ax, y=ax, z=z, colorscale="Viridis", zmin=0, zmax=np.nanmax(z), showscale=(col == 2),
                             colorbar=dict(title=dict(text="Pairs", font=dict(color=K, size=12)), tickvals=[0, 1, 2, 3, 4], ticktext=["1", "10", "100", "1,000", "10,000"], len=0.42, y=0.80, x=1.0, tickfont=dict(color=K, size=12))), 1, col)
    fig.add_trace(go.Scatter(x=[-lim, lim], y=[-lim, lim], mode="lines", line=dict(color="#888888", dash="dash", width=1), showlegend=False), 1, col)
    rr = get("k562_two_screens", "active_pairs")
    fig.add_annotation(x=-0.9 * lim, y=0.88 * lim, xref=f"x{col if col > 1 else ''}", yref=f"y{col if col > 1 else ''}", showarrow=False, xanchor="left", text=f"r = {rr['r_S' if part == 'symmetric' else 'r_A']:.2f}")
    fig.update_xaxes(title_text="Genome-wide screen, held-out cells", range=[-lim, lim], row=1, col=col); fig.update_yaxes(title_text="Essential-gene screen", range=[-lim, lim], row=1, col=col)
    source(f"fig2{'ab'[col - 1]}_{part}_density", p[["x_lo", "x_hi", "y_lo", "y_hi", "count"]])
fig.add_trace(go.Scatter(x=[0, 0.55], y=[0, 0.55], mode="lines", line=dict(color="#888888", dash="dash", width=1), showlegend=False), 2, 1)
fig.add_trace(go.Scatter(x=cal.true_skew_share, y=cal.single_screen_share_mean, mode="markers+lines", marker=dict(color=C[1], size=9, symbol="square"), line=dict(color=C[1], width=1), name="Within one screen", legend="legend"), 2, 1)
fig.add_trace(go.Scatter(x=cal.true_skew_share, y=cal.cross_screen_estimate_mean, mode="markers", marker=dict(color=C[0], size=10, symbol="diamond"), name="Between two screens", legend="legend"), 2, 1)
fig.update_xaxes(title_text="True skew-symmetric share", range=[-0.03, 0.56], row=2, col=1); fig.update_yaxes(title_text="Estimated share", range=[-0.03, 0.56], row=2, col=1)
source("fig2c_calibration", cal)
sets = [("All pairs", "all_pairs", 0.25, "scaled"), ("Active pairs", "active_pairs", 0.25, "scaled"), ("Active, neighbours<br>excluded", "active_pairs_without_neighbours", 0.25, "scaled"),
        ("Active, unscaled<br>response", "active_pairs", 0.25, "raw"), ("Active, 0.5 UMI<br>floor", "active_pairs", 0.5, "scaled")]
rs = [get("k562_two_screens", s, fl, kd) for _, s, fl, kd in sets]
v = [r["f_asym"] for r in rs]
fig.add_trace(go.Scatter(y=[s[0] for s in sets], x=v, mode="markers", marker=dict(color=C[0], size=10, symbol="diamond"), error_x=err([r["f_asym_ci"][0] for r in rs], [r["f_asym_ci"][1] for r in rs], v), showlegend=False), 2, 2)
fig.add_vline(x=0.5, line=dict(color="#888888", dash="dash", width=1), row=2, col=2)
fig.update_xaxes(title_text="Skew-symmetric share", range=[0, 0.6], row=2, col=2); fig.update_yaxes(autorange="reversed", row=2, col=2)
source("fig2d_share_k562", pd.DataFrame({"set": [s[0].replace("<br>", " ") for s in sets], "pairs": [r["n_pairs"] for r in rs], "share": v, "ci_low": [r["f_asym_ci"][0] for r in rs], "ci_high": [r["f_asym_ci"][1] for r in rs]}))
style(fig, 820, 800).update_layout(legend=dict(x=0.19, y=0.03, yanchor="bottom", bgcolor="rgba(255,255,255,0.7)"), margin=dict(l=70, r=70, t=40, b=60))
save(fig, "fig2_k562_two_screens")

# ---------------------------------------------------------------- Figure 3: six cell lines, within and between
within = [("within_" + l, NAME[l]) for l in LINES]
between = [(f"between_{a}_{b}", f"{NAME[a]} and {NAME[b]}") for i, a in enumerate(LINES) for b in LINES[i + 1:]]
between += [(f"between_gwps_{b}", f"K562 and {NAME[b]}, genome-wide") for b in ("hct116", "hek293t")]
comp = [("k562_two_screens", "K562, two screens")] + within + between
rs = [get(f, "active_pairs") for f, _ in comp]
labs = [n for _, n in comp]
fig = make_subplots(rows=1, cols=2, subplot_titles=["a  Agreement between screens", "b  Skew-symmetric share"], horizontal_spacing=0.08, shared_yaxes=True)
for key, nm, c, sym in (("r_S", "Symmetric part", C[2], "circle"), ("r_A", "Skew-symmetric part", C[0], "diamond")):
    v = [r[key] for r in rs]
    fig.add_trace(go.Scatter(y=labs, x=v, mode="markers", marker=dict(color=c, size=10, symbol=sym), name=nm, error_x=err([r[key + "_ci"][0] for r in rs], [r[key + "_ci"][1] for r in rs], v)), 1, 1)
v = [r["f_asym"] for r in rs]
fig.add_trace(go.Scatter(y=labs, x=v, mode="markers", marker=dict(color=C[0], size=10, symbol="diamond"), showlegend=False, error_x=err([r["f_asym_ci"][0] for r in rs], [r["f_asym_ci"][1] for r in rs], v)), 1, 2)
fig.add_vline(x=0.5, line=dict(color="#888888", dash="dash", width=1), row=1, col=2)
for y0 in (0.5, 0.5 + len(within), 0.5 + len(within) + 15): fig.add_hline(y=y0, line=dict(color="#BBBBBB", width=1))
fig.update_xaxes(title_text="Pearson correlation", range=[0, 1], row=1, col=1); fig.update_xaxes(title_text="Skew-symmetric share", range=[0, 0.6], row=1, col=2)
fig.update_yaxes(autorange="reversed")
style(fig, 820, 860).update_layout(legend=dict(x=0.20, y=0.01, yanchor="bottom", bgcolor="rgba(255,255,255,0.7)"), margin=dict(l=215, r=20, t=40, b=60))
save(fig, "fig3_cell_lines")
source("fig3_cell_lines", pd.DataFrame({"comparison": labs, "genes": [r["n_genes"] for r in rs], "active_pairs": [r["n_pairs"] for r in rs], "r_symmetric": [r["r_S"] for r in rs], "r_symmetric_ci_low": [r["r_S_ci"][0] for r in rs], "r_symmetric_ci_high": [r["r_S_ci"][1] for r in rs],
                                        "r_skew": [r["r_A"] for r in rs], "r_skew_ci_low": [r["r_A_ci"][0] for r in rs], "r_skew_ci_high": [r["r_A_ci"][1] for r in rs], "skew_share": v, "share_ci_low": [r["f_asym_ci"][0] for r in rs], "share_ci_high": [r["f_asym_ci"][1] for r in rs]}))

# ---------------------------------------------------------------- Figure 4: pair classes, direction, counts
five = [("k562_two_screens", "K562,<br>two screens")] + [("within_" + l, NAME[l]) for l in LINES]
fig = make_subplots(rows=2, cols=2, specs=[[{"colspan": 2}, None], [{}, {}]], subplot_titles=["a  Skew-symmetric share by class of gene pair", "b  Transcription factor and target", "c  Pairs with a replicated effect"], horizontal_spacing=0.17, vertical_spacing=0.20)
x = [n for _, n in five]; src = []
for (lab, fsuf, s), c in zip((("Active pairs outside complexes", "", "active_pairs_outside_complexes"), ("Subunits of one complex, CORUM", "", "complex_active"), ("Subunits of one complex, six named complexes", "_alt", "complex_active")), (C[0], C[4], C[3])):
    rs = [get(f + fsuf, s) for f, _ in five]; v = [r["f_asym"] for r in rs]
    fig.add_trace(go.Bar(x=x, y=v, name=lab, marker_color=c, legend="legend", error_y=err([r["f_asym_ci"][0] for r in rs], [r["f_asym_ci"][1] for r in rs], v)), 1, 1)
    src += [dict(screen=n.replace("<br>", " "), pair_class=lab, pairs=r["n_pairs"], share=r["f_asym"], ci_low=r["f_asym_ci"][0], ci_high=r["f_asym_ci"][1]) for (_, n), r in zip(five, rs)]
fig.update_yaxes(title_text="Skew-symmetric share", range=[0, 0.72], row=1, col=1); source("fig4a_pair_classes", src)
src = []
for (lab, fsuf), c in zip((("DoRothEA A, B and TRRUST", ""), ("CollecTRI", "_alt")), (C[0], C[2])):
    k = n = None; ks, ns = [], []
    for f, _ in five:
        t = [get(f + fsuf, "tf_direction_select_screen1_read_screen2"), get(f + fsuf, "tf_direction_select_screen2_read_screen1")]
        ks.append(sum(r["concordant_active"] for r in t)); ns.append(sum(r["n_active"] for r in t))
    v = [100 * a / b for a, b in zip(ks, ns)]; ci = [wilson(a, b) for a, b in zip(ks, ns)]
    fig.add_trace(go.Bar(x=x, y=v, name=lab, marker_color=c, legend="legend2", error_y=err([100 * a for a, _ in ci], [100 * b for _, b in ci], v)), 2, 1)
    src += [dict(screen=nm.replace("<br>", " "), database=lab, pairs_read_in_both_directions_of_selection=b, larger_effect_from_factor_to_target=a, percent=p) for (_, nm), a, b, p in zip(five, ks, ns, v)]
fig.add_hline(y=50, line=dict(color="#888888", dash="dash", width=1), row=2, col=1)
fig.update_yaxes(title_text="Larger effect from factor to target (%)", range=[0, 150], tickvals=[0, 25, 50, 75, 100], row=2, col=1); source("fig4b_tf_direction", src)
dp = pd.concat([pd.read_csv(R + f + ".directed.counts.tsv", sep="\t") for f, _ in five]); tot = dp.pairs_with_replicated_effect
base = np.zeros(len(dp))
for col, lab, c in (("both_directions", "Both directions", C[1]), ("one_direction_reverse_intermediate", "One direction, reverse uncertain", C[5]), ("one_direction_reverse_absent", "One direction, reverse absent", C[0])):
    h = 100 * dp[col].values / tot.values
    fig.add_trace(go.Bar(x=x, y=h, base=base, name=lab, marker_color=c, legend="legend3", offsetgroup="pairs", width=0.6), 2, 2)
    base = base + h
fig.update_yaxes(title_text="Pairs (%)", range=[0, 150], tickvals=[0, 25, 50, 75, 100], row=2, col=2)
dp.insert(0, "screen", [n.replace("<br>", " ") for _, n in five]); source("fig4c_directed_pairs", dp.drop(columns=["screen1", "screen2"]))
style(fig, 820, 900).update_layout(barmode="group", legend=dict(x=0.01, y=1.0, bgcolor="rgba(255,255,255,0.7)"), legend2=dict(x=0.01, y=0.385, yanchor="top", font=dict(color=K, size=12), bgcolor="rgba(255,255,255,0.7)"),
                                   legend3=dict(x=0.60, y=0.385, yanchor="top", font=dict(color=K, size=12), bgcolor="rgba(255,255,255,0.7)"), margin=dict(l=70, r=20, t=40, b=70))
save(fig, "fig4_pair_classes")

# ---------------------------------------------------------------- Figure 5: gene level
g = pd.read_csv(R + "k562_two_screens.gene_level.genes.tsv", sep="\t"); g = g[g.matrix == "scaled"]
summ = {f: pd.read_csv(R + f + ".gene_level.summary.tsv", sep="\t", nrows=2) for f, _ in five}
fig = make_subplots(rows=1, cols=2, subplot_titles=["a  K562, effect sent and received by each gene", "b  Gene-level and pair-specific parts"], horizontal_spacing=0.14, column_widths=[0.5, 0.5])
pos = g[(g.replicable_effect_sent > 0) & (g.replicable_effect_received > 0)]
fig.add_trace(go.Scatter(x=pos.replicable_effect_received, y=pos.replicable_effect_sent, mode="markers", marker=dict(color=C[0], size=4, opacity=0.5), showlegend=False), 1, 1)
top = pd.concat([pos.sort_values("replicable_effect_sent", ascending=False).head(2), pos.sort_values("replicable_effect_received", ascending=False).head(2)])
for k, r in enumerate(top.itertuples()):
    fig.add_annotation(x=math.log10(r.replicable_effect_received), y=math.log10(r.replicable_effect_sent), text=r.gene, showarrow=True, arrowhead=0, arrowwidth=1, arrowcolor=K, ax=[-45, 45, -30, 10][k], ay=[-14, -22, -40, -45][k], xref="x", yref="y")
tv = [0.01, 0.1, 1, 10, 100, 1000]; tt = ["0.01", "0.1", "1", "10", "100", "1,000"]
fig.update_xaxes(type="log", title_text="Effect received as a readout", tickvals=tv, ticktext=tt, row=1, col=1); fig.update_yaxes(type="log", title_text="Effect sent as a knockdown", tickvals=tv, ticktext=tt, row=1, col=1)
source("fig5a_sent_received", g[["gene", "replicable_effect_sent", "replicable_effect_received", "gene_level_term_screen1", "gene_level_term_screen2"]])
sc = [summ[f][summ[f].matrix == "scaled"].iloc[0] for f, _ in five]
fig.add_trace(go.Bar(x=x, y=[100 * float(s.share_of_V_asym_gene_level) for s in sc], name="Gene-level share of the<br>skew-symmetric part (%)", marker_color=C[4]), 1, 2)
fig.add_trace(go.Bar(x=x, y=[100 * float(s.f_asym_pair_specific) for s in sc], name="Skew-symmetric share after<br>removal of gene-level terms (%)", marker_color=C[0]), 1, 2)
fig.update_yaxes(title_text="Percent", range=[0, 82], row=1, col=2)
source("fig5b_gene_level", pd.DataFrame([dict(screen=n.replace("<br>", " "), **s.to_dict()) for (_, n), s in zip(five, sc)]))
style(fig, 820, 470).update_layout(barmode="group", legend=dict(x=0.56, y=1.0, bgcolor="rgba(255,255,255,0.7)"), margin=dict(l=70, r=20, t=40, b=70))
save(fig, "fig5_gene_level")

# ---------------------------------------------------------------- Figure 6: linear-response prediction
lr = pd.read_csv(R + "linear_response.tsv", sep="\t")
lr["screen"] = [("K562,<br>two screens" if "gwps" in b else NAME[a.split("/")[-1]]) for a, b in zip(lr.screen1, lr.screen2)]
fig = make_subplots(rows=1, cols=3, subplot_titles=["a  Whole response", "b  Symmetric part", "c  Skew-symmetric part"], horizontal_spacing=0.07, shared_yaxes=True)
for col, (pc, cc) in enumerate((("r_prediction", "r_ceiling"), ("r_symmetric_part", "r_ceiling_symmetric"), ("r_skew_part", "r_ceiling_skew")), 1):
    for pred, lab, c in (("symmetric", "Co-expression, one scale", C[2]), ("on_target_scale", "Scaled by on-target knockdown", C[4]), ("fitted_scale", "One fitted scale per knockdown", C[0])):
        s = lr[lr.predictor == pred]
        fig.add_trace(go.Bar(x=s.screen, y=s[pc], name=lab, marker_color=c, showlegend=(col == 1)), 1, col)
    s = lr[lr.predictor == "symmetric"]
    fig.add_trace(go.Scatter(x=s.screen, y=s[cc], mode="markers", marker=dict(color=K, size=11, symbol="line-ew", line=dict(color=K, width=3)), name="Second screen (ceiling)", showlegend=(col == 1)), 1, col)
    fig.update_xaxes(tickangle=-90, row=1, col=col)
fig.update_yaxes(title_text="Correlation with the second screen", row=1, col=1); fig.update_yaxes(range=[-0.3, 0.8])
fig.add_hline(y=0, line=dict(color=K, width=1))
style(fig, 820, 600).update_layout(barmode="group", legend=dict(orientation="h", x=0.0, y=-0.36, yanchor="top"), margin=dict(l=70, r=20, t=40, b=190))
save(fig, "fig6_linear_response")
source("fig6_linear_response", lr.drop(columns=["screen1", "screen2"]).assign(screen=lr.screen.str.replace("<br>", " ")))

# ---------------------------------------------------------------- Figure 7: held-out knockdowns
ho = pd.read_csv(R + "heldout_prediction.tsv", sep="\t")
ho["screen"] = [("K562,<br>two screens" if "gwps" in b else NAME[a.split("/")[-1]]) for a, b in zip(ho.screen1, ho.screen2)]
std, raw = ho[ho.response == "standardised"], ho[ho.response == "unnormalised_counts"]
fig = make_subplots(rows=1, cols=3, subplot_titles=["a  Whole response", "b  Skew-symmetric part", "c  Unnormalised counts"], horizontal_spacing=0.07, column_widths=[0.36, 0.36, 0.28])
P7 = (("coexpression_global_scale", "Co-expression, one scale", C[2]), ("response_scale_from_screen2", "Co-expression, scale fitted<br>to the scored response", C[0]),
      ("mean_expression_scale_from_screen2", "Mean expression, scale fitted<br>to the scored response", C[5]), ("training_mean", "Mean of training knockdowns", C[4]),
      ("stacked_all", "Model trained on other knockdowns", C[1]), ("stacked_with_reference_lines", "With the response in other cell lines", C[3]))
for col, (d, y) in enumerate(((std, "r_whole"), (std, "r_skew_part"), (raw, "r_whole")), 1):
    for pred, lab, c in P7:
        if (pred == "mean_expression_scale_from_screen2") != (col == 3) and pred in ("mean_expression_scale_from_screen2", "stacked_with_reference_lines"): continue
        s = d[d.predictor == pred]
        fig.add_trace(go.Bar(x=s.screen, y=s[y], name=lab, marker_color=c, legendgroup=pred, showlegend=(col == 1 or pred == "mean_expression_scale_from_screen2")), 1, col)
    s = d[d.predictor == "ceiling_screen1"]
    fig.add_trace(go.Scatter(x=s.screen, y=s[y], mode="markers", marker=dict(color=K, size=11, symbol="line-ew", line=dict(color=K, width=3)), name="Second screen (ceiling)", showlegend=(col == 1)), 1, col)
    fig.update_xaxes(tickangle=-60, row=1, col=col); fig.update_yaxes(range=[-0.05, 1.0], row=1, col=col)
fig.update_yaxes(title_text="Correlation with the second screen", row=1, col=1)
style(fig, 820, 640).update_layout(barmode="group", bargap=0.18, legend=dict(orientation="h", x=0.0, y=-0.33, yanchor="top"), margin=dict(l=70, r=20, t=40, b=230))
save(fig, "fig7_heldout_prediction")
source("fig7_heldout_prediction", ho.drop(columns=["screen1", "screen2", "fitted"]).assign(screen=ho.screen.str.replace("<br>", " ")))
print("figures written to", os.path.abspath(F))
