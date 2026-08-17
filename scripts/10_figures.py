"""Publication figures at 600 dpi.

Colour means the same thing in every panel: red is a significant association in
the positive direction, blue a significant association in the negative
direction, grey not significant (see figstyle). Panels 1a, 1c, 2a, 5 and 6 all
use that encoding and each carries its own legend, so no panel depends on the
reader remembering a key from an earlier figure.

Every legend is checked against the data after drawing (figstyle.audit_legends);
the build fails rather than shipping a legend that sits on top of the plot.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import data_io as D
from figstyle import (UP, DOWN, NS, GREY, ACCENT, sig_handles, bar_handles,
                      legend_below, legend_right, audit_legends)

DPI = 600
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7,
    "axes.linewidth": 0.6, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 110, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})

PROBLEMS = []


def T(name):
    p = os.path.join(D.TABLES, name)
    return pd.read_csv(p, sep="\t", low_memory=False) if os.path.exists(p) else None


def save(fig, name):
    PROBLEMS.extend(audit_legends(fig, name))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(D.FIGURES, f"{name}.{ext}"), dpi=DPI)
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


def sig_colour(delta, q, thr=0.05):
    if not np.isfinite(q) or q >= thr:
        return NS
    return UP if delta > 0 else DOWN


# --------------------------------------------------------------- Figure 1
def figure1():
    de, st = T("T1_differential_expression.tsv"), T("T2_stage_association.tsv")
    if de is None:
        return
    fig = plt.figure(figsize=(7.2, 8.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1], hspace=0.55, wspace=0.34)

    # (a) forest plot vs GTEx
    ax = fig.add_subplot(gs[0, :])
    g = de[(de.comparator == "GTEx") & de.q.notna()].sort_values("cliffs_delta")
    for i, (_, r) in enumerate(g.iterrows()):
        c = sig_colour(r.cliffs_delta, r.q)
        ax.plot([r.delta_lo, r.delta_hi], [i, i], color=c, lw=1.1, zorder=2)
        # One visual channel per variable: colour carries significance and
        # direction. Marker shape previously carried tissue match, which meant
        # the legend had to show a shape in some arbitrary colour; readers took
        # the grey legend square to be a fourth category rather than a shape
        # key. Tissue match now sits on the cancer label instead.
        ax.scatter(r.cliffs_delta, i, s=17, color=c, zorder=3, marker="o",
                   edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="#222222", lw=0.7, ls="--", zorder=1)
    ax.set_yticks(np.arange(len(g)))
    ax.set_yticklabels([f"{r.cohort}{'*' if r.tissue_match != 'good' else ''}"
                        f"  (n={int(r.n_tumour)}/{int(r.n_normal)})"
                        for _, r in g.iterrows()], fontsize=5.6)
    ax.set_xlabel("Cliff's delta, tumour vs GTEx normal (95% CI)")
    ax.set_xlim(-0.75, 1.05)
    ax.set_ylim(-1, len(g))
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
    legend_below(ax, sig_handles(
        "higher in tumour (q<0.05)", "lower in tumour (q<0.05)",
        extra=[Line2D([], [], ls="", marker="",
                      label="* approximate TCGA-GTEx tissue match")]),
        ncol=2, pad=0.075)

    # (b) GTEx vs adjacent-normal comparator
    ax = fig.add_subplot(gs[1, 0])
    a = de[de.comparator == "GTEx"].set_index("cohort")
    b = de[de.comparator == "TCGA_adjacent"].set_index("cohort")
    m = a[["cliffs_delta", "q"]].join(b[["cliffs_delta", "q"]],
                                      lsuffix="_g", rsuffix="_t", how="inner").dropna()
    flip = np.sign(m.cliffs_delta_g) != np.sign(m.cliffs_delta_t)
    both = (m.q_g < 0.05) & (m.q_t < 0.05)
    reverse, partial = flip & both, flip & ~both
    ax.axhline(0, color=GREY, lw=0.5); ax.axvline(0, color=GREY, lw=0.5)
    ax.plot([-1, 1], [-1, 1], color=GREY, lw=0.6, ls=":")
    ax.scatter(m.cliffs_delta_g[~flip], m.cliffs_delta_t[~flip], s=16,
               color=GREY, edgecolor="white", linewidth=0.3, label="same direction")
    ax.scatter(m.cliffs_delta_g[partial], m.cliffs_delta_t[partial], s=22,
               color="#F4A582", edgecolor="white", linewidth=0.4, zorder=3,
               label="sign differs, one comparator ns")
    ax.scatter(m.cliffs_delta_g[reverse], m.cliffs_delta_t[reverse], s=30,
               color=UP, edgecolor="white", linewidth=0.4, zorder=4,
               label="reverses (significant both ways)")
    for c in m.index[flip]:
        ax.annotate(c, (m.cliffs_delta_g[c], m.cliffs_delta_t[c]),
                    fontsize=5.6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Cliff's delta vs GTEx normal")
    ax.set_ylabel("Cliff's delta vs adjacent normal")
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
    legend_below(ax, ncol=1)

    # (c) stage trend
    ax = fig.add_subplot(gs[1, 1])
    if st is not None and len(st):
        s = st.sort_values("trend_z")
        cols = [UP if (r.trend_q < 0.05 and r.trend_z > 0) else
                DOWN if (r.trend_q < 0.05 and r.trend_z < 0) else NS
                for _, r in s.iterrows()]
        ax.barh(np.arange(len(s)), s.trend_z, color=cols, height=0.68)
        ax.axvline(0, color="#222222", lw=0.7)
        for v in (-1.96, 1.96):
            ax.axvline(v, color=GREY, lw=0.5, ls=":")
        ax.set_yticks(np.arange(len(s)))
        ax.set_yticklabels([f"{r.cohort} (n={int(r.n_total)})" for _, r in s.iterrows()],
                           fontsize=5.6)
        ax.set_xlabel("Jonckheere-Terpstra z (stage I to IV trend)")
        ax.set_title("c", loc="left", fontweight="bold", fontsize=10)
        # Same encoding as panel a, restated so the panel stands alone.
        legend_below(ax, bar_handles("increases with stage (q<0.05)",
                                     "decreases with stage (q<0.05)",
                                     "no significant trend"), ncol=1)
    save(fig, "Figure1_expression")


# --------------------------------------------------------------- Figure 2
def figure2():
    sv = T("T3_survival.tsv")
    if sv is None:
        return
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.32, height_ratios=[1.35, 1])

    ax = fig.add_subplot(gs[0, :])
    r = sv[(sv.endpoint == "OS") & sv.adj_HR.notna()].sort_values("adj_HR")
    flagged = False
    for i, (_, x) in enumerate(r.iterrows()):
        c = UP if (x.adj_q < 0.05 and x.adj_HR > 1) else (
            DOWN if (x.adj_q < 0.05 and x.adj_HR < 1) else NS)
        ax.plot([x.adj_HR_lo, x.adj_HR_hi], [i, i], color=c, lw=1.1)
        ax.scatter(x.adj_HR, i, s=18, color=c, zorder=3, edgecolor="white", linewidth=0.3)
        if pd.notna(x.get("ph_p_NCL")) and x.ph_p_NCL < 0.05:
            ax.annotate("†", (x.adj_HR_hi, i), fontsize=6, color=GREY,
                        xytext=(2, -2), textcoords="offset points")
            flagged = True
    ax.axvline(1, color="#222222", lw=0.7, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(r)))
    ax.set_yticklabels([f"{x.cohort} (n={int(x.n_adj)}, e={int(x.events_adj)})"
                        for _, x in r.iterrows()], fontsize=5.6)
    ax.set_xlabel("Adjusted hazard ratio per SD of NCL, overall survival (95% CI)")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
    extra = [Line2D([], [], marker=r"$\dagger$", ls="", color=GREY,
                    label="proportional hazards violated")] if flagged else None
    legend_below(ax, sig_handles("adverse (q<0.05)", "protective (q<0.05)",
                                 extra=extra), ncol=2)

    # (b) univariate vs adjusted
    ax = fig.add_subplot(gs[1, 0])
    d = sv[sv.adj_HR.notna() & sv.uni_HR.notna()]
    ax.plot([0.4, 3.2], [0.4, 3.2], color=GREY, lw=0.6, ls=":")
    for ep, mk in [("OS", "o"), ("DSS", "s"), ("PFS", "^")]:
        s = d[d.endpoint == ep]
        ax.scatter(s.uni_HR, s.adj_HR, s=13, marker=mk, alpha=0.75,
                   edgecolor="white", linewidth=0.25, label=ep)
    ax.axhline(1, color=GREY, lw=0.5); ax.axvline(1, color=GREY, lw=0.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Univariate HR"); ax.set_ylabel("Covariate-adjusted HR")
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
    legend_below(ax, ncol=3)

    # (c) how many survive adjustment
    ax = fig.add_subplot(gs[1, 1])
    eps = ["OS", "DSS", "PFS"]
    uni = [int((sv[(sv.endpoint == e)].uni_q < 0.05).sum()) for e in eps]
    adj = [int((sv[(sv.endpoint == e)].adj_q < 0.05).sum()) for e in eps]
    x = np.arange(3)
    ax.bar(x - 0.19, uni, 0.38, color=NS, label="univariate")
    ax.bar(x + 0.19, adj, 0.38, color=ACCENT, label="covariate-adjusted")
    for i, (u, a) in enumerate(zip(uni, adj)):
        ax.annotate(str(u), (i - 0.19, u), ha="center", va="bottom", fontsize=6)
        ax.annotate(str(a), (i + 0.19, a), ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels(eps)
    ax.set_ylabel("cancers significant (q<0.05)")
    ax.set_ylim(0, max(uni + adj) * 1.25)
    ax.set_title("c", loc="left", fontweight="bold", fontsize=10)
    legend_below(ax, ncol=2)
    save(fig, "Figure2_survival")


# --------------------------------------------------------------- Figure 3
def figure3():
    inf, conc, gs_ = (T("T4_immune_infiltration.tsv"),
                      T("T5_algorithm_concordance.tsv"),
                      T("T8_genomic_scores.tsv"))
    if inf is None:
        return
    fig = plt.figure(figsize=(7.2, 7.4))
    grid = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.34, height_ratios=[1, 1.15])

    # (a) concordance summary per cell type
    ax = fig.add_subplot(grid[0, :])
    if conc is not None and len(conc):
        agg = (conc.groupby("cell_type")
               .agg(concordant=("concordant", "sum"),
                    conflict=("direction", lambda v: (v == "mixed").sum()),
                    total=("cohort", "size")).sort_values("total"))
        y = np.arange(len(agg))
        ax.barh(y, agg.concordant, color=ACCENT, height=0.6, label="concordant")
        ax.barh(y, agg.conflict, left=agg.concordant, color=UP, height=0.6,
                label="algorithms disagree in direction")
        ax.barh(y, agg.total - agg.concordant - agg.conflict,
                left=agg.concordant + agg.conflict, color=NS, height=0.6,
                label="no or insufficient signal")
        ax.set_yticks(y); ax.set_yticklabels(agg.index, fontsize=6)
        ax.set_xlabel("number of cancers")
        ax.set_xlim(0, 33)
        ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
        # Bars span the full width, so the legend must sit outside the axes.
        legend_below(ax, ncol=3)

    # (b) immune / stromal / microenvironment scores
    ax = fig.add_subplot(grid[1, 0])
    if gs_ is not None:
        sub = gs_[gs_.measure.isin(["Immune score", "Stromal score",
                                    "Microenvironment score"])]
        piv = sub.pivot_table(index="cohort", columns="measure", values="rho")
        piv = piv.sort_values("Microenvironment score")
        y = np.arange(len(piv))
        for m, mk in [("Immune score", "o"), ("Stromal score", "s"),
                      ("Microenvironment score", "^")]:
            if m in piv:
                ax.scatter(piv[m], y, s=11, marker=mk, alpha=0.85, label=m,
                           edgecolor="white", linewidth=0.2)
        ax.axvline(0, color="#222222", lw=0.7, ls="--")
        ax.set_yticks(y); ax.set_yticklabels(piv.index, fontsize=5.0)
        ax.set_xlabel("Spearman rho with NCL")
        ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
        # Points reach both edges of the panel; keep the key outside.
        legend_below(ax, ncol=1, pad=0.075)

    # (c) distribution of purity-adjusted rho by algorithm
    ax = fig.add_subplot(grid[1, 1])
    algs = sorted(inf.algorithm.unique())
    data = [inf[inf.algorithm == a].rho_purity_adj.dropna() for a in algs]
    bp = ax.boxplot(data, vert=False, widths=0.6, patch_artist=True,
                    flierprops=dict(marker=".", ms=1, alpha=0.25))
    for p in bp["boxes"]:
        p.set(facecolor="#DDDDDD", edgecolor="#333333", linewidth=0.5)
    for p in bp["medians"]:
        p.set(color=UP, linewidth=1.0)
    ax.axvline(0, color="#222222", lw=0.7, ls="--")
    ax.set_yticklabels(algs, fontsize=6)
    ax.set_xlabel("purity-adjusted Spearman rho with NCL")
    ax.set_title("c", loc="left", fontweight="bold", fontsize=10)
    save(fig, "Figure3_immune")


# --------------------------------------------------------------- Figure 4
def figure4():
    cp = T("T6_checkpoints.tsv")
    if cp is None:
        return
    fig = plt.figure(figsize=(7.2, 6.0))
    grid = fig.add_gridspec(1, 2, wspace=0.42, width_ratios=[1.15, 1])

    ax = fig.add_subplot(grid[0, 0])
    rows = []
    for a, g in cp.groupby("alias"):
        both = (g.q_purity_adj < 0.05) & (g.q_prolif_adj < 0.05)
        rows.append((a, int(((g.rho_prolif_adj > 0) & both).sum()),
                     int(((g.rho_prolif_adj < 0) & both).sum()), len(g)))
    r = pd.DataFrame(rows, columns=["alias", "pos", "neg", "n"]).set_index("alias")
    r = r.sort_values("pos")
    y = np.arange(len(r))
    ax.barh(y, r.pos, color=UP, height=0.62, label="positive")
    ax.barh(y, -r.neg, color=DOWN, height=0.62, label="negative")
    ax.axvline(0, color="#222222", lw=0.7)
    ax.set_yticks(y); ax.set_yticklabels(r.index, fontsize=6.4)
    ax.set_xlabel("cancers with a robust association\n"
                  "(significant after purity AND proliferation adjustment)")
    ax.set_title("a", loc="left", fontweight="bold", fontsize=10)
    legend_below(ax, ncol=2)

    ax = fig.add_subplot(grid[0, 1])
    b = cp[cp.alias == "B7-H3"].sort_values("rho")
    y = np.arange(len(b))
    ax.scatter(b.rho, y, s=13, color=NS, label="unadjusted",
               edgecolor="white", linewidth=0.2)
    ax.scatter(b.rho_purity_adj, y, s=13, color="#7FBC91", label="purity-adjusted",
               edgecolor="white", linewidth=0.2)
    ax.scatter(b.rho_prolif_adj, y, s=15, color=ACCENT, label="proliferation-adjusted",
               edgecolor="white", linewidth=0.2)
    ax.axvline(0, color="#222222", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{x.cohort} (n={int(x.n)})" for _, x in b.iterrows()],
                       fontsize=5.2)
    ax.set_xlabel("Spearman rho, NCL vs CD276 (B7-H3)")
    ax.set_title("b", loc="left", fontweight="bold", fontsize=10)
    legend_below(ax, ncol=1)
    save(fig, "Figure4_checkpoints")


# --------------------------------------------------------------- Figure 5
def figure5():
    cv = T("T7_cptac_validation.tsv")
    if cv is None or "cliffs_delta" not in cv:
        return
    ok = cv[cv.cliffs_delta.notna()].sort_values("cliffs_delta")
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    y = np.arange(len(ok))
    for i, (_, r) in enumerate(ok.iterrows()):
        c = sig_colour(r.cliffs_delta, r.q_ranksum)
        ax.plot([r.delta_lo, r.delta_hi], [i, i], color=c, lw=1.2)
        ax.scatter(r.cliffs_delta, i, s=22, color=c, zorder=3,
                   edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="#222222", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.cohort} ({int(r.n_tumour)}/{int(r.n_normal)}"
                        + (f", {int(r.n_paired)} paired)" if pd.notna(r.n_paired) else ")")
                        for _, r in ok.iterrows()], fontsize=6)
    ax.set_xlabel("Cliff's delta, NCL protein tumour vs normal (95% CI)")
    ax.set_title("CPTAC protein-level validation", fontsize=8, loc="left")
    # Same encoding as Figures 1a and 2a, restated so the panel stands alone.
    legend_below(ax, sig_handles("higher in tumour (q<0.05)",
                                 "lower in tumour (q<0.05)"), ncol=3)
    save(fig, "Figure5_cptac")


# --------------------------------------------------------------- Figure 6
def figure6():
    g = T("T10_gsea_per_cancer.tsv.gz")
    if g is None:
        print("  Figure 6 skipped: GSEA results not present")
        return
    qcol = next((c for c in ["FDR_q-val", "FDR_qval"] if c in g.columns), None)
    if qcol is None:
        return
    g[qcol] = pd.to_numeric(g[qcol], errors="coerce")
    g["NES"] = pd.to_numeric(g["NES"], errors="coerce")
    sig = g[(g[qcol] < 0.05) & g.NES.notna()]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 5.6))
    any_neg = False
    for ax, coll in zip(axes, ["Hallmark", "Reactome"]):
        s = sig[sig.collection == coll]
        if not len(s):
            ax.set_visible(False); continue
        n_c = g[g.collection == coll].cohort.nunique()
        agg = (s.groupby("Term").agg(n=("cohort", "nunique"),
                                     up=("NES", lambda v: (v > 0).sum()),
                                     dn=("NES", lambda v: (v < 0).sum()),
                                     med=("NES", "median")))
        agg = agg[(agg.up == 0) | (agg.dn == 0)].sort_values("n", ascending=False).head(16)
        agg = agg.sort_values("med")
        any_neg |= bool((agg.med < 0).any())
        y = np.arange(len(agg))
        ax.barh(y, agg.med, color=[UP if v > 0 else DOWN for v in agg.med], height=0.65)
        ax.axvline(0, color="#222222", lw=0.7)
        labs = [t.replace(f"{coll.upper()}_", "").replace("_", " ").lower()[:52]
                for t in agg.index]
        ax.set_yticks(y)
        ax.set_yticklabels([f"{l}  ({int(n)}/{n_c})" for l, n in zip(labs, agg.n)],
                           fontsize=5.2)
        ax.set_xlabel("median NES across cancers")
        ax.set_title(coll, fontsize=8, loc="left", fontweight="bold")

    # Red and blue carry the same directional meaning as in the other figures.
    handles = bar_handles("positively enriched with high NCL",
                          "negatively enriched with high NCL")
    if not any_neg:
        handles = [handles[0],
                   Rectangle((0, 0), 1, 1, fc="white", ec="none",
                             label="(no negatively enriched set met the display threshold)")]
    legend_below(axes[0], handles, ncol=1, fontsize=5.6)
    fig.suptitle("Pathways consistently associated with NCL expression",
                 fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    save(fig, "Figure6_gsea")


if __name__ == "__main__":
    print("writing figures ->", D.FIGURES)
    for fn in (figure1, figure2, figure3, figure4, figure5, figure6):
        try:
            fn()
        except Exception as e:
            print(f"  {fn.__name__} FAILED: {type(e).__name__}: {e}")
    print()
    if PROBLEMS:
        print("LEGEND COLLISIONS DETECTED:")
        for p in PROBLEMS:
            print("  " + p)
        raise SystemExit(1)
    print("legend audit: no legend overlaps any data mark")
