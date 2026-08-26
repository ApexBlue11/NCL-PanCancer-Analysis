"""Shared figure conventions and a genuine legend-collision test.

Colour carries one meaning across every figure in this manuscript:

    red   a statistically significant association in the positive direction
          (higher in tumour, adverse hazard, positive enrichment)
    blue  a statistically significant association in the negative direction
    grey  not significant at q<0.05

Because the same encoding recurs in Figures 1a, 1c, 2a, 5 and 6, each of those
panels carries its own legend rather than relying on the reader to carry one
across pages.

`audit_legends` is the reason this module exists. Placing a legend with
loc="best" only avoids the *lines* matplotlib knows about; it happily lands on
scatter points and bar segments. The audit transforms every data mark to
display coordinates and reports any that fall under a legend box, so a
collision fails loudly instead of shipping.
"""
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

UP, DOWN, NS, GREY = "#B2182B", "#2166AC", "#BBBBBB", "#666666"
ACCENT = "#1B7837"


def sig_handles(pos_label, neg_label, ns_label="Not significant", extra=None):
    """Legend handles for the shared red / blue / grey significance encoding."""
    h = [Line2D([], [], marker="o", ls="", color=UP, label=pos_label),
         Line2D([], [], marker="o", ls="", color=DOWN, label=neg_label),
         Line2D([], [], marker="o", ls="", color=NS, label=ns_label)]
    if extra:
        h.extend(extra)
    return h


def bar_handles(pos_label, neg_label, ns_label=None):
    """Legend handles drawn as filled swatches, for bar-based panels."""
    h = [Rectangle((0, 0), 1, 1, fc=UP, ec="none", label=pos_label),
         Rectangle((0, 0), 1, 1, fc=DOWN, ec="none", label=neg_label)]
    if ns_label:
        h.append(Rectangle((0, 0), 1, 1, fc=NS, ec="none", label=ns_label))
    return h


def _legend(ax, handles, labels, kw):
    """Call ax.legend correctly whether or not handles/labels were supplied."""
    if handles is not None and labels is not None:
        return ax.legend(handles=handles, labels=labels, **kw)
    if handles is not None:
        return ax.legend(handles=handles, **kw)
    return ax.legend(**kw)


def legend_below(ax, handles=None, labels=None, ncol=3, y=None, fontsize=5.6,
                 pad=0.035):
    """Legend beneath the axes, clear of the tick labels and the axis label.

    The vertical offset is measured rather than guessed: a hand-picked constant
    that clears the x-axis label in one panel will collide with it in another
    whose tick labels are taller or whose axes box is shorter. Pass `y` to
    override.
    """
    if y is None:
        fig = ax.figure
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        ab = ax.get_window_extent(r)
        lows = [ab.y0]
        if ax.xaxis.label.get_text():
            lows.append(ax.xaxis.label.get_window_extent(r).y0)
        for t in ax.get_xticklabels():
            if t.get_text():
                lows.append(t.get_window_extent(r).y0)
        # convert the lowest occupied point into axes fraction, then step below
        y = (min(lows) - ab.y0) / ab.height - pad
    return _legend(ax, handles, labels,
                   dict(loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol,
                        frameon=False, fontsize=fontsize, handletextpad=0.5,
                        columnspacing=1.4, borderaxespad=0.0))


def legend_right(ax, handles=None, labels=None, y=1.0, fontsize=5.6):
    """Legend to the right of the axes."""
    return _legend(ax, handles, labels,
                   dict(loc="upper left", bbox_to_anchor=(1.01, y), frameon=False,
                        fontsize=fontsize, handletextpad=0.5, borderaxespad=0.0))


def _mark_boxes(ax, renderer):
    """Display-coordinate boxes for the individual data marks on an axes.

    Scatter collections are expanded to their per-point offsets: asking a
    PathCollection for its window extent returns one box spanning every point,
    which would flag any legend anywhere inside the axes.
    """
    boxes = []
    for coll in ax.collections:
        try:
            offs = coll.get_offsets()
        except Exception:
            continue
        if offs is None or len(offs) == 0:
            continue
        pts = ax.transData.transform(np.asarray(offs, float))
        for x, y in pts:
            if np.isfinite(x) and np.isfinite(y):
                boxes.append((x - 3, y - 3, x + 3, y + 3))
    for patch in ax.patches:
        try:
            b = patch.get_window_extent(renderer)
        except Exception:
            continue
        if b.width > 0 and b.height > 0:
            boxes.append((b.x0, b.y0, b.x1, b.y1))
    for line in ax.lines:
        xy = line.get_xydata()
        if xy is None or len(xy) == 0:
            continue
        # Reference rules (a single horizontal or vertical line) are not data.
        if len(xy) == 2 and (xy[0][0] == xy[1][0] or xy[0][1] == xy[1][1]):
            continue
        pts = ax.transData.transform(np.asarray(xy, float))
        for x, y in pts:
            if np.isfinite(x) and np.isfinite(y):
                boxes.append((x - 2, y - 2, x + 2, y + 2))
    return boxes


def _text_boxes(ax, renderer):
    """Display boxes for the axis labels, tick labels and title of an axes."""
    out = []
    for label, what in ((ax.xaxis.label, "x-axis label"),
                        (ax.yaxis.label, "y-axis label"),
                        (ax.title, "panel title")):
        if label.get_text():
            b = label.get_window_extent(renderer)
            if b.width > 0 and b.height > 0:
                out.append((b, what))
    for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        if not t.get_text():
            continue
        b = t.get_window_extent(renderer)
        if b.width > 0 and b.height > 0:
            out.append((b, "tick label"))
    return out


def audit_legends(fig, figname):
    """Report legends overlapping data marks or axis text.

    Both matter: a legend sitting on the x-axis label is as unusable as one
    sitting on the points, and only checking data marks misses it.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    problems = []
    for i, ax in enumerate(fig.axes):
        lg = ax.get_legend()
        if lg is None:
            continue
        lb = lg.get_window_extent(renderer)
        panel = (ax.get_title() or "axes %d" % i).strip() or "axes %d" % i

        hits = sum(1 for x0, y0, x1, y1 in _mark_boxes(ax, renderer)
                   if not (x1 < lb.x0 or x0 > lb.x1 or y1 < lb.y0 or y0 > lb.y1))
        if hits:
            problems.append("%s panel '%s': legend covers %d data mark(s)"
                            % (figname, panel, hits))

        for b, what in _text_boxes(ax, renderer):
            if not (b.x1 < lb.x0 or b.x0 > lb.x1 or b.y1 < lb.y0 or b.y0 > lb.y1):
                problems.append("%s panel '%s': legend overlaps the %s"
                                % (figname, panel, what))
    return problems
