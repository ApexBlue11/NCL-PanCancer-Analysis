"""Negative control for the legend-collision audit.

An audit that always passes is worthless. These cases put a legend directly on
top of data and confirm the audit reports it, then move the same legend outside
the axes and confirm the report clears.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from figstyle import audit_legends, legend_below

ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print("  [%s] %s %s" % ("PASS" if cond else "FAIL", name, detail))


rng = np.random.default_rng(0)

# --- scatter: legend dropped in the middle of the cloud must be caught
fig, ax = plt.subplots()
ax.scatter(rng.normal(0, 1, 400), rng.normal(0, 1, 400), s=14, label="points")
ax.legend(loc="center")
p = audit_legends(fig, "scatter-centre")
check("scatter with centred legend is flagged", len(p) == 1, p[0] if p else "")
plt.close(fig)

# --- same data, legend moved below the axes: must come back clean
fig, ax = plt.subplots()
ax.scatter(rng.normal(0, 1, 400), rng.normal(0, 1, 400), s=14, label="points")
legend_below(ax)
p = audit_legends(fig, "scatter-below")
check("same scatter with legend below is clean", len(p) == 0, str(p))
plt.close(fig)

# --- bars spanning the axis: legend inside must be caught
fig, ax = plt.subplots()
ax.barh(np.arange(10), np.full(10, 33), color="grey", label="bars")
ax.set_xlim(0, 33)
ax.legend(loc="lower right")
p = audit_legends(fig, "bars-inside")
check("full-width bars with inside legend flagged", len(p) == 1, p[0] if p else "")
plt.close(fig)

fig, ax = plt.subplots()
ax.barh(np.arange(10), np.full(10, 33), color="grey", label="bars")
ax.set_xlim(0, 33)
legend_below(ax)
p = audit_legends(fig, "bars-below")
check("same bars with legend below is clean", len(p) == 0, str(p))
plt.close(fig)

# --- a legend in genuinely empty space should not be flagged
fig, ax = plt.subplots()
ax.scatter(rng.uniform(0, 0.3, 100), rng.uniform(0, 0.3, 100), s=14, label="corner")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(loc="upper right")
p = audit_legends(fig, "empty-corner")
check("legend in empty region not flagged", len(p) == 0, str(p))
plt.close(fig)

print("\n%s" % ("AUDIT VERIFIED" if ok else "AUDIT IS UNRELIABLE"))
raise SystemExit(0 if ok else 1)
