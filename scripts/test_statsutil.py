"""Verify statsutil against brute-force / known references before it is used."""
import numpy as np
from scipy import stats
import statsutil as S

rng = np.random.default_rng(7)
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}")


print("cliffs_delta vs brute force")
for (n1, n2, shift) in [(30, 25, 0.0), (40, 60, 1.0), (15, 15, -0.7), (50, 50, 0.3)]:
    a = rng.normal(shift, 1, n1); b = rng.normal(0, 1, n2)
    brute = (np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])) / (n1 * n2)
    d, lo, hi = S.cliffs_delta(a, b, n_boot=400)
    check(f"n=({n1},{n2}) shift={shift}", abs(d - brute) < 1e-9,
          f"delta={d:+.4f} brute={brute:+.4f} CI=[{lo:+.3f},{hi:+.3f}]")

print("cliffs_delta handles ties")
a = np.array([1, 1, 2, 2, 3.]); b = np.array([1, 2, 2, 3, 3.])
brute = (np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])) / 25
d, _, _ = S.cliffs_delta(a, b, n_boot=200)
check("tied data", abs(d - brute) < 1e-9, f"delta={d:+.4f} brute={brute:+.4f}")

print("hedges_g against a hand-computed case")
a = rng.normal(1, 1, 200); b = rng.normal(0, 1, 200)
g, lo, hi = S.hedges_g(a, b)
n1, n2 = 200, 200
sp = np.sqrt(((n1-1)*a.var(ddof=1) + (n2-1)*b.var(ddof=1)) / (n1+n2-2))
expect = (a.mean()-b.mean())/sp * (1 - 3/(4*(n1+n2)-9))
check("g matches formula", abs(g-expect) < 1e-12, f"g={g:.4f}")
check("CI brackets g", lo < g < hi, f"[{lo:.3f},{hi:.3f}]")

print("bh_fdr against statsmodels")
from statsmodels.stats.multitest import multipletests
p = np.concatenate([rng.uniform(0, 1, 60), rng.uniform(0, 1e-4, 15)])
mine = S.bh_fdr(p)
ref = multipletests(p, method="fdr_bh")[1]
check("matches statsmodels", np.allclose(mine, ref), f"maxdiff={np.abs(mine-ref).max():.2e}")
pn = np.array([0.01, np.nan, 0.2, 0.001])
check("NaN passthrough", np.isnan(S.bh_fdr(pn)[1]) and np.isfinite(S.bh_fdr(pn)[0]))

print("jonckheere_terpstra")
# Strong increasing trend must give a large positive z and tiny p.
g_inc = [rng.normal(i * 0.9, 1, 40) for i in range(4)]
_, z, p = S.jonckheere_terpstra(g_inc)
check("increasing -> z>0, p small", z > 5 and p < 1e-6, f"z={z:.2f} p={p:.2e}")
g_dec = [rng.normal(-i * 0.9, 1, 40) for i in range(4)]
_, z2, p2 = S.jonckheere_terpstra(g_dec)
check("decreasing -> z<0", z2 < -5, f"z={z2:.2f}")
# No trend: p should be uniform-ish, so calibration check over many draws.
zs = []
for s in range(300):
    r2 = np.random.default_rng(1000 + s)
    zs.append(S.jonckheere_terpstra([r2.normal(0, 1, 30) for _ in range(4)])[2])
zs = np.array(zs)
check("null p calibrated", 0.02 < (zs < 0.05).mean() < 0.10,
      f"type-I rate={(zs<0.05).mean():.3f}")
# Non-monotonic (up then down) should not be called a trend.
g_hump = [rng.normal(m, 1, 40) for m in (0, 2, 2, 0)]
_, z3, p3 = S.jonckheere_terpstra(g_hump)
check("hump not a trend", p3 > 0.05, f"z={z3:.2f} p={p3:.3f}")
# Permutation path agrees with the asymptotic one.
_, _, p_perm = S.jonckheere_terpstra(g_inc, n_perm=2000)
check("permutation agrees", p_perm < 0.01, f"p_perm={p_perm:.4f}")

print("partial_spearman")
n = 400
c = rng.normal(0, 1, n)
x = c + rng.normal(0, .5, n)
y = c + rng.normal(0, .5, n)          # correlated only through c
r_raw, _ = stats.spearmanr(x, y)
r_adj, p_adj, nn = S.partial_spearman(x, y, c)
check("confound removed", abs(r_adj) < 0.15 < r_raw, f"raw={r_raw:.3f} adj={r_adj:.3f}")
y2 = x * 2 + rng.normal(0, .3, n)     # genuine association independent of c
r2_adj, p2_adj, _ = S.partial_spearman(x, y2, c)
check("real signal kept", r2_adj > 0.6 and p2_adj < 1e-10, f"adj={r2_adj:.3f}")

print("spearman_ci")
lo, hi = S.spearman_ci(0.5, 100)
check("CI brackets rho", lo < 0.5 < hi, f"[{lo:.3f},{hi:.3f}]")
check("CI narrows with n", (S.spearman_ci(0.5, 1000)[1]-S.spearman_ci(0.5, 1000)[0])
      < (hi-lo))

print("\nALL PASS" if ok else "\nFAILURES PRESENT")
raise SystemExit(0 if ok else 1)
