"""Effect sizes, confidence intervals and trend tests.

Reviewer 2 asked for effect sizes with confidence intervals and for multiple
testing correction; the submitted manuscript reported neither. Everything the
revision claims should come through these functions so the reporting is uniform.
"""
import numpy as np
from scipy import stats


def hedges_g(a, b):
    """Standardised mean difference with small-sample correction, plus 95% CI."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan, np.nan
    s_pool = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if s_pool == 0:
        return np.nan, np.nan, np.nan
    d = (a.mean() - b.mean()) / s_pool
    J = 1 - 3 / (4 * (n1 + n2) - 9)          # Hedges' bias correction
    g = J * d
    se = np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2 - 2)))
    return g, g - 1.96 * se, g + 1.96 * se


def _delta_point(a, b):
    """Cliff's delta from the Mann-Whitney U: O(n log n), ties handled by U."""
    n1, n2 = len(a), len(b)
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2 * u / (n1 * n2) - 1


def cliffs_delta(a, b, n_boot=2000, seed=0):
    """Non-parametric effect size in [-1,1] with a bootstrap percentile 95% CI.

    A stratified bootstrap is used instead of Cliff's asymptotic variance so the
    interval needs no distributional assumption -- worth the cost here because
    several TCGA cohorts have very few adjacent normals (PAAD n=4, SKCM n=1).
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan, np.nan
    delta = _delta_point(a, b)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for k in range(n_boot):
        boots[k] = _delta_point(rng.choice(a, n1, replace=True),
                                rng.choice(b, n2, replace=True))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return delta, lo, hi


def spearman_ci(rho, n, alpha=0.05):
    """Fisher z 95% CI for a Spearman correlation."""
    if not np.isfinite(rho) or n < 4:
        return np.nan, np.nan
    rho = np.clip(rho, -0.999999, 0.999999)
    z = np.arctanh(rho)
    se = 1.06 / np.sqrt(n - 3)               # Bonett-Wright SE for Spearman
    crit = stats.norm.ppf(1 - alpha / 2)
    return np.tanh(z - crit * se), np.tanh(z + crit * se)


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values; NaNs pass through untouched."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def _jt_stat(groups):
    """Jonckheere-Terpstra statistic: ordered-pair dominance summed over i<j."""
    total = 0.0
    for i in range(len(groups) - 1):
        for j in range(i + 1, len(groups)):
            total += stats.mannwhitneyu(groups[j], groups[i],
                                        alternative="two-sided").statistic
    return total


def jonckheere_terpstra(groups, n_perm=0, seed=0):
    """Test for a monotonic trend across ordered groups.

    This is the test the submitted manuscript needed and did not run. It claimed
    NCL "progressively increased as the cancer advanced through stages I-IV"
    from UALCAN plots whose asterisks compare each stage against normal tissue,
    not against the preceding stage -- so a stage-ordered trend was never tested.

    The tie-corrected normal approximation is the reported p-value; set n_perm>0
    for a permutation check on small cohorts.

    Returns (JT statistic, standardised z, two-sided p).
    """
    groups = [np.asarray(g, float) for g in groups]
    groups = [g[np.isfinite(g)] for g in groups]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 3:
        return np.nan, np.nan, np.nan

    obs = _jt_stat(groups)
    ns = [len(g) for g in groups]
    N = sum(ns)
    mu = (N ** 2 - sum(n ** 2 for n in ns)) / 4

    # Tie-corrected variance (Lehmann); reduces to the classic form without ties.
    _, tie_counts = np.unique(np.concatenate(groups), return_counts=True)
    t_term = sum(t * (t - 1) * (2 * t + 5) for t in tie_counts)
    n_term = sum(n * (n - 1) * (2 * n + 5) for n in ns)
    var = (N * (N - 1) * (2 * N + 5) - n_term - t_term) / 72
    if var <= 0:
        return obs, np.nan, np.nan
    z = (obs - mu) / np.sqrt(var)
    p = 2 * stats.norm.sf(abs(z))

    if n_perm:
        rng = np.random.default_rng(seed)
        pooled = np.concatenate(groups)
        hits = 0
        for _ in range(n_perm):
            rng.shuffle(pooled)
            idx, perm = 0, []
            for n in ns:
                perm.append(pooled[idx:idx + n]); idx += n
            if abs(_jt_stat(perm) - mu) >= abs(obs - mu):
                hits += 1
        p = (hits + 1) / (n_perm + 1)
    return obs, z, p


def partial_spearman(x, y, covar):
    """Spearman correlation of x and y adjusted for one covariate.

    Used for immune-infiltration correlations adjusted for tumour purity, which
    is how TIMER2 reports them and how they must be interpreted.
    """
    x, y, covar = (np.asarray(v, float) for v in (x, y, covar))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(covar)
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    rx, ry, rc = (stats.rankdata(v[ok]) for v in (x, y, covar))
    resid = []
    for v in (rx, ry):
        A = np.column_stack([np.ones_like(rc), rc])
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        resid.append(v - A @ beta)
    n = int(ok.sum())
    r, _ = stats.pearsonr(resid[0], resid[1])
    # One covariate consumed: t has n-3 degrees of freedom.
    if abs(r) >= 1:
        return r, 0.0, n
    t = r * np.sqrt((n - 3) / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), n - 3)
    return r, p, n
