#!/usr/bin/env python3
"""
BIB revision memo item D (2026-09-13 memo, Tier 2): package the paper's
resolution argument ("~0.10 AUC at 25 events", "~10^3 patients for 0.03-0.05
AUC") as a reusable function + lookup table + worked example, using the
paper's own measured inter-model score correlation (Supplementary
Table tab:score_correlation, r=0.71-0.90) rather than an assumed independence
value -- the memo's point: "the paired case is far more favourable than the
independent case, and your ensemble-score correlations give a real rho
instead of an assumed one."

Method: Hanley-McNeil variance of a single AUC estimate, combined into the
variance of a PAIRED AUC difference via the standard covariance term
(Var(AUC1) + Var(AUC2) - 2*rho*sqrt(Var(AUC1)Var(AUC2))). rho here is the
inter-model score correlation as a proxy for the correlation between the two
models' AUC estimates on the same patients -- a documented approximation,
not the exact DeLong placement-value covariance (which needs the raw scores,
not just a target design), flagged in resolvable_delta_auc's docstring so it
isn't mistaken for the exact per-pair DeLong SE already reported in
Supplementary Table tab:se_check for the ACTUAL observed contrasts.

This is a pure statistics module -- no data access, runs locally.
"""
import math


def _hanley_mcneil_var(auc, n_pos, n_neg):
    """Var(AUC) via Hanley & McNeil (1982)."""
    q1 = auc / (2 - auc)
    q2 = 2 * auc ** 2 / (1 + auc)
    return (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2) + (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg)


def _z(p):
    """Inverse standard normal CDF (no scipy dependency)."""
    # Acklam's algorithm, standard rational approximation, good to ~1e-9.
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    p_low, p_high = 0.02425, 1 - 0.02425
    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    elif p <= p_high:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)


def resolvable_delta_auc(n_pos, n_neg, rho=0.8, alpha=0.05, power=0.8, auc=0.75):
    """Minimum pairwise AUC difference resolvable at the given power, for a
    PAIRED comparison of two models scored on the same n_pos/n_neg patients
    with inter-model score correlation rho (proxy for AUC-estimate
    correlation -- see module docstring; use the measured pairwise value
    from Supplementary Table tab:score_correlation for a specific model
    pair, or a conservative value from that table's range (0.71-0.90) for a
    generic "any two models" estimate).

    auc: the AUC level to evaluate the Hanley-McNeil variance at (the
    variance is fairly flat near 0.7-0.85, where this paper's numbers sit;
    default 0.75 matches the paper's own headline range).
    """
    var1 = _hanley_mcneil_var(auc, n_pos, n_neg)
    var2 = var1  # same n_pos/n_neg, same auc level -> symmetric for a generic two-model comparison
    var_diff = var1 + var2 - 2 * rho * math.sqrt(var1 * var2)
    se_diff = math.sqrt(max(var_diff, 0.0))
    z_alpha = _z(1 - alpha / 2)
    z_power = _z(power)
    return (z_alpha + z_power) * se_diff


def worked_example():
    print('Worked example (memo item D): a 3-cohort, ~60-event study', flush=True)
    # roughly this paper's own COVID event rate (25/240 ~ 10.4%) scaled to 60 events
    n_pos, n_neg = 60, 520
    for rho in (0.71, 0.80, 0.90):
        d = resolvable_delta_auc(n_pos, n_neg, rho=rho)
        print(f'  n_fatal={n_pos}, n_survived={n_neg}, rho={rho}: '
              f'resolvable delta-AUC = {d:.3f} at alpha=0.05, power=0.80', flush=True)


def lookup_table():
    """Events x detectable delta-AUC at 80% power, event rate fixed at this
    paper's own within-COVID rate (25/240 = 10.4%), rho from the measured
    range in tab:score_correlation."""
    print('\nLookup table: events -> resolvable delta-AUC (80% power, alpha=0.05, '
          'event rate fixed at 10.4%, this paper\'s own within-COVID rate)', flush=True)
    print(f'{"n_fatal":>8} {"n_total":>8} {"rho=0.71":>10} {"rho=0.80":>10} {"rho=0.90":>10}', flush=True)
    rows = []
    for n_fatal in (10, 17, 25, 30, 50, 75, 100, 150, 250, 500, 1000):
        n_total = round(n_fatal / 0.104)
        n_neg = n_total - n_fatal
        vals = [resolvable_delta_auc(n_fatal, n_neg, rho=r) for r in (0.71, 0.80, 0.90)]
        rows.append((n_fatal, n_total, *vals))
        marker = ''
        if n_fatal == 25:
            marker = '  <- this paper, within-COVID (240 patients, 25 fatal)'
        elif n_fatal == 17:
            marker = '  <- this paper, within-sepsis (91 patients, 17 fatal)'
        print(f'{n_fatal:>8} {n_total:>8} {vals[0]:>10.3f} {vals[1]:>10.3f} {vals[2]:>10.3f}{marker}', flush=True)
    return rows


if __name__ == '__main__':
    lookup_table()
    print()
    worked_example()
    print()
    print('Sanity check against the paper\'s own stated numbers:', flush=True)
    d25 = resolvable_delta_auc(25, round(240 * 0.896), rho=0.8)
    print(f'  25 fatal events (this paper\'s within-COVID design), rho=0.80: '
          f'resolvable delta-AUC = {d25:.3f} (paper states "roughly 0.10")', flush=True)
    # ~10^3 patients for 0.03-0.05 AUC, same event rate
    for n_fatal_target in (100, 104):
        n_total = round(n_fatal_target / 0.104)
        d = resolvable_delta_auc(n_fatal_target, n_total - n_fatal_target, rho=0.8)
        print(f'  {n_fatal_target} fatal events (~{n_total} patients at this paper\'s event rate), '
              f'rho=0.80: resolvable delta-AUC = {d:.3f} (paper states "0.03-0.05 requires ~10^3 patients")',
              flush=True)
