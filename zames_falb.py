"""
Zames-Falb multipliers: a tighter rigorous bound by using more of what we
actually know about freeplay.

Popov only used "sector [0,k], memoryless, time-invariant". But the freeplay
restoring moment is also:
    ODD                 M(-beta) = -M(beta)   (symmetric gap, no preload)
    MONOTONE            dM/dbeta >= 0 everywhere
    SLOPE-RESTRICTED    dM/dbeta in [0, K_beta]
Zames-Falb multipliers are exactly the class that exploits monotonicity, so
they should - in principle - certify a higher speed than Popov.

CRITERION. For phi slope-restricted in [0,k], the loop is stable if there
exists a multiplier M(s) = 1 - Hhat(s), with h in L1 and ||h||_1 <= 1, s.t.

    Re[ M(jw) (1/k + G(jw)) ] > 0   for all w,

and h(t) >= 0 for all t. For an ODD nonlinearity the sign restriction on h is
dropped (Turner, "Zames-Falb multipliers: don't panic", arXiv 2106.15913,
Theorems 4 and 5). Symmetric freeplay is odd. Freeplay about a PRELOADED
operating point is not (see preload.py), so odd=False is the only valid class
there: every basis coefficient must then be non-negative, since each basis
term has a non-negative impulse response.

PARAMETRIZATION. Basis  psi_i(s) = (rho/(s+rho))^i   (causal)
                        psi_i(s) = (rho/(-s+rho))^i  (anticausal)
each of which has impulse response with unit L1 norm, so ||h||_1 <= sum|c_i|
and the constraint sum|c_i| <= 1 is sufficient. The certification condition is
AFFINE in c, so feasibility is a linear program on a frequency grid.
"""
import numpy as np
from scipy.optimize import linprog, brentq
from conner3dof import ConnerSection
from absolute_stability import loop, G


def basis(w, rho, N):
    """Causal and anticausal first-order-lag powers evaluated on jw."""
    s = 1j*w
    out = []
    for i in range(1, N + 1):
        out.append((rho/(s + rho))**i)
        out.append((rho/(-s + rho))**i)
    return np.array(out)                       # (2N, nw)


def zf_margin(sec, U, w, rhos=(6., 20., 60., 200.), N=4, sd=True, odd=True):
    """Max over multipliers of the worst-case certification margin.
    > 0 means Zames-Falb certifies absolute stability at this speed."""
    A, B, C, k, _ = loop(sec, U, sd)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    F = 1/k + G(A, B, C, w)
    best = -np.inf
    for rho in rhos:
        P = basis(w, rho, N)                   # (m, nw)
        m = P.shape[0]
        # maximize t  s.t.  Re[F] - sum_i c_i Re[psi_i F] >= t ,  sum|c_i| <= 1
        # variables: [t, c+ (m), c- (m)],  c = c+ - c-
        R = (P*F).real                         # (m, nw)
        nw = len(w)
        Aub = np.zeros((nw + 1, 1 + 2*m))
        Aub[:nw, 0] = 1.0
        Aub[:nw, 1:1+m] = R.T
        Aub[:nw, 1+m:] = -R.T
        bub = np.concatenate([F.real, [1.0]])
        Aub[nw, 1:] = 1.0
        cobj = np.zeros(1 + 2*m); cobj[0] = -1.0
        bnds = [(None, None)] + [(0, None)]*m + [(0, None) if odd else (0, 0)]*m
        r = linprog(cobj, A_ub=Aub, b_ub=bub, bounds=bnds, method="highs")
        if r.status == 0:
            best = max(best, float(r.x[0]))
    return best


if __name__ == "__main__":
    sec = ConnerSection(); UF = 23.970
    w = np.geomspace(1e-3, 1e3, 1200)
    print("Zames-Falb margin vs speed (positive = certified):")
    for r in (0.10, 0.14, 0.16, 0.18, 0.20, 0.24, 0.27):
        print(f"   U/U_f={r:.2f}   margin = {zf_margin(sec, r*UF, w):+.5f}")
    lo, hi = 0.05*UF, 0.284*UF
    if zf_margin(sec, hi, w) > 0:
        Uz = hi
    else:
        Uz = brentq(lambda U: zf_margin(sec, U, w), lo, hi, xtol=1e-3)
    print(f"\n  ZAMES-FALB certified speed: {Uz:.3f} m/s = {Uz/UF:.4f} U_f")
    print(f"  Popov (previous rigorous bound):            0.16115 U_f")
    print(f"  earliest observed LCO (2.12 deg gap):       0.18    U_f")


# ---------------------------------------------------------------------------
def combined_margin(sec, U, w, rhos=(6., 20., 60., 200.), N=4,
                    etas=None, sd=True, odd=True):
    """Popov PLUS Zames-Falb,  M(s) = 1 + eta s - Hhat(s).

    Zames-Falb ALONE fails on this plant, and the reason is structural: with
    ||h||_1 <= 1 the multiplier obeys |M(jw)| <= 2, so it can only supply
    bounded phase advance. Popov's eta s is unbounded, and that is what this
    loop needs. The sum is a conic combination of two separately valid
    conditions, which is how the literature adds a Popov term to a Zames-Falb
    search (Carrasco, Maya-Gonzalez, Lanzon & Heath, SCL 2014). eta = 0 gives
    Zames-Falb, h = 0 gives Popov, so it can never do worse than either.

    The August version used the PRODUCT (1 + eta s)(1 - Hhat(s)). No theorem
    I could find covers that product, and it is falsified on this plant: it
    certified 0.2371 U_f (0.2356 U_f with h >= 0) while the plant sustains a
    limit cycle at 0.18 U_f. The sum certifies 0.1613 U_f, which is
    consistent with every observed onset.

    For fixed eta the condition is affine in c, so it stays an LP; eta is
    gridded outside. Search eta over many decades - a narrow range silently
    returns the endpoint, which is the bug that bit the plain Popov search.
    """
    A, B, C, k, _ = loop(sec, U, sd)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    F = 1/k + G(A, B, C, w)
    etas = etas if etas is not None else np.concatenate(
        [[0.0], np.geomspace(1e-5, 1e6, 45)])
    nw = len(w)
    best = -np.inf
    for eta in etas:
        Fe = (1 + 1j*w*eta)*F
        for rho in rhos:
            P = basis(w, rho, N)
            m = P.shape[0]
            R = (P*F).real
            Aub = np.zeros((nw + 1, 1 + 2*m))
            Aub[:nw, 0] = 1.0
            Aub[:nw, 1:1+m] = R.T
            Aub[:nw, 1+m:] = -R.T
            bub = np.concatenate([Fe.real, [1.0]])
            Aub[nw, 1:] = 1.0
            cobj = np.zeros(1 + 2*m); cobj[0] = -1.0
            r = linprog(cobj, A_ub=Aub, b_ub=bub,
                        bounds=[(None, None)] + [(0, None)]*m
                        + [(0, None) if odd else (0, 0)]*m,
                        method="highs")
            if r.status == 0:
                best = max(best, float(r.x[0]))
    return best
