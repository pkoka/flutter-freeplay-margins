"""
Freeplay under preload: the biased equilibrium, what absolute stability
actually certifies about it, and a describing-function prediction of the limit
cycle checked against time-domain simulation.

WHY. The memo's Section 6 certified stability of the ORIGIN of a symmetric
deadband. In practice a control surface is under load, so it is pushed to one
side of the gap and only sees the freeplay when it swings back through it
(Ting, Sep 2026). Under preload the origin is not an equilibrium, and a
certificate that does not state its operating point says nothing.

PRELOAD. A constant forcing f0 in zdot, through the existing integrator's
u_func/B input, so the Henon switching is unchanged. Two physical sources:
  hinge   a bench preload M0 on the flap hinge (spring, cable, gravity),
          f0 = B_phi*M0, written as p = M0/(K_beta*delta): how far M0 would
          push the flap past the stop, in gap widths, in still air.
  alpha0  a mounting angle of attack alpha0. The steady aero load on that
          rigid offset enters every aerodynamic term alpha does (including the
          Jones lag states) but not the pitch spring. It scales with U^2.

EQUILIBRIUM. z* solves 0 = A_r z + a_r + f0 in a region r that contains
beta*. Each region gives one candidate; exactly one is consistent here.

SHIFTED LOOP. With z = z* + x,
    xdot = A_free x - B_phi * phi_s(C_beta x),
    phi_s(s) = phi(beta* + s) - phi(beta*).
phi is monotone with slope in [0, K_beta], so phi_s is too, and phi_s(0) = 0,
which puts phi_s in the sector [0, K_beta]. The linear part G is unchanged.
So the Popov certificate on [0, K_beta] is the same NUMBER for every preload,
but what it certifies is global asymptotic stability of z*(U, preload), not of
the origin. phi_s is not odd unless beta* = 0, so only Zames-Falb multipliers
with h >= 0 apply (zames_falb.py, odd=False).

THE SECTOR NARROWS ON A STOP. If the flap rests on a stop, |beta*| > delta,
the smallest ratio phi_s(s)/s is reached when the swing just reaches the far
edge of the gap:
    k1 = K_beta (|beta*| - delta)/(|beta*| + delta)  > 0.
So phi_s lies in [k1, K_beta], and Popov after the loop transformation
(absolute_stability.popov_margin(k1=...)) certifies a preload-dependent speed
that rises steeply with preload. Inside the gap k1 = 0 and nothing is gained.

DESCRIBING FUNCTION. Dual-input: beta = B0 + A sin(wt) through the freeplay
gives a mean N0(B0, A) and a first-harmonic gain N1(B0, A). Harmonic balance
on the mean and on the first harmonic:
    B0 + G(0) N0 = beta_f,     1 + N1 G(jw) = 0,
with beta_f = -C_beta A_free^-1 f0 the free-region static flap response to the
preload. N1 is real, so w is a phase crossover of G, and N1 = gamma*K_beta
with gamma = -1/(K_beta Re G(jw)) in (0, 1). At zero preload B0 = 0 solves the
mean equation and this is the memo's gamma family exactly.

TIME DOMAIN. Every point starts from z* plus 16 random full-state
perturbations (displacements 1 to 12 gap widths, with velocities), because one
initial condition reports which attractor that one trajectory fell into, not
what the plant does. Flap-position kicks alone are not enough: the
subcritical limit cycle at 0.18 U_f is only reached from starts that carry
velocity, and a position-only sweep reports its onset at 0.30 U_f.
"""
import os
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")   # one BLAS thread per worker
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import inv
from scipy.optimize import brentq
from concurrent.futures import ProcessPoolExecutor
from conner3dof import ConnerSection, state_space, aero_matrices, ss_flutter
from freeplay import regions, integrate
from absolute_stability import loop, G, popov_margin, certified_speed
from zames_falb import zf_margin, combined_margin

DELTA = np.deg2rad(2.12)
RATIOS = np.round(np.arange(0.10, 0.701, 0.02), 3)
N_IC = 16                                          # random starts per point
CASES = (("hinge", 0.0), ("hinge", 0.25), ("hinge", 0.5), ("hinge", 1.0),
         ("hinge", 2.0), ("alpha0", 2.0), ("alpha0", 4.0), ("alpha0", 6.0))
T_SIM, T_REC = 40.0, 8.0
TH = (np.arange(4096) + 0.5)*2*np.pi/4096
SIN = np.sin(TH)


def hinge_direction(sec, U):
    """How a flap hinge moment enters zdot (the B_phi of absolute_stability)."""
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    b = np.zeros(8); b[3:6] = inv(sec.Ms + A_non) @ np.array([0.0, 1.0, 0.0])
    return b


def preload_forcing(sec, U, kind, level, delta=DELTA):
    """f0 for a bench hinge preload (level = p) or a mounting angle (level in
    degrees). Zero level gives zero forcing for either kind."""
    if kind == "hinge":
        return hinge_direction(sec, U)*sec.Cb*delta*level
    if kind == "alpha0":
        A = state_space(sec, U, True, "free")
        A_non, _, _ = aero_matrices(sec, U, 0.0)
        col = A[:, 0].copy()
        col[3:6] += inv(sec.Ms + A_non) @ sec.Ks[:, 0]   # remove the pitch spring
        return col*np.deg2rad(level)
    raise ValueError(kind)


def freeplay_moment(beta, delta, k):
    return k*(beta - np.clip(beta, -delta, delta))


def field(sec, U, delta, f0, z):
    """The true piecewise vector field, for residual checks."""
    As, Af, ap, an = regions(sec, U, delta)
    b = z[1]
    if b > delta:
        return As @ z + ap + f0
    if b < -delta:
        return As @ z + an + f0
    return Af @ z + f0


def equilibria(sec, U, delta, f0):
    """Every region-consistent equilibrium, as (region, z*, residual)."""
    As, Af, ap, an = regions(sec, U, delta)
    out = []
    for reg, A, a, inside in ((2, Af, 0.0, lambda b: abs(b) <= delta),
                              (3, As, ap, lambda b: b > delta),
                              (1, As, an, lambda b: b < -delta)):
        z = -np.linalg.solve(A, a + f0)
        if inside(z[1]):
            out.append((reg, z, float(np.abs(field(sec, U, delta, f0, z)).max())))
    return out


def equilibrium(sec, U, delta, f0):
    eq = equilibria(sec, U, delta, f0)
    assert len(eq) == 1, f"{len(eq)} equilibria at U={U:.3f}"
    return eq[0]


def sector_lower(beta_star, delta, k):
    """Exact lower sector bound of phi_s, derived in the module docstring."""
    b = abs(beta_star)
    return k*(b - delta)/(b + delta) if b > delta else 0.0


def certified_preloaded(sec, UF, kind, lev, delta=DELTA, w=None):
    """Highest speed below which Popov on [k1(U), K_beta] certifies z*(U)."""
    w = np.geomspace(1e-4, 1e4, 4000) if w is None else w

    def m(U):
        f0 = preload_forcing(sec, U, kind, lev, delta)
        _, zs, _ = equilibrium(sec, U, delta, f0)
        return popov_margin(sec, U, w, k1=sector_lower(zs[1], delta, sec.Cb))
    Us = np.linspace(0.05, 0.99, 95)*UF
    for a, b in zip(Us[:-1], Us[1:]):
        if m(b) <= 0:
            return brentq(m, a, b, xtol=1e-4)
    return Us[-1]


def shifted_sector(beta_star, delta, k, span=50.0, n=200001):
    """min and max of phi_s(s)/s over s in [-span, span]*delta, s != 0."""
    s = np.linspace(-span, span, n)*delta
    s = s[np.abs(s) > 1e-9*delta]
    r = (freeplay_moment(beta_star + s, delta, k)
         - freeplay_moment(beta_star, delta, k))/s
    odd = np.abs(freeplay_moment(beta_star + s, delta, k)
                 - freeplay_moment(beta_star, delta, k)
                 + freeplay_moment(beta_star - s, delta, k)
                 - freeplay_moment(beta_star, delta, k)).max()
    return r.min()/k, r.max()/k, odd/(k*delta)


# --------------------------------------------------------------------------
# dual-input describing function and harmonic balance
# --------------------------------------------------------------------------
def didf(B0, A, delta, k):
    """(N0, N1): mean output and first-harmonic gain of the freeplay for the
    input B0 + A sin(theta). Midpoint quadrature on 4096 points; the moment is
    piecewise linear, so this is accurate to about 1e-7 relative."""
    phi = freeplay_moment(B0 + A*SIN, delta, k)
    return phi.mean(), 2*np.mean(phi*SIN)/A


def phase_crossovers(A, B, C, w=None):
    """Frequencies where Im G = 0 and Re G < 0, refined by root finding."""
    w = np.geomspace(1.0, 400.0, 6000) if w is None else w
    g = G(A, B, C, w)
    out = []
    for i in np.nonzero(np.diff(np.sign(g.imag)) != 0)[0]:
        f = lambda wi: G(A, B, C, [wi])[0].imag
        wc = brentq(f, w[i], w[i + 1], xtol=1e-10)
        gc = G(A, B, C, [wc])[0]
        if gc.real < 0:
            out.append((wc, gc.real))
    return out


def bias_for_gain(A, gam, delta, k):
    """The B0 >= 0 with N1(B0, A) = gam*k, or None. N1 is even in B0 and
    non-decreasing in |B0| (more of the swing on the stiff side), reaching k
    once B0 >= A + delta, so the root is unique when N1(0, A) <= gam*k."""
    f = lambda B0: didf(B0, A, delta, k)[1]/k - gam
    f0 = f(0.0)
    if f0 > 1e-12:
        return None
    if f0 >= 0.0:                     # at A_sym to rounding: symmetric cycle
        return 0.0
    return brentq(f, 0.0, A + delta, xtol=1e-14)


def harmonic_balance(sec, U, delta, f0, n=240):
    """All (f_Hz, gamma, B0, A) satisfying the dual-input harmonic balance.

    For each phase crossover the gain fixes gamma. Along the amplitude axis,
    bias_for_gain gives B0 = +/- b(A), and the mean equation's residual is
    scanned for sign changes, then refined. A ranges up to the symmetric
    amplitude A_sym(gamma), where b = 0."""
    Af, B, C, k, _ = loop(sec, U)
    G0 = (C @ np.linalg.solve(-Af, B)).item()
    beta_f = (C @ np.linalg.solve(-Af, f0)).item()
    sols = []
    for wc, reG in phase_crossovers(Af, B, C):
        gam = -1.0/(k*reG)
        if not 0.0 < gam < 1.0:
            continue
        A_sym = brentq(lambda A: didf(0.0, A, delta, k)[1]/k - gam,
                       delta*(1 + 1e-9), 1e5*delta, xtol=1e-14)
        As = np.concatenate([np.geomspace(1e-3*delta, A_sym, n)[:-1], [A_sym]])
        for sgn in (1.0, -1.0):
            def R(A):
                if A == A_sym:                # b = 0 there by definition, and
                    return -beta_f/delta      # N0(0, A) = 0 because phi is odd
                b = bias_for_gain(A, gam, delta, k)
                if b is None:
                    return np.nan
                b *= sgn
                return (b + G0*didf(b, A, delta, k)[0] - beta_f)/delta
            Rv = np.array([R(A) for A in As])
            hits = [brentq(R, As[i], As[i + 1], xtol=1e-14)
                    for i in np.nonzero(np.sign(Rv[:-1])*np.sign(Rv[1:]) < 0)[0]]
            if abs(Rv[-1]) < 1e-10:
                hits.append(A_sym)
            for A in hits:
                B0 = 0.0 if A == A_sym else sgn*bias_for_gain(A, gam, delta, k)
                cand = (wc/(2*np.pi), gam, B0, A)
                if not any(abs(c[2] - B0) < 1e-6*delta and abs(c[3] - A) < 1e-6*delta
                           and abs(c[0] - cand[0]) < 1e-6 for c in sols):
                    sols.append(cand)
    return sols


# --------------------------------------------------------------------------
# time domain
# --------------------------------------------------------------------------
def random_start(zs, delta, b, seed):
    """z* plus a random displacement-and-velocity perturbation, scaled by a
    random 1 to 12 gap widths (plunge scaled by the semichord)."""
    rng = np.random.default_rng(seed)
    sc = rng.uniform(1.0, 12.0)*delta*np.array([1.0, 1.0, b])
    return zs + np.concatenate([rng.standard_normal(3)*sc,
                                rng.standard_normal(3)*sc*30.0, [0.0, 0.0]])


def _simulate(args):
    """One trajectory from a random start about z*; report the recorded window."""
    sec, U, delta, f0, zs, kick = args
    z0 = random_start(zs, delta, sec.b, kick)
    t, z, _ = integrate(sec, U, delta, T=T_SIM, z0=z0, u_func=lambda s: 1.0,
                        B=f0[:, None], record_from=T_SIM - T_REC)
    b = z[:, 1]
    n2 = len(b)//2
    amp = lambda x: 0.5*(x.max() - x.min())
    reg = np.where(b > delta, 3, np.where(b < -delta, 1, 2))
    return dict(kick=kick, mean=b.mean()/delta, amp=amp(b)/delta,
                amp_first=amp(b[:n2])/delta, amp_last=amp(b[n2:])/delta,
                switches=int(np.count_nonzero(np.diff(reg))),
                alpha_amp=amp(z[:, 0])/delta, diverged=bool(np.abs(z[:, :3]).max() > 0.5))


def sustained(r, tol=1e-3):
    """A limit cycle: finite, not decaying across the recorded window."""
    return (not r["diverged"]) and r["amp_last"] > tol and \
        r["amp_last"] > 0.95*r["amp_first"]


def time_domain(pool, sec, UF, delta=DELTA, cases=CASES, ratios=RATIOS,
                kicks=tuple(range(N_IC))):
    jobs, keys = [], []
    for kind, lev in cases:
        for r in ratios:
            U = r*UF
            f0 = preload_forcing(sec, U, kind, lev, delta)
            _, zs, _ = equilibrium(sec, U, delta, f0)
            for kk in kicks:
                jobs.append((sec, U, delta, f0, zs, kk)); keys.append((kind, lev, r))
    res = {}
    for key, out in zip(keys, pool.map(_simulate, jobs, chunksize=4)):
        res.setdefault(key, []).append(out)
    return res


# --------------------------------------------------------------------------
def main():
    sec = ConnerSection()
    UF, _ = ss_flutter(sec, np.linspace(1.0, 40.0, 3000), True)
    d = DELTA; k = sec.Cb
    print(f"U_f = {UF:.3f} m/s, gap delta = 2.12 deg, K_beta = {k:.5f} N m/rad\n")

    # ------------------------------------------------------------- 1
    print("=" * 78); print("1. OPERATING POINT z*(U, preload)"); print("=" * 78)
    print(f"{'case':>12} " + " ".join(f"{r:>7.2f}" for r in (0.1, 0.2, 0.3, 0.5, 0.7))
          + "    <- U/U_f;  entries are beta*/delta [region]")
    worst = 0.0
    for kind, lev in CASES:
        row = []
        for r in (0.1, 0.2, 0.3, 0.5, 0.7):
            f0 = preload_forcing(sec, r*UF, kind, lev)
            reg, zs, res = equilibrium(sec, r*UF, d, f0)
            worst = max(worst, res)
            row.append(f"{zs[1]/d:+6.2f}{'fsS'[[2, 1, 3].index(reg)] if reg != 2 else ' '}")
        print(f"{kind+' '+f'{lev:g}':>12} " + " ".join(f"{c:>7}" for c in row))
    print(f"  exactly one region-consistent equilibrium at every point; "
          f"max residual {worst:.1e}")
    print("  (s/S: flap resting on the negative/positive stop; blank: inside the gap)")
    print("  alpha0 floats the flap at a fixed fraction of the gap at every speed, "
          "because\n  the aero hinge moment and the aero hinge stiffness both "
          "scale with U^2.")

    # ------------------------------------------------------------- 2
    print("\n" + "=" * 78); print("2. WHAT ABSOLUTE STABILITY CERTIFIES ABOUT z*"); print("=" * 78)
    for kind, lev in CASES:
        f0 = preload_forcing(sec, 0.15*UF, kind, lev)
        _, zs, _ = equilibrium(sec, 0.15*UF, d, f0)
        lo, hi, odd = shifted_sector(zs[1], d, k)
        print(f"  {kind+' '+f'{lev:g}':>12}: beta*/delta {zs[1]/d:+6.2f}  "
              f"phi_s/s in [{lo:.4f}, {hi:.4f}] K_beta,  "
              f"oddness defect {odd:.3f} K_beta delta")
    w = np.geomspace(1e-4, 1e4, 4000)
    Up = certified_speed(sec, popov_margin, w, lo=0.03*UF, hi=0.25*UF)
    print(f"\n  Popov certified speed, any constant preload: {Up/UF:.4f} U_f")
    print("  The number is unchanged because G and the sector are unchanged. What")
    print("  changed is the statement: below it, z*(U, preload) is globally")
    print("  asymptotically stable. It is NOT the origin unless the preload is zero.")
    print("\n  On a stop the sector narrows to [k1, K_beta]; Popov after the loop")
    print("  transformation then certifies a preload-dependent speed:")
    cert = {}
    for kind, lev in CASES:
        cert[(kind, lev)] = certified_preloaded(sec, UF, kind, lev)/UF
        f0 = preload_forcing(sec, cert[(kind, lev)]*UF, kind, lev)
        _, zs, _ = equilibrium(sec, cert[(kind, lev)]*UF, d, f0)
        print(f"    {kind+' '+f'{lev:g}':>12}: k1 = {sector_lower(zs[1], d, k)/k:.4f} K_beta"
              f"   certified {cert[(kind, lev)]:.4f} U_f")
    wl = np.geomspace(1e-3, 1e3, 1500)
    out = []
    for name, fn in (("Zames-Falb alone, odd class", lambda U: zf_margin(sec, U, wl)),
                     ("Zames-Falb alone, h >= 0", lambda U: zf_margin(sec, U, wl, odd=False)),
                     ("Popov + Zames-Falb, odd class", lambda U: combined_margin(sec, U, wl)),
                     ("Popov + Zames-Falb, h >= 0", lambda U: combined_margin(sec, U, wl, odd=False))):
        if fn(0.10*UF) <= 0:
            out.append((name, None)); continue
        hi = 0.284*UF
        Uc = hi if fn(hi) > 0 else brentq(fn, 0.10*UF, hi, xtol=1e-3)
        out.append((name, Uc/UF))
    for name, v in out:
        print(f"  {name:32s} " + (f"{v:.4f} U_f" if v else "nothing certified at 0.10 U_f"))
    print("  Only the h >= 0 rows apply under preload. The odd rows apply only at")
    print("  zero preload, where the freeplay about z* = 0 is odd.")

    # ------------------------------------------------------------- 3, 4
    with ProcessPoolExecutor(max_workers=12) as pool:
        TD = time_domain(pool, sec, UF)

    print("\n" + "=" * 78); print("3. CHECK: BELOW THE CERTIFIED SPEED EVERY START RETURNS TO z*"); print("=" * 78)
    print(f"  a valid certificate cannot be beaten by any simulated trajectory ({N_IC} random")
    print("  starts per speed, displacements up to 12 gap widths, with velocity).")
    for key in CASES:
        below = [x for (kd, lv, r), rs in TD.items() for x in rs
                 if (kd, lv) == key and r < cert[key]]
        worst = max((x["amp_last"] for x in below), default=0.0)
        print(f"  {key[0]+' '+f'{key[1]:g}':>12}: {len(below):4d} trajectories below "
              f"{cert[key]:.4f} U_f, largest motion left after {T_SIM:.0f} s "
              f"{worst:.1e} delta " + ("ok" if worst < 1e-3 else "VIOLATION"))

    print("\n" + "=" * 78); print("4. LIMIT CYCLES: DESCRIBING FUNCTION AGAINST SIMULATION"); print("=" * 78)
    print("  amplitudes are half peak-to-peak flap motion about its mean, in gap widths.")
    print(f"  sim: the largest sustained cycle over {N_IC} random starts, and how many of")
    print("  them end on a cycle. DF: every harmonic-balance solution (amplitude/mean).")
    onsets = {}
    for kind, lev in CASES:
        print(f"\n  {kind} {lev:g}")
        print(f"  {'U/Uf':>6} {'beta*/d':>8} {'sim amp':>8} {'sim mean':>9} {'cyc':>4}"
              f" {'f sim':>6}   DF solutions (A/d @ B0/d, f Hz)")
        on_sim = on_df = None
        for r in RATIOS:
            U = r*UF
            f0 = preload_forcing(sec, U, kind, lev)
            _, zs, _ = equilibrium(sec, U, d, f0)
            rs = TD[(kind, lev, r)]
            cyc = [x for x in rs if sustained(x)]
            div = any(x["diverged"] for x in rs)
            hb = harmonic_balance(sec, U, d, f0)
            if cyc and on_sim is None:
                on_sim = r
            if hb and on_df is None:
                on_df = r
            if cyc:
                best = max(cyc, key=lambda x: x["amp_last"])
                simtxt = f"{best['amp_last']:>8.3f} {best['mean']:>+9.3f} {len(cyc):>4d}"
            else:
                simtxt = f"{'diverge' if div else '-':>8} {'':>9} {0:>4d}"
            hbtxt = "  ".join(f"{A/d:.2f}@{B0/d:+.2f},{f:.1f}" for f, g, B0, A in
                              sorted(hb, key=lambda s: s[3]))
            if r <= 0.40 or abs(r*100 % 10) < 1e-6:
                print(f"  {r:>6.2f} {zs[1]/d:>+8.2f} {simtxt}        {hbtxt or '-'}")
        onsets[(kind, lev)] = (on_sim, on_df)

    print("\n" + "=" * 78); print("5. SUMMARY: ONSET OF SUSTAINED MOTION"); print("=" * 78)
    print(f"  {'case':>12} {'Popov':>8} {'DF onset':>9} {'sim onset':>10}")
    for (kind, lev), (on_sim, on_df) in onsets.items():
        print(f"  {kind+' '+f'{lev:g}':>12} {cert[(kind, lev)]:>8.4f} "
              f"{(f'{on_df:.2f}' if on_df else '-'):>9} {(f'{on_sim:.2f}' if on_sim else '-'):>10}")
    print(f"  grid step {RATIOS[1]-RATIOS[0]:.2f} U_f; onsets are the first grid speed with")
    print("  a solution (DF) or a sustained cycle from some random start (sim).")

    print("\n" + "=" * 78); print("6. IF THE PRELOAD IS MISESTIMATED BY 25 PERCENT"); print("=" * 78)
    print("  DF amplitude predicted with the preload level scaled, against the true")
    print("  plant's simulated cycle. This is the bench-versus-sweep question: the")
    print("  DF needs the preload as an input, so its error carries straight through.")
    for kind, lev in (("hinge", 1.0), ("alpha0", 4.0)):
        for r in (0.30, 0.50):
            U = r*UF
            rs = [x for x in TD[(kind, lev, r)] if sustained(x)]
            sim = max(x["amp_last"] for x in rs) if rs else np.nan
            row = []
            for fac in (0.75, 1.0, 1.25):
                hb = harmonic_balance(sec, U, d, preload_forcing(sec, U, kind, lev*fac))
                A = [a/d for _, _, _, a in hb]
                row.append(min(A, key=lambda a: abs(a - sim)) if A and np.isfinite(sim) else np.nan)
            print(f"  {kind} {lev:g} at {r:.2f} U_f: sim {sim:.3f}   DF x0.75 {row[0]:.3f}"
                  f"   x1.00 {row[1]:.3f}   x1.25 {row[2]:.3f}   (closest branch)")
    np.savez("preload_results.npz", ratios=RATIOS, Up=Up/UF,
             cert=np.array([cert[c] for c in CASES]),
             cases=np.array([f"{k}:{l}" for k, l in CASES]),
             sim_amp=np.array([[max([x["amp_last"] for x in TD[(k, l, r)] if sustained(x)],
                                    default=0.0) for r in RATIOS] for k, l in CASES]),
             sim_mean=np.array([[max([x for x in TD[(k, l, r)] if sustained(x)],
                                     key=lambda x: x["amp_last"], default={"mean": np.nan})["mean"]
                                 for r in RATIOS] for k, l in CASES]))


if __name__ == "__main__":
    main()
