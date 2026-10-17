"""
Closed-loop back-calculation of the open-loop plant on the Conner rig.

WHY. Near and past flutter the open-loop plant cannot be identified open loop:
the input stops correlating with the output as the structure goes unstable
(Ting's limitation #3), and stage 3 showed that the dwell an open-loop test
needs grows as 1/(zeta w) on the way there. The field's answer (Ting, Sep
2026): run the flutter suppression system in the loop, sweep, identify the
closed loop, and back-calculate the open-loop plant. This does that in
simulation, where the truth is known, and asks the question Ting was asked:
is the controller treated as exactly known, or divided out from measurement?

SETUP
  truth   perturbed plant (stage2a_harness.perturbed_section, 4% scatter),
          U_f = 25.83 m/s = 1.0775 U_f(nominal); flap hinge-moment actuator
          with a 30 Hz first-order lag, held between samples; unmeasured
          turbulence force at the shaker station; RVDT noise; 16-bit ADC.
  design  nominal Conner model, no actuator lag. The controller never sees the
          truth, as on a real AFS.
  AFS     discrete LQG at 800 Hz, designed at each test speed on the nominal
          ZOH model with the designer's actuator estimate (40 Hz; truth 30).
          Strictly proper, so the computation delay is inherent.
          The stage-2 H-infinity central controllers carry poles near 1e5
          rad/s and cannot run at 800 Hz, which is why this is not them.
          The AS-BUILT controller differs from the DOCUMENTED one by a 5% gain
          error and a 100 Hz output smoothing filter.
  test    periodic multisine r added to the controller output (u = K y + r);
          u and y logged at 800 Hz. Integer periods, transient periods
          dropped, spectra at the excited lines only (stage 2a estimator).

ESTIMATORS, per excited line
  direct      H1 from u and y alone: sum Y U^* / sum |U|^2. r is not used.
  ind exact   T_yr measured; T_ur = 1 + K T_yr with the AS-BUILT controller
  ind doc     the same with the DOCUMENTED controller
  joint I/O   T_yr / T_ur, both measured against r. No controller model.

The open-loop plant being estimated is command u to measured y: the actuator
is part of it, the controller is not. Truth is the exact sampled-data FRF.
"""
import os
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")   # one BLAS thread per worker
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import expm, solve_discrete_are, inv
from concurrent.futures import ProcessPoolExecutor
from conner3dof import ConnerSection, state_space, ss_flutter
from phase1_conner import input_matrix, output_matrices
from freeplay import regions, step
from stage2a_harness import perturbed_section
from preload import hinge_direction, preload_forcing, equilibrium, DELTA

FS = 800.0; TS = 1/FS
DT_SIM = 1.25e-4; NSUB = 10                   # integrator steps per sample
F_ACT, F_ACT_DESIGN = 30.0, 40.0             # actuator bandwidth: truth, designer's
F_SMOOTH, GAIN_ERR = 100.0, 0.05              # what the as-built controller adds
RHO = 0.01                                    # LQR control weight
P, NPER, DROP = 8.0, 6, 2
BAND = (0.5, 30.0)
YSCALE = np.array([0.02, 0.05, 0.002])        # rad, rad, m (as stage2_hinf)
M_REF = 0.1                                   # N m per unit span
SIG_V = np.array([5e-5, 1e-4, 5e-6])          # RVDT noise: rad, rad, m
SIG_W = 0.5                                   # turbulence force, N per unit span
FSR = np.array([0.17, 0.17, 0.02])            # ADC full scale: +/-10 deg, +/-2 cm
R_AMP = 0.04                                  # multisine amplitude, N m per unit span
RATIOS = (0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20)
SEEDS = (11, 23, 37, 51, 67, 83)
NAMES = ("direct", "ind exact", "ind doc", "joint I/O")
OUT = ("alpha", "beta", "h")


# --------------------------------------------------------------------------
# sampled plant, controller
# --------------------------------------------------------------------------
def zoh(A, B):
    n, m = B.shape
    M = np.zeros((n + m, n + m)); M[:n, :n] = A; M[:n, n:] = B
    E = expm(M*TS)
    return E[:n, :n], E[:n, n:]


def sampled_plant(sec, U, f_act=F_ACT, region="stiff"):
    """(Phi, Gu, Gw, C) from command u_k and turbulence w_k to y_k. With an
    actuator (f_act in Hz, None for none) the state carries the held hinge
    moment m_{k-1}, and m_k = a m_{k-1} + (1 - a) u_k acts over [t_k, t_k+1)."""
    A = state_space(sec, U, True, region)
    Phi, Gm = zoh(A, hinge_direction(sec, U)[:, None])
    _, Gw = zoh(A, input_matrix(sec, U))
    C, _ = output_matrices()
    if f_act is None:
        return Phi, Gm, Gw, C
    a = np.exp(-2*np.pi*f_act*TS)
    Pa = np.zeros((9, 9)); Pa[:8, :8] = Phi; Pa[:8, 8:] = a*Gm; Pa[8, 8] = a
    return (Pa, np.vstack([(1 - a)*Gm, [[1 - a]]]), np.vstack([Gw, [[0.0]]]),
            np.hstack([C, np.zeros((3, 1))]))


def lqg(sec, U):
    """Documented controller: LQR on normalised outputs plus a steady-state
    Kalman predictor, both on the nominal ZOH model. (Ac, Bc, Cc) with
    xc' = Ac xc + Bc y, u_fb = Cc xc."""
    Phi, Gu, Gw, C = sampled_plant(sec, U, f_act=F_ACT_DESIGN)
    n = Phi.shape[0]
    Cw = C/YSCALE[:, None]
    R = np.array([[RHO/M_REF**2]])
    X = solve_discrete_are(Phi, Gu, Cw.T @ Cw + 1e-12*np.eye(n), R)
    F = np.linalg.solve(R + Gu.T @ X @ Gu, Gu.T @ X @ Phi)
    W = (Gw @ Gw.T)*SIG_W**2 + (Gu @ Gu.T)*M_REF**2 + 1e-14*np.eye(n)
    V = np.diag(SIG_V**2)
    Y = solve_discrete_are(Phi.T, C.T, W, V)
    L = Phi @ Y @ C.T @ inv(C @ Y @ C.T + V)
    return Phi - Gu @ F - L @ C, L, -F


def as_built(K):
    """The documented controller with a 5% gain error and a first-order
    smoothing filter on its output, as implemented on the rig."""
    Ac, Bc, Cc = K
    n = Ac.shape[0]
    b, g = np.exp(-2*np.pi*F_SMOOTH*TS), 1 + GAIN_ERR
    A2 = np.zeros((n + 1, n + 1)); A2[:n, :n] = Ac
    A2[n, :n] = (1 - b)*Cc[0]; A2[n, n] = b
    return A2, np.vstack([Bc, np.zeros((1, 3))]), g*np.hstack([(1 - b)*Cc, [[b]]])


def frf(A, B, C, f):
    """Discrete-time frequency response C (zI - A)^-1 B at z = exp(j 2 pi f TS)."""
    n = A.shape[0]
    return np.stack([C @ np.linalg.solve(np.exp(2j*np.pi*fi*TS)*np.eye(n) - A, B)
                     for fi in f], axis=-1)


def closed_loop_matrix(sec_t, U, K):
    Pa, Ga, _, Ca = sampled_plant(sec_t, U)
    Ac, Bc, Cc = K
    return np.block([[Pa, Ga @ Cc], [Bc @ Ca, Ac]])


def least_damped(eigs_z):
    s = np.log(np.asarray(eigs_z, complex))/TS
    osc = s[np.abs(s.imag) > 1e-3]
    return float(np.min(-osc.real/np.abs(osc)))


# --------------------------------------------------------------------------
# experiment
# --------------------------------------------------------------------------
def multisine(rng):
    n, df = int(P*FS), 1.0/P
    ks = np.arange(int(np.ceil(BAND[0]/df)), int(np.floor(BAND[1]/df)) + 1)
    ph = rng.uniform(0, 2*np.pi, len(ks))
    t = np.arange(n*NPER)/FS
    r = R_AMP*sum(np.cos(2*np.pi*k*df*t + p) for k, p in zip(ks, ph))/np.sqrt(len(ks))
    return r, ks


def quantise(y):
    q = 2*FSR/2**16
    return np.round(y/q[:, None])*q[:, None]


def run(args):
    """One test point. mode 'linear' marches the exact sampled plant; mode
    'freeplay' integrates the true switching plant with Henon between samples,
    with optional preload forcing f0 about its equilibrium. closed=False is an
    open-loop test (u = r). Returns spectra at the excited lines."""
    sec_t, U, K, seed, closed, mode, f0, nper, drop = args
    rng = np.random.default_rng(seed)
    r, ks = multisine(rng)
    r = np.tile(r[:int(P*FS)], nper)
    N = len(r)
    w = rng.standard_normal(N)*SIG_W
    v = rng.standard_normal((3, N))*SIG_V[:, None]
    Ac, Bc, Cc = K
    xc = np.zeros(Ac.shape[0])
    Y = np.zeros((3, N)); Uu = np.zeros(N)
    C, _ = output_matrices()
    a = np.exp(-2*np.pi*F_ACT*TS)
    if mode == "linear":
        Pa, Ga, Gwa, Ca = sampled_plant(sec_t, U)
        xi = np.zeros(9)
        for k in range(N):
            y = Ca @ xi + v[:, k]
            u = (Cc @ xc).item() + r[k] if closed else r[k]
            Y[:, k], Uu[k] = y, u
            xc = Ac @ xc + Bc @ y
            xi = Pa @ xi + Ga[:, 0]*u + Gwa[:, 0]*w[k]
        switches, bmax = 0, np.nan
    else:
        fields = regions(sec_t, U, DELTA)
        bm, bw = hinge_direction(sec_t, U), input_matrix(sec_t, U)[:, 0]
        _, zs, _ = equilibrium(sec_t, U, DELTA, f0)
        z, t, m, switches = zs.copy(), 0.0, 0.0, 0
        bdev = []
        for k in range(N):
            y = C @ (z - zs) + v[:, k]            # controller acts on deviation from trim
            u = (Cc @ xc).item() + r[k] if closed else r[k]
            Y[:, k], Uu[k] = y, u
            xc = Ac @ xc + Bc @ y
            m = a*m + (1 - a)*u
            force = bm*m + bw*w[k] + f0
            ext = lambda s, fv=force: fv
            for _ in range(NSUB):
                z, t, sw = step(z, t, DT_SIM, DELTA, fields, ext)
                switches += sw
            bdev.append(z[1])
            if not np.all(np.isfinite(z)) or np.abs(z[:3]).max() > 0.5:
                return dict(diverged=True)
        bmax = float(np.max(np.abs(np.array(bdev[drop*int(P*FS):]) - zs[1]))/DELTA)
    Y = quantise(Y)
    n = int(P*FS)
    Yp = np.array([np.fft.rfft(Y[:, p*n:(p + 1)*n], axis=1)[:, ks] for p in range(drop, nper)])
    Up = np.array([np.fft.rfft(Uu[p*n:(p + 1)*n])[ks] for p in range(drop, nper)])
    Rk = np.fft.rfft(r[:n])[ks]
    return dict(diverged=False, f=ks/P, Yp=Yp, Up=Up, R=Rk, switches=switches,
                bmax=bmax, snr=np.std(Y, axis=1)/SIG_V)


def estimates(res, Kas, Kdoc):
    """The four open-loop estimates (3, nf) from one run's spectra."""
    Yp, Up, Rk, f = res["Yp"], res["Up"], res["R"], res["f"]
    Yb, Ub = Yp.mean(0), Up.mean(0)
    Tyr = Yb/Rk
    G = {"direct": (Yp*np.conj(Up)[:, None, :]).sum(0)/(np.abs(Up)**2).sum(0),
         "joint I/O": Yb/Ub}
    for name, K in (("ind exact", Kas), ("ind doc", Kdoc)):
        Kf = frf(*K, f)[0]                                   # (3, nf)
        G[name] = Tyr/(1 + (Kf*Tyr).sum(0))
    return G


# --------------------------------------------------------------------------
# modal extraction
# --------------------------------------------------------------------------
def modal_fit(f, H, band=(1.5, 14.0), nd=6, iters=40):
    """Sanathanan-Koerner common-denominator rational fit in s over a band,
    channels normalised. Returns continuous poles (rad/s)."""
    sel = (f >= band[0]) & (f <= band[1])
    w = 2*np.pi*f[sel]; wm = w.max(); s = 1j*w/wm
    Hs = H[:, sel]/np.median(np.abs(H[:, sel]), axis=1, keepdims=True)
    m = len(w)
    S = np.vander(s, nd + 1, increasing=True)                # 1, s, ..., s^nd
    Dp = np.ones(m)
    for _ in range(iters):
        rows, rhs = [], []
        for i in range(3):
            blk = np.zeros((m, nd + 3*(nd + 1)), complex)
            blk[:, :nd] = -Hs[i][:, None]*S[:, 1:]           # -H * (d1 s + ... )
            blk[:, nd + i*(nd + 1): nd + (i + 1)*(nd + 1)] = S
            rows.append(blk/Dp[:, None]); rhs.append(Hs[i]/Dp)
        A = np.vstack(rows); b = np.concatenate(rhs)
        x, *_ = np.linalg.lstsq(np.vstack([A.real, A.imag]),
                                np.concatenate([b.real, b.imag]), rcond=None)
        d = x[:nd]
        Dn = 1 + S[:, 1:] @ d
        if np.allclose(Dn, Dp, rtol=1e-10, atol=0):
            break
        Dp = Dn
    return np.roots(np.concatenate([d[::-1], [1.0]]))*wm


def flutter_mode(poles, fband=(3.0, 9.0)):
    """(zeta, f_Hz) of the least-damped oscillatory pole in the band."""
    p = poles[(poles.imag > 1e-6) & (poles.imag/(2*np.pi) > fband[0])
              & (poles.imag/(2*np.pi) < fband[1])]
    if len(p) == 0:
        return np.nan, np.nan
    z = -p.real/np.abs(p)
    i = np.argmin(z)
    return float(z[i]), float(abs(p[i])/(2*np.pi))


def truth_mode(sec, U):
    ev = np.linalg.eigvals(state_space(sec, U, True, "stiff"))
    return flutter_mode(ev)


def zero_crossing(rs, zs):
    rs, zs = np.asarray(rs), np.asarray(zs)
    for i in range(len(rs) - 1):
        if np.isfinite(zs[i]) and np.isfinite(zs[i + 1]) and zs[i] > 0 >= zs[i + 1]:
            return rs[i] + zs[i]/(zs[i] - zs[i + 1])*(rs[i + 1] - rs[i])
    return np.nan


# --------------------------------------------------------------------------
def main():
    sec_n = ConnerSection(); sec_t, _ = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1.0, 40.0, 800))
    UFt, _ = ss_flutter(sec_t, np.linspace(1.0, 40.0, 800))
    print(f"nominal U_f {UF:.3f} m/s; TRUE plant U_f {UFt:.3f} m/s = {UFt/UF:.4f} U_f(nom)")
    print("all speeds below are U/U_f(nominal), the only U_f a test team would know\n")

    # ------------------------------------------------------------- 0
    print("=" * 86); print("0. THE AFS ON THE TRUE PLANT"); print("=" * 86)
    Ks, drops = {}, {}
    print(f"{'U/Uf':>6} {'zeta nominal':>13} {'zeta TRUE':>10} {'CL documented':>14}"
          f" {'CL as-built':>12} {'|z|max':>8} {'settle s':>9}")
    for r in RATIOS:
        U = r*UF
        Kd = lqg(sec_n, U); Ka = as_built(Kd); Ks[r] = (Ka, Kd)
        zn, _ = truth_mode(sec_n, U); zt, _ = truth_mode(sec_t, U)
        evd = np.linalg.eigvals(closed_loop_matrix(sec_t, U, Kd))
        eva = np.linalg.eigvals(closed_loop_matrix(sec_t, U, Ka))
        assert np.abs(eva).max() < 1, f"as-built AFS does not stabilise the truth at {r}"
        tau = TS/(1 - np.abs(eva).max())                 # slowest closed-loop decay
        drops[r] = max(DROP, int(np.ceil(5*tau/P)))
        print(f"{r:>6.2f} {zn:>+13.4f} {zt:>+10.4f} {least_damped(evd):>14.4f}"
              f" {least_damped(eva):>12.4f} {np.abs(eva).max():>8.5f} {drops[r]*P:>9.0f}")
    fine = np.arange(0.50, 1.2001, 0.01)
    worst = max(np.abs(np.linalg.eigvals(closed_loop_matrix(
        sec_t, r*UF, as_built(lqg(sec_n, r*UF))))).max() for r in fine)
    print(f"  as-built loop on the TRUE plant, 0.50-1.20 U_f in 0.01 steps: max |z| {worst:.5f}"
          + ("  stable everywhere" if worst < 1 else "  UNSTABLE SOMEWHERE"))
    print("  Near 1.05 the nominal model is past its flutter speed and the true plant")
    print("  is not, so a controller designed on the model fights an instability that")
    print("  is not there and REMOVES true damping. A measured plant is what fixes that.")

    # ------------------------------------------------------------- 1
    f_chk = np.linspace(BAND[0], BAND[1], 400)
    print("\n" + "=" * 86); print("1. IDENTITY CHECKS, NO NOISE"); print("=" * 86)
    for r in (0.8, 1.1):
        U = r*UF
        Pa, Ga, _, Ca = sampled_plant(sec_t, U)
        Gt = frf(Pa, Ga, Ca, f_chk)[:, 0, :]
        Ka, Kd = Ks[r]
        Kf = frf(*Ka, f_chk)[0]
        Tyr = Gt/(1 - (Kf*Gt).sum(0))                        # y = G (K y + r)
        Tur = 1/(1 - (Kf*Gt).sum(0))
        e_ind = np.abs(Tyr/(1 + (Kf*Tyr).sum(0)) - Gt).max()/np.abs(Gt).max()
        e_jio = np.abs(Tyr/Tur - Gt).max()/np.abs(Gt).max()
        zf, ff = flutter_mode(modal_fit(f_chk, Gt))
        zt, ft = truth_mode(sec_t, U)
        print(f"  U/Uf {r:.2f}: indirect(exact K) {e_ind:.1e}, joint I/O {e_jio:.1e} "
              f"relative; modal fit on exact FRF zeta {zf:+.5f} vs eig {zt:+.5f}, "
              f"f {ff:.3f} vs {ft:.3f} Hz")

    # ------------------------------------------------------------- 2
    print("\n" + "=" * 86); print("2. CLOSED-LOOP BACK-CALCULATION, LINEAR TRUTH"); print("=" * 86)
    jobs = [(sec_t, r*UF, Ks[r][0], s, True, "linear", None, drops[r] + 4, drops[r])
            for r in RATIOS for s in SEEDS]
    with ProcessPoolExecutor(max_workers=12) as pool:
        out = list(pool.map(run, jobs, chunksize=2))
    res = {(r, s): o for (_, U, _, s, *_), o, r in
           zip(jobs, out, [r for r in RATIOS for _ in SEEDS])}
    snr = np.median([o["snr"] for o in out], axis=0)
    print(f"  median achieved SNR (signal rms / sensor noise): "
          + ", ".join(f"{o} {20*np.log10(v):.0f} dB" for o, v in zip(OUT, snr)))
    print("  error = median over excited lines of |G_est - G_true|/|G_true|, median over"
          f" {len(SEEDS)} replicates")
    print(f"{'U/Uf':>6} " + " ".join(f"{n:>22}" for n in NAMES))
    print(f"{'':>6} " + " ".join(f"{'a      b      h':>22}" for _ in NAMES))
    modal = {n: [] for n in NAMES}
    for r in RATIOS:
        U = r*UF
        Ka, Kd = Ks[r]
        f = res[(r, SEEDS[0])]["f"]
        Pa, Ga, _, Ca = sampled_plant(sec_t, U)
        Gt = frf(Pa, Ga, Ca, f)[:, 0, :]
        errs = {n: [] for n in NAMES}; zeta = {n: [] for n in NAMES}
        for s in SEEDS:
            G = estimates(res[(r, s)], Ka, Kd)
            for n in NAMES:
                errs[n].append(np.median(np.abs(G[n] - Gt)/np.abs(Gt), axis=1))
                zeta[n].append(flutter_mode(modal_fit(f, G[n]))[0])
        print(f"{r:>6.2f} " + " ".join(
            " ".join(f"{v:6.3f}" for v in np.median(errs[n], axis=0)) + " " for n in NAMES))
        for n in NAMES:
            modal[n].append((np.median(zeta[n]), np.std(zeta[n], ddof=1)))

    print("\n  flutter-mode damping of the back-calculated plant (median, sd over replicates)")
    print(f"{'U/Uf':>6} {'truth':>8} " + " ".join(f"{n:>18}" for n in NAMES))
    truth_z = [truth_mode(sec_t, r*UF)[0] for r in RATIOS]
    for i, r in enumerate(RATIOS):
        print(f"{r:>6.2f} {truth_z[i]:>+8.4f} " + " ".join(
            f"{modal[n][i][0]:>+10.4f} ({modal[n][i][1]:.4f})" for n in NAMES))
    print(f"\n  flutter speed from the zeta = 0 crossing, U/U_f(nom); truth {UFt/UF:.4f}")
    print(f"    {'eig of the true plant':>22}: {zero_crossing(RATIOS, truth_z):.4f}")
    for n in NAMES:
        print(f"    {n:>22}: {zero_crossing(RATIOS, [m for m, _ in modal[n]]):.4f}")

    # ------------------------------------------------------------- 3
    print("\n" + "=" * 86); print("3. THE SAME TEST OPEN LOOP (u = r), FOR COMPARISON"); print("=" * 86)
    print("  settle five decay times of the true flutter mode (stage 3), where possible")
    ol = []
    for r in RATIOS:
        zt, ft = truth_mode(sec_t, r*UF)
        if zt <= 0:
            ol.append((r, None)); continue
        drop = max(DROP, int(np.ceil(5/(zt*2*np.pi*ft)/P)))
        ol.append((r, drop))
    jobs = [(sec_t, r*UF, Ks[r][0], s, False, "linear", None, d + 4, d)
            for r, d in ol if d is not None and d <= 40 for s in SEEDS]
    with ProcessPoolExecutor(max_workers=12) as pool:
        out = list(pool.map(run, jobs, chunksize=1))
    print(f"{'U/Uf':>6} {'settle s':>9} {'err a':>7} {'err b':>7} {'err h':>7}"
          f" {'zeta est':>10} {'truth':>8}")
    k = 0
    for r, d in ol:
        if d is None or d > 40:
            print(f"{r:>6.2f} {'-':>9}   open-loop test not possible: "
                  + ("the true plant is unstable" if d is None else "settling exceeds 5 min"))
            continue
        U = r*UF
        Pa, Ga, _, Ca = sampled_plant(sec_t, U)
        f = out[k]["f"]; Gt = frf(Pa, Ga, Ca, f)[:, 0, :]
        e, zz = [], []
        for _ in SEEDS:
            o = out[k]; k += 1
            G = o["Yp"].mean(0)/o["Up"].mean(0)
            e.append(np.median(np.abs(G - Gt)/np.abs(Gt), axis=1))
            zz.append(flutter_mode(modal_fit(f, G))[0])
        e = np.median(e, axis=0)
        print(f"{r:>6.2f} {d*P:>9.0f} {e[0]:>7.3f} {e[1]:>7.3f} {e[2]:>7.3f}"
              f" {np.median(zz):>+10.4f} {truth_mode(sec_t, U)[0]:>+8.4f}")

    # ------------------------------------------------------------- 4
    print("\n" + "=" * 86); print("4. WITH FREEPLAY IN THE TRUTH"); print("=" * 86)
    print("  the AFS was designed on the linear nominal model and never saw the gap.")
    print("  truth for the error is the TRUE STIFF linear plant.")
    fp_cases = (("symmetric gap", "hinge", 0.0), ("hinge preload p = 1", "hinge", 1.0))
    fp_ratios = (0.50, 0.80, 1.00, 1.10)
    jobs, keys = [], []
    for label, kind, lev in fp_cases:
        for r in fp_ratios:
            U = r*UF
            f0 = preload_forcing(sec_t, U, kind, lev)
            for s in SEEDS[:4]:
                jobs.append((sec_t, U, Ks[r][0], s, True, "freeplay", f0,
                             drops[r] + 4, drops[r]))
                keys.append((label, r))
    with ProcessPoolExecutor(max_workers=12) as pool:
        out = list(pool.map(run, jobs, chunksize=1))
    print(f"{'case':>20} {'U/Uf':>6} {'max|db|/d':>9} {'switch':>7} "
          + " ".join(f"{n:>11}" for n in NAMES) + "    <- beta-channel error")
    for label, kind, lev in fp_cases:
        for r in fp_ratios:
            os_ = [o for kk, o in zip(keys, out) if kk == (label, r)]
            if any(o["diverged"] for o in os_):
                print(f"{label:>20} {r:>6.2f}   diverged"); continue
            U = r*UF
            Ka, Kd = Ks[r]
            f = os_[0]["f"]
            Pa, Ga, _, Ca = sampled_plant(sec_t, U)
            Gt = frf(Pa, Ga, Ca, f)[:, 0, :]
            eb = {n: np.median([np.median(np.abs(estimates(o, Ka, Kd)[n][1] - Gt[1])/np.abs(Gt[1]))
                                for o in os_]) for n in NAMES}
            print(f"{label:>20} {r:>6.2f} {np.median([o['bmax'] for o in os_]):>9.2f}"
                  f" {int(np.median([o['switches'] for o in os_])):>7d} "
                  + " ".join(f"{eb[n]:>11.3f}" for n in NAMES))
    print("  with preload the flap stays on its stop (no switching) and the plant is")
    print("  exactly the stiff linear one; with a symmetric gap it is not.")


if __name__ == "__main__":
    main()
