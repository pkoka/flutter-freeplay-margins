"""
Stage 2. Measure model uncertainty from many simulated tunnel runs, regress it
against operating condition, and test whether the memo's absolute-stability
certificate can consume the result.

  1  perturbed nominal, so dA and dB are nonzero for the first time
  2  grid of runs across airspeed and shaker amplitude, one run = one entry
  3  in-run screen using per-period scatter, no ground truth required
  4  mismatch regressed against condition
  5  W built from the in-gap runs against the FREE-region nominal
  6  robust Popov certificate, and whether it holds on the true plant

WHY THE FREE NOMINAL IN STEP 5. The memo's sector bound already covers the
deadband exactly. A W measured against the STIFF nominal also contains the
freeplay, so the freeplay is counted twice and the second count destroys the
certificate. Below LCO onset and below the gap the flap never switches, the
true plant is EXACTLY the free-region linear system, and the comparison
isolates the linear model error the Popov test actually needs.

Run once to build grid.pkl, thereafter it is reused.
"""
import numpy as np, pickle, os, time
from scipy.optimize import brentq
from conner3dof import ConnerSection, ss_flutter
from absolute_stability import loop, G
from freeplay import integrate
from conner3dof import state_space
from stage2a_harness import perturbed_section, periodic_run, frf_model, OUT

CACHE = "grid.pkl"
RATIOS = (0.10, 0.16, 0.22, 0.28, 0.35, 0.45, 0.60)
AMPS = (0.05, 0.2, 0.8, 3.0, 12.0, 50.0)
W_RAD = np.geomspace(1e-2, 1e3, 2000)
DELTA = np.deg2rad(2.12)


def build_grid(sec_n, sec_t, UF):
    rows = []
    t0 = time.time()
    for r in RATIOS:
        U = r*UF
        for a in AMPS:
            R = periodic_run(sec_t, U, DELTA, a, P=16.0, nper=3, drop=1, seed=101)
            f, H = R["f"], R["H"]
            rows.append(dict(
                r=r, U=U, amp=a, f=f, H=H, sw=R["switches"],
                bm=R["beta_max_over_delta"], scatter=R["scatter"],
                Hn_stiff=frf_model(sec_n, U, f, "stiff"),
                Hn_free=frf_model(sec_n, U, f, "free"),
                Ht_stiff=frf_model(sec_t, U, f, "stiff"),
                Ht_free=frf_model(sec_t, U, f, "free")))
    print(f"  {len(rows)} runs in {time.time()-t0:.0f} s")
    return rows


# --------------------------------------------------------------------------
# a usable weight: floor at antiresonances, then a fitted low-order bound
# --------------------------------------------------------------------------
FLOOR = 0.05


def mismatch_floor(H, Hr, floor=FLOOR):
    """|H - Hr| relative to |Hr|, except where |Hr| falls below floor*max|Hr|
    (an antiresonance), where it is measured against that floor instead.
    Pure relative error is meaningless where |Hr| -> 0, and those bins alone
    set the infinity norm of a raw per-bin W (max 6.9 in September)."""
    den = np.maximum(np.abs(Hr), floor*np.abs(Hr).max(axis=-1, keepdims=True))
    return np.abs(H - Hr)/den


def fit_weight(f, Wdata, n=36):
    """Lowest first-order W(s) = k (s/w1 + 1)/(s/w2 + 1), w1, w2 > 0 (stable,
    minimum phase), bounding Wdata at every measured line, chosen to minimise
    the mean log magnitude over the band. Returns W evaluated on rad/s."""
    w = 2*np.pi*f
    best = None
    for w1 in np.geomspace(2*np.pi*0.2, 2*np.pi*300, n):
        for w2 in np.geomspace(2*np.pi*0.2, 2*np.pi*300, n):
            shape = np.abs((1j*w/w1 + 1)/(1j*w/w2 + 1))
            k = np.max(Wdata/shape)
            cost = np.mean(np.log(k*shape))
            if best is None or cost < best[0]:
                best = (cost, k, w1, w2)
    _, k, w1, w2 = best
    return lambda wr: k*np.abs((1j*wr/w1 + 1)/(1j*wr/w2 + 1)), (k, w1/(2*np.pi), w2/(2*np.pi))


def hinge_frf(sec, U, f, region="free"):
    """FRF from flap hinge moment to the three RVDTs: the channel the
    robust-stability test is about."""
    A = state_space(sec, U, True, region)
    b = np.zeros((8, 1)); b[3:6, 0] = loop(sec, U)[1][3:6, 0]
    Cm = np.zeros((3, 8)); Cm[0, 0] = Cm[1, 1] = Cm[2, 2] = 1.0
    return np.stack([(Cm @ np.linalg.solve(2j*np.pi*fi*np.eye(8) - A, b))[:, 0]
                     for fi in f], axis=-1)


def measure_weight(pct, UF, conds=((0.10, 0.05), (0.10, 0.2), (0.16, 0.05),
                                   (0.16, 0.2), (0.22, 0.05), (0.22, 0.2)), seed=101):
    """W from in-gap runs on the pct-perturbed truth against the FREE nominal,
    max over channels and runs, with the antiresonance floor. Also returns the
    true hinge-channel error at the same lines, to check the shaker-channel
    proxy, and how many runs stayed in the gap."""
    sec_n = ConnerSection(); sec_t, _ = perturbed_section(pct)
    Ws, Ts, kept = [], [], 0
    for r, a in conds:
        R = periodic_run(sec_t, r*UF, DELTA, a, P=16.0, nper=3, drop=1, seed=seed)
        if R["switches"] > 0:
            continue
        kept += 1
        f = R["f"]
        Ws.append(mismatch_floor(R["H"], frf_model(sec_n, r*UF, f, "free")).max(0))
        Ts.append(mismatch_floor(hinge_frf(sec_t, r*UF, f), hinge_frf(sec_n, r*UF, f)).max(0))
    return f, np.max(Ws, axis=0), np.max(Ts, axis=0), kept


def popov_robust(sec, U, Wv):
    """max over eta>=0 of  min_w { Re[(1+j w eta)G] - |1+j w eta||G|W(w) } + 1/k.
    W=0 recovers the memo's Popov test exactly."""
    A, B, C, k, _ = loop(sec, U)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    g = G(A, B, C, W_RAD)
    def m(eta):
        s = 1 + 1j*W_RAD*eta
        return float(np.min((s*g).real - np.abs(s)*np.abs(g)*Wv) + 1/k)
    return float(max(m(e) for e in
                     np.concatenate([[0.0], np.geomspace(1e-6, 1e8, 200)])))


def certified(sec, Wv, lo=0.05, hi=6.9):
    if popov_robust(sec, lo, Wv) <= 0: return np.nan
    if popov_robust(sec, hi, Wv) > 0: return hi
    return brentq(lambda U: popov_robust(sec, U, Wv), lo, hi, xtol=1e-4)


def main():
    sec_n = ConnerSection()
    sec_t, fac = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1.0, 40.0, 800))
    if os.path.exists(CACHE):
        rows = pickle.load(open(CACHE, "rb"))["rows"]
        print(f"loaded {len(rows)} cached runs")
    else:
        rows = build_grid(sec_n, sec_t, UF)
        pickle.dump(dict(rows=rows, UF=UF, fac=fac), open(CACHE, "wb"))

    rel = lambda H, Hr: np.abs(H - Hr)/np.abs(Hr)
    for w in rows:
        w["scat"] = float(np.median(np.median(
            np.abs(w["scatter"])/np.abs(w["H"]), axis=1)))

    # ---------------------------------------------------------------- 3
    print("\n" + "="*74); print("3. IN-RUN SCREEN (per-period scatter only)"); print("="*74)
    keep = [w for w in rows if w["scat"] <= 0.05]
    bad = [w for w in rows if w["scat"] > 0.05]
    print(f"  {len(keep)} of {len(rows)} runs pass at scatter <= 0.05")
    print(f"  median alpha mismatch, kept {np.median([np.median(rel(w['H'],w['Hn_stiff'])[0]) for w in keep]):.2f}"
          f"  rejected {np.median([np.median(rel(w['H'],w['Hn_stiff'])[0]) for w in bad]):.2f}")

    # ---------------------------------------------------------------- 4
    print("\n" + "="*74); print("4. MISMATCH REGRESSED AGAINST CONDITION"); print("="*74)
    print("  log W = c0 + c1 log(beta_max/delta) + c2 (U/U_f)")
    X = np.column_stack([np.ones(len(keep)), np.log([w["bm"] for w in keep]),
                         [w["r"] for w in keep]])
    for i, o in enumerate(OUT):
        yv = np.log([np.median(rel(w["H"], w["Hn_stiff"])[i]) for w in keep])
        c, *_ = np.linalg.lstsq(X, yv, rcond=None)
        r2 = 1 - np.sum((yv - X@c)**2)/np.sum((yv - yv.mean())**2)
        print(f"  {o:6s} c0 {c[0]:+.3f}  c1 {c[1]:+.3f}  c2 {c[2]:+.3f}   R2 {r2:.3f}")

    # ---------------------------------------------------------------- 5
    print("\n" + "="*74); print("5. W FROM THE IN-GAP RUNS, FREE NOMINAL"); print("="*74)
    ingap = [w for w in rows if w["sw"] == 0]
    print(f"  {len(ingap)} runs with zero switching events: the true plant is")
    print(f"  exactly the free-region linear system in these.")
    print(f"{'chan':>7} {'vs FREE nom':>13} {'TRUTH free pair':>17} {'ratio':>7}"
          f" {'vs STIFF nom (double counts)':>30}")
    for i, o in enumerate(OUT):
        m = np.median([np.median(rel(w["H"], w["Hn_free"])[i]) for w in ingap])
        t = np.median([np.median(rel(w["Ht_free"], w["Hn_free"])[i]) for w in ingap])
        s = np.median([np.median(rel(w["H"], w["Hn_stiff"])[i]) for w in ingap])
        print(f"{o:>7} {m:>13.4f} {t:>17.4f} {m/t:>7.2f} {s:>30.4f}")
    print("  ratio > 1 means the measurement is conservative, the safe direction.")
    print("  it is the 40 dB noise floor: at 120 dB the beta channel lands exactly.")
    fgrid = ingap[0]["f"]
    Wenv = np.max(np.array([rel(w["H"], w["Hn_free"]) for w in ingap]), axis=0)[0]

    def W_of(w_rad, chan):
        f = w_rad/(2*np.pi)
        out = np.interp(f, fgrid, chan)
        out[(f < fgrid[0]) | (f > fgrid[-1])] = chan.max()   # no data outside band
        return out

    # ---------------------------------------------------------------- 6
    print("\n" + "="*74); print("6. ROBUST POPOV, AND WHAT IT CAN ABSORB"); print("="*74)
    z = np.zeros_like(W_RAD)
    Un, Ut = certified(sec_n, z), certified(sec_t, z)
    print(f"  nominal model, W = 0            {Un/UF:.4f} U_f   (memo's number)")
    print(f"  TRUE model,    W = 0            {Ut/UF:.4f} U_f   (unknowable in practice)")
    Ur = certified(sec_n, W_of(W_RAD, Wenv))
    print(f"  nominal + measured W(w)         "
          + (f"{Ur/UF:.4f} U_f" if np.isfinite(Ur) else "nothing certified"))
    print("\n  how much flat multiplicative uncertainty can the certificate take?")
    for Wf in (0.0, 0.01, 0.02, 0.03, 0.04, 0.05):
        U = certified(sec_n, Wf*np.ones_like(W_RAD))
        print(f"    W = {Wf:.2f}  " + (f"{U/UF:.4f} U_f" if np.isfinite(U) else "nothing certified"))
    A, B, C, k, _ = loop(sec_n, Un)
    g = G(A, B, C, W_RAD); j = np.argmin(g.real)
    print(f"\n  binding frequency {W_RAD[j]/2/np.pi:.3f} Hz (inside the 0.5-30 Hz tested band)")
    print(f"  |G| there {np.abs(g[j]):.2f} against 1/k = {1/k:.4f}, so W = 0.05 alone")
    print(f"  contributes {np.abs(g[j])*0.05:.2f}, six times the whole budget.")

    print("\n" + "="*74); print("END TO END"); print("="*74)
    onset = None
    for r in np.arange(0.25, 0.30, 0.01):
        _, zz, _ = integrate(sec_t, r*UF, DELTA, T=30.0,
                             z0=[0, 3*DELTA, 0, 0, 0, 0, 0, 0], record_from=25.0)
        if np.std(zz[:, 1])/DELTA > 1e-3 and onset is None:
            onset = r
    print(f"  TRUE plant LCO onset          {onset:.3f} U_f")
    print(f"  nominal certificate, no W     {Un/UF:.4f} U_f  "
          f"{'VALID' if Un/UF < onset else 'VIOLATED'}")
    print(f"  robust certificate, measured W " +
          (f"{Ur/UF:.4f} U_f" if np.isfinite(Ur) else "nothing certified"))
    print("\n  sweep the perturbation: is the NON-robust certificate ever violated?")
    print(f"{'pct':>5} {'true cert':>10} {'true onset':>11} {'verdict':>10}")
    for pct in (2.0, 8.0, 32.0):
        st, _ = perturbed_section(pct)
        on = np.nan
        for r in np.arange(0.10, 0.45, 0.02):
            _, zz, _ = integrate(st, r*UF, DELTA, T=25.0,
                                 z0=[0, 3*DELTA, 0, 0, 0, 0, 0, 0], record_from=20.0)
            if np.std(zz[:, 1])/DELTA > 1e-3:
                on = r; break
        print(f"{pct:>5.0f} {certified(st, z)/UF:>10.4f} {on:>11.3f} "
              f"{'valid' if Un/UF < on else 'VIOLATED':>10}")
    print("\n  CONCLUSION. Popov's whole robustness budget is about 3 percent")
    print("  multiplicative uncertainty. Real model error here is 3-11 percent.")
    print("  So any honest W destroys the certificate, while the non-robust")
    print("  certificate survives 32 percent parameter error because Popov is")
    print("  already conservative by a wide margin. The measured W and this")
    print("  certificate are not in series. W belongs to an H-infinity design.")


if __name__ == "__main__":
    main()
