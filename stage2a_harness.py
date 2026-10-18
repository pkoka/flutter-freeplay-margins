"""
Stage 2a. Perturbed-nominal multi-run harness and its gate test.

WHY. Every earlier phase used a nominal model that was exact by construction,
so the only model error was the freeplay. Beat 3 of the Ting call sheet says
that is the case that does not hold near flutter: A and B are themselves wrong
there. Here the nominal is deliberately wrong, so the measured mismatch
contains both linear model error and the nonlinearity, which is the situation
the method has to survive.

GATE for this stage. Drive hard enough that the freeplay is negligible, then
check the harness recovers the KNOWN linear model error. If it cannot do that,
no mismatch it reports at any other condition means anything.

TWO BUGS THIS SCRIPT FIXES, both of which were inflating the measured mismatch:

  1. SAMPLE RATE. step = int(round(1/(FS*dt))) with FS=800, dt=2e-4 gives 6,
     so the true rate is 833.33 Hz and the frequency axis is 4% wrong. On
     lightly damped modes that alone produces large errors near resonance.
     Same construction appears in phase1_conner.simulate_freeplay and
     margin_amplitude.rap_record in the repo. Fixed by dt = 1.25e-4 -> step 10.

  2. ESTIMATOR. A multisine is deterministic and periodic. Welch/H1 with
     arbitrary segments smears the excitation lines and was the DOMINANT
     artifact, inflating mismatch 3-4x. Correct practice: integer periods,
     drop the transient period, DFT at the excited lines only, average the
     complex spectra over the kept periods.
"""
import numpy as np
from conner3dof import ConnerSection, state_space, ss_flutter
from phase1_conner import input_matrix, output_matrices
from freeplay import integrate

FS, DT_SIM = 800.0, 1.25e-4          # step = 10 exactly
BAND = (0.5, 30.0)
OUT = ["alpha", "beta", "h"]
PERTURBED = ["Ka", "Kb", "Kh", "Ia", "Ib", "Sa", "Sb", "M", "za", "zb", "zh"]


def perturbed_section(pct=4.0, seed=7):
    """A physically meaningful wrong model: stiffnesses, inertias and damping
    a few percent off the drawing. Rebuilds every dependent quantity the
    constructor derives from them."""
    rng = np.random.default_rng(seed)
    sec = ConnerSection()
    fac = {k: 1.0 + pct/100*rng.standard_normal() for k in PERTURBED}
    for k, v in fac.items():
        setattr(sec, k, getattr(sec, k)*v)
    b, M = sec.b, sec.M
    sec.Mh = M + 2*0.47485/sec.span
    sec.Ca, sec.Cb, sec.Ch = sec.Ka*M*b**2, sec.Kb*M*b**2, sec.Kh*M
    sec.Ms = np.array([
        [sec.Ia, sec.Ib + b*(sec.c - sec.a)*sec.Sb, sec.Sa],
        [sec.Ib + b*(sec.c - sec.a)*sec.Sb, sec.Ib, sec.Sb],
        [sec.Sa, sec.Sb, sec.Mh]])
    sec.Ks = np.diag([sec.Ca, sec.Cb, sec.Ch])
    sec.kappa = np.pi*sec.rho*b**2/M
    return sec, fac


def frf_model(sec, U, f, region="stiff"):
    A = state_space(sec, U, True, region)
    B = input_matrix(sec, U)
    C, _ = output_matrices()
    H = np.empty((3, len(f)), complex)
    for j, wj in enumerate(2*np.pi*f):
        H[:, j] = (C @ np.linalg.solve(1j*wj*np.eye(8) - A, B))[:, 0]
    return H


def periodic_run(sec, U, delta, amp, P=16.0, nper=4, drop=1, seed=1,
                 snr_db=40.0, band=BAND):
    """One simulated tunnel entry. Periodic random-phase multisine, integer
    periods, FRF from the DFT at the excited lines. Returns the per-period
    scatter as well, which is the within-run repeatability."""
    rng = np.random.default_rng(seed)
    B = input_matrix(sec, U)
    C, _ = output_matrices()
    n, df = int(P*FS), 1.0/P
    ks = np.arange(int(np.ceil(band[0]/df)), int(np.floor(band[1]/df)) + 1)
    ph = rng.uniform(0, 2*np.pi, len(ks))
    tt = np.arange(n*nper)/FS
    uu = amp*sum(np.cos(2*np.pi*k*df*tt + p)
                 for k, p in zip(ks, ph))/np.sqrt(len(ks))
    uf = lambda s: np.interp(s, tt, uu)
    t, z, ns = integrate(sec, U, delta, T=P*nper, dt=DT_SIM,
                         u_func=uf, B=B, record_from=0.0)
    step = int(round(1/(FS*DT_SIM)))
    assert abs(1/(step*DT_SIM) - FS) < 1e-9, "sample rate is not FS"
    idx = np.arange(0, len(t), step)[:len(tt)]
    y = (C @ z.T)[:, idx]
    u = uu[:len(idx)]
    y = y + rng.standard_normal(y.shape)*np.atleast_2d(
        np.std(y, axis=1)*10**(-snr_db/20)).T
    Hs = []
    for p in range(drop, nper):
        sl = slice(p*n, (p + 1)*n)
        Uk = np.fft.rfft(u[sl])[ks]
        Yk = np.fft.rfft(y[:, sl], axis=1)[:, ks]
        Hs.append(Yk/Uk)
    H = np.mean(Hs, axis=0)
    return dict(f=ks*df, H=H, scatter=np.std(Hs, axis=0), U=U, amp=amp,
                switches=ns, beta_max_over_delta=np.abs(z[:, 1]).max()/delta)


def gate(pct=4.0, ratios=(0.40, 0.85), amps=(16.0, 256.0), delta_deg=2.12):
    sec_nom = ConnerSection()
    sec_true, fac = perturbed_section(pct)
    UF, _ = ss_flutter(sec_nom, np.linspace(1.0, 40.0, 800))
    UFt, _ = ss_flutter(sec_true, np.linspace(1.0, 40.0, 800))
    delta = np.deg2rad(delta_deg)
    print("PERTURBATION  " + "  ".join(f"{k}:{v:.3f}" for k, v in fac.items()))
    print(f"  U_f nominal {UF:.3f}, true {UFt:.3f}  ({100*(UFt-UF)/UF:+.2f}%)")
    Un, Ut = 0.73*UF, 0.73*UF
    dA = np.linalg.norm(state_space(sec_true, Ut, True, "stiff")
                        - state_space(sec_nom, Un, True, "stiff"))
    print(f"  ||dA||_F/||A||_F {dA/np.linalg.norm(state_space(sec_nom,Un,True,'stiff')):.4f}"
          f"   ||dB||/||B|| "
          f"{np.linalg.norm(input_matrix(sec_true,Ut)-input_matrix(sec_nom,Un))/np.linalg.norm(input_matrix(sec_nom,Un)):.4f}\n")
    for r in ratios:
        U = r*UF
        for amp in amps:
            R = periodic_run(sec_true, U, delta, amp)
            f, H = R["f"], R["H"]
            Ht, Hn = frf_model(sec_true, U, f), frf_model(sec_nom, U, f)
            row = lambda M, D: "  ".join(
                f"{o} {np.median((np.abs(M)/np.abs(D))[i]):.3f}"
                for i, o in enumerate(OUT))
            print(f"U/Uf {r:.2f}  amp {amp:6.0f}  max|beta|/delta {R['beta_max_over_delta']:5.1f}"
                  f"  ({R['switches']} switches)")
            print("   meas vs TRUE model : " + row(H - Ht, Ht))
            print("   meas vs NOM  model : " + row(H - Hn, Hn))
            print("   TRUE vs NOM  target: " + row(Ht - Hn, Hn))
            print("   per-period scatter : " + row(R["scatter"], H))


if __name__ == "__main__":
    gate()
