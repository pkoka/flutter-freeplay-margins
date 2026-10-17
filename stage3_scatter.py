"""
Stage 3. Run-to-run scatter floor.

Model error and run-to-run scatter are indistinguishable in a single run. The
spread across replicates at a fixed condition is the floor; anything below it
is not model error. Replicates here re-draw the multisine phases and the
measurement noise, which is the simulation analogue of re-running the same
tunnel point. A real campaign would also break and remake the mount, which
this cannot capture and which will be larger.

Conditions are chosen in-gap and below LCO onset, so the true plant is exactly
the free-region linear system and the only thing moving between replicates is
the measurement chain.

TWO CORRECTIONS to the September version of the collapse sweep:

  1. ZETA OF THE RIGHT PLANT. The data come from the perturbed TRUE plant, but
     zeta_min was taken from the NOMINAL model. Near onset the two differ by
     more than 2x (0.284 U_f: nominal 0.0011, true 0.0024), so the cliff was
     quoted against a plant that was not the one being measured. Both are
     printed now; the true-plant zeta is the one that describes the test.
  2. SETTLING. The estimator drops one 16 s period before averaging. Near onset
     the flutter mode's decay time 1/(zeta w) passes 10 s, so start-up
     transients are still in the averaged periods. collapse_sweep(settle_taus=5)
     drops enough periods to cover five decay times of the true plant, which
     separates "the measurement chain cannot resolve the plant" from "the test
     did not wait long enough".

Four replicates was too thin a basis for a standard deviation; twelve now.
"""
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from conner3dof import ConnerSection, ss_flutter, state_space
from stage2a_harness import perturbed_section, periodic_run, frf_model, OUT

CONDS = ((0.10, 0.2), (0.16, 0.8), (0.22, 0.2), (0.28, 0.8))
SEEDS = (11, 23, 37, 51, 67, 83, 97, 109, 127, 139, 151, 163)
P = 16.0


def slowest_mode(sec, U):
    """(zeta, f_Hz, tau_s) of the least-damped free-region mode at U."""
    ev = np.linalg.eigvals(state_space(sec, U, True, "free"))
    osc = ev[np.abs(ev.imag) > 1e-6]
    i = np.argmin(-osc.real/np.abs(osc))
    w = abs(osc[i]); z = -osc[i].real/w
    return float(z), w/(2*np.pi), (1/(z*w) if z > 0 else np.inf)


def _run(args):
    sec, U, d, amp, nper, drop, seed = args
    return periodic_run(sec, U, d, amp, P=P, nper=nper, drop=drop, seed=seed)


def replicates(pool, sec_t, U, d, amp, nper=3, drop=1, seeds=SEEDS):
    return list(pool.map(_run, [(sec_t, U, d, amp, nper, drop, s) for s in seeds]))


def mismatch(R, sec_n, U):
    Hfn = frf_model(sec_n, U, R["f"], "free")
    return [np.median((np.abs(R["H"] - Hfn)/np.abs(Hfn))[i]) for i in range(3)]


def truth_of(sec_t, sec_n, U, f):
    Hft, Hfn = frf_model(sec_t, U, f, "free"), frf_model(sec_n, U, f, "free")
    return [np.median((np.abs(Hft - Hfn)/np.abs(Hfn))[i]) for i in range(3)]


def main(pool):
    sec_n = ConnerSection(); sec_t, _ = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1.0, 40.0, 800))
    d = np.deg2rad(2.12)
    print(f"{len(CONDS)} conditions x {len(SEEDS)} replicates, 40 dB, 16-bit\n")
    print(f"{'U/Uf':>6} {'amp':>5} {'chan':>6} {'mean W':>8} {'sd':>8} {'sd/mean':>8}"
          f" {'truth':>8} {'signal/floor':>13}")
    allrows = []
    for r, a in CONDS:
        U = r*UF
        Rs = replicates(pool, sec_t, U, d, a)
        meas = np.array([mismatch(R, sec_n, U) for R in Rs])
        truth = truth_of(sec_t, sec_n, U, Rs[0]["f"])
        for i, o in enumerate(OUT):
            m, sd = meas[:, i].mean(), meas[:, i].std(ddof=1)
            print(f"{r:>6.2f} {a:>5.1f} {o:>6} {m:>8.4f} {sd:>8.4f} {sd/m:>8.3f}"
                  f" {truth[i]:>8.4f} {truth[i]/sd if sd > 0 else np.inf:>13.1f}")
            allrows.append((r, a, o, m, sd, truth[i]))
        zt, _, tau = slowest_mode(sec_t, U)
        bm = max(R["beta_max_over_delta"] for R in Rs)
        sw = max(R["switches"] for R in Rs)
        print(f"{'':>6} true-plant zeta {zt:.5f}, decay time {tau:.1f} s, "
              f"max|beta|/delta {bm:.2f}, max switches {sw}"
              + ("   <- flap left the gap: freeplay, not linear" if sw > 0 else ""))
        print()
    print("signal/floor is the true model error divided by the replicate scatter.")
    print("Above ~3 the measured W is reporting model error. Near 1 it is reporting")
    print("the measurement chain, and no amount of averaging over conditions fixes")
    print("that, because the scatter is what sets the resolution of the whole method.")
    sig = [t/sd for *_, sd, t in allrows if sd > 0]
    print(f"\nmedian signal/floor across all condition-channel pairs: {np.median(sig):.1f}")


def collapse_sweep(pool, speeds=(0.240, 0.255, 0.265, 0.272, 0.278, 0.281, 0.284,
                                 0.2855, 0.2865, 0.2875),
                   seeds=SEEDS, target_bm=0.6, settle_taus=None):
    """Where does measurement quality collapse as the free-region instability
    is approached? Reported against the least-damped pole's zeta, which
    transfers to any rig, not just against U/U_f.

    settle_taus=None keeps the September protocol (drop one period). A number
    drops enough periods to cover that many decay times of the TRUE plant's
    slowest mode, which an experimenter would approximate by waiting until the
    per-period scatter stops falling."""
    sec_n = ConnerSection(); sec_t, _ = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1.0, 40.0, 800))
    d = np.deg2rad(2.12)
    print("\n" + "=" * 100)
    print("COLLAPSE OF MEASUREMENT QUALITY APPROACHING THE INSTABILITY   "
          + ("drop 1 period (September protocol)" if settle_taus is None
             else f"settle {settle_taus:g} true-plant decay times"))
    print("=" * 100)
    print("  amplitude is auto-scaled at each speed to keep the flap inside the")
    print("  gap, so the true plant stays exactly free-region linear throughout.")
    print(f"{'U/Uf':>7} {'z_nom':>8} {'z_true':>8} {'tau s':>6} {'drop':>4} {'bmax/d':>6} "
          + " ".join(f"{'S/F '+o:>9}" for o in OUT)
          + "  " + " ".join(f"{'m/t '+o:>9}" for o in OUT))
    for r in speeds:
        U = r*UF
        zn, _, _ = slowest_mode(sec_n, U)
        zt, _, tau = slowest_mode(sec_t, U)
        drop = 1 if settle_taus is None else max(1, int(np.ceil(settle_taus*tau/P)))
        nper = drop + 2
        amp = 0.5
        for _ in range(3):                       # auto-scale into the gap
            R = _run((sec_t, U, d, amp, nper, drop, seeds[0]))
            bm = R["beta_max_over_delta"]
            if bm < 0.9: break
            amp *= target_bm/bm
        Rs = replicates(pool, sec_t, U, d, amp, nper, drop, seeds)
        bm = max(R["beta_max_over_delta"] for R in Rs)
        meas = np.array([mismatch(R, sec_n, U) for R in Rs])
        truth = truth_of(sec_t, sec_n, U, Rs[0]["f"])
        sf = [truth[i]/meas[:, i].std(ddof=1) for i in range(3)]
        mt = [meas[:, i].mean()/truth[i] for i in range(3)]
        print(f"{r:>7.4f} {zn:>8.5f} {zt:>8.5f} {tau:>6.1f} {drop:>4d} {bm:>6.2f} "
              + " ".join(f"{v:>9.2f}" for v in sf) + "  "
              + " ".join(f"{v:>9.2f}" for v in mt))
    print("\n  S/F above 3 means the replicate scatter resolves model error.")
    print("  m/t is mean measured mismatch over true mismatch: 1 is unbiased, and")
    print("  the 40 dB noise floor alone puts it near 2 to 3 (stage 2a gate).")


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=12) as pool:
        main(pool)
        collapse_sweep(pool)
        collapse_sweep(pool, settle_taus=5.0)
