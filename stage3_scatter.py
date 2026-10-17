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
"""
import numpy as np
from conner3dof import ConnerSection, ss_flutter, state_space
from stage2a_harness import perturbed_section, periodic_run, frf_model, OUT

CONDS = ((0.10, 0.2), (0.16, 0.8), (0.22, 0.2), (0.28, 0.8))
SEEDS = (11, 23, 37, 51, 67, 83)


def main():
    sec_n = ConnerSection(); sec_t, _ = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1.0, 40.0, 800))
    d = np.deg2rad(2.12)
    print(f"{len(CONDS)} conditions x {len(SEEDS)} replicates, 40 dB, 16-bit\n")
    print(f"{'U/Uf':>6} {'amp':>5} {'chan':>6} {'mean W':>8} {'sd':>8} {'sd/mean':>8}"
          f" {'truth':>8} {'signal/floor':>13}")
    allrows = []
    for r, a in CONDS:
        U = r*UF
        meas = []
        for s in SEEDS:
            R = periodic_run(sec_t, U, d, a, P=16.0, nper=3, drop=1, seed=s)
            Hfn = frf_model(sec_n, U, R["f"], "free")
            meas.append([np.median((np.abs(R["H"]-Hfn)/np.abs(Hfn))[i]) for i in range(3)])
            Hft = frf_model(sec_t, U, R["f"], "free")
            truth = [np.median((np.abs(Hft-Hfn)/np.abs(Hfn))[i]) for i in range(3)]
        meas = np.array(meas)
        for i, o in enumerate(OUT):
            m, sd = meas[:, i].mean(), meas[:, i].std(ddof=1)
            print(f"{r:>6.2f} {a:>5.1f} {o:>6} {m:>8.4f} {sd:>8.4f} {sd/m:>8.3f}"
                  f" {truth[i]:>8.4f} {truth[i]/sd if sd > 0 else np.inf:>13.1f}")
            allrows.append((r, a, o, m, sd, truth[i]))
        print()
    print("signal/floor is the true model error divided by the replicate scatter.")
    print("Above ~3 the measured W is reporting model error. Near 1 it is reporting")
    print("the measurement chain, and no amount of averaging over conditions fixes")
    print("that, because the scatter is what sets the resolution of the whole method.")
    sig = [t/sd for *_, sd, t in [(r, a, o, m, sd, t) for r, a, o, m, sd, t in allrows] if sd > 0]
    print(f"\nmedian signal/floor across all condition-channel pairs: {np.median(sig):.1f}")


def collapse_sweep(speeds=(0.240, 0.255, 0.265, 0.272, 0.278, 0.281, 0.284, 0.2855),
                   seeds=(11, 23, 37, 51), target_bm=0.6):
    """Where does measurement quality collapse as the free-region instability
    at 0.286 U_f is approached? Reported against the least-damped pole's zeta,
    which transfers to any rig, not just against U/U_f."""
    sec_n = ConnerSection(); sec_t, _ = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1.0, 40.0, 800))
    d = np.deg2rad(2.12)
    print("\n" + "=" * 82)
    print("COLLAPSE OF MEASUREMENT QUALITY APPROACHING THE INSTABILITY")
    print("=" * 82)
    print("  amplitude is auto-scaled at each speed to keep the flap inside the")
    print("  gap, so the true plant stays exactly free-region linear throughout.")
    print(f"{'U/Uf':>7} {'zeta_min':>9} {'amp':>8} {'bmax/d':>7} "
          + " ".join(f"{'S/F '+o:>10}" for o in OUT))
    for r in speeds:
        U = r*UF
        ev = np.linalg.eigvals(state_space(sec_n, U, True, "free"))
        osc = ev[np.abs(ev.imag) > 1e-6]
        zmin = float(np.min(-osc.real/np.abs(osc)))
        amp = 0.5
        for _ in range(3):                       # auto-scale into the gap
            R = periodic_run(sec_t, U, d, amp, P=16.0, nper=3, drop=1, seed=seeds[0])
            bm = R["beta_max_over_delta"]
            if bm < 0.9: break
            amp *= target_bm/bm
        meas = []
        for sd_ in seeds:
            R = periodic_run(sec_t, U, d, amp, P=16.0, nper=3, drop=1, seed=sd_)
            Hfn = frf_model(sec_n, U, R["f"], "free")
            meas.append([np.median((np.abs(R["H"]-Hfn)/np.abs(Hfn))[i]) for i in range(3)])
        Hft = frf_model(sec_t, U, R["f"], "free")
        truth = [np.median((np.abs(Hft-Hfn)/np.abs(Hfn))[i]) for i in range(3)]
        meas = np.array(meas)
        sf = [truth[i]/meas[:, i].std(ddof=1) if meas[:, i].std(ddof=1) > 0 else np.inf
              for i in range(3)]
        print(f"{r:>7.4f} {zmin:>9.5f} {amp:>8.3f} {bm:>7.3f} "
              + " ".join(f"{v:>10.2f}" for v in sf))
    print("\n  signal/floor above 3 means the measurement resolves model error.")


if __name__ == "__main__":
    main()
    collapse_sweep()
