"""Verification + dataset generation for the Conner-anchored Phase 1 harness."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.linalg import solve_lyapunov
from conner3dof import ConnerSection, state_space, ss_flutter
from phase1_conner import (simulate_linear, simulate_freeplay, excite_sweep,
                           excite_raps, build_dataset, input_matrix,
                           output_matrices, FS, NSAMP)

sec = ConnerSection()
UF, fF = ss_flutter(sec, np.linspace(1., 40., 3000), True)
print(f"Conner rig: U_f = {UF:.3f} m/s, f_f = {fF:.3f} Hz;  fs = {FS} Hz, "
      f"{NSAMP} pts = {NSAMP/FS:.2f} s")

# --- Hankel singular values with the REAL sensor set (3 RVDTs) -------------
print("\nHankel singular values, shaker -> 3 RVDTs, 8-state plant (normalized):")
for r in (0.30, 0.65, 0.90, 0.95):
    U = r*UF
    A, B = state_space(sec, U, True, "stiff"), input_matrix(sec, U)
    C, _ = output_matrices()
    Wc = solve_lyapunov(A, -B @ B.T)
    Wo = solve_lyapunov(A.T, -C.T @ C)
    h = np.sort(np.sqrt(np.abs(np.linalg.eigvals(Wc @ Wo))))[::-1]
    h /= h[0]
    print(f"  U/U_f={r:.2f}  " + "  ".join(f"{v:.2e}" for v in h))
print("  40 dB SNR corresponds to a relative floor of 1.0e-2")

# --- datasets --------------------------------------------------------------
ds_sweep = build_dataset(sec, UF, mode="sweep", seed=1, snr_db=40.0)
ds_raps = build_dataset(sec, UF, mode="raps", seed=1, snr_db=40.0)
flat = {}
for tag, ds in (("sweep", ds_sweep), ("raps", ds_raps)):
    for i, d in enumerate(ds):
        for k, v in d.items():
            if isinstance(v, (np.ndarray, float, int, str)):
                flat[f"{tag}_{i}_{k}"] = v
np.savez_compressed("phase1_conner_linear.npz", n=len(ds_sweep), **flat)

fp = {}
for gi, dg in enumerate((1.15, 1.83, 2.12)):
    for r in (0.35, 0.49, 0.73, 0.90):
        d = simulate_freeplay(sec, r*UF, np.deg2rad(dg), beta0_mult=3.0,
                              rng=np.random.default_rng(gi*10 + int(r*100)))
        for k, v in d.items():
            if isinstance(v, (np.ndarray, float, int, str)):
                fp[f"g{dg}_r{r}_{k}"] = v
np.savez_compressed("phase1_conner_freeplay.npz", **fp)

# --- figure ----------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
d = ds_sweep[-1]
ax[0].plot(d["t"], np.rad2deg(d["y"][0]), lw=.5, label=r"$\alpha$ [deg]")
ax[0].plot(d["t"], 100*d["y"][2], lw=.5, label="h [cm]")
ax[0].set_xlabel("t [s]"); ax[0].set_title(f"swept sine, $U/U_f$={d['ratio']}, 800 Hz")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

for dd in ds_sweep:
    f, P = welch(dd["y"][0], fs=FS, nperseg=1024)
    ax[1].semilogy(f, P, lw=.9, label=f"{dd['ratio']:.2f}")
ax[1].set_xlim(0, 30); ax[1].set_xlabel("f [Hz]"); ax[1].set_ylabel("PSD, pitch")
ax[1].set_title("pitch/plunge coalescence"); ax[1].legend(fontsize=7, title="$U/U_f$")
ax[1].grid(alpha=.3)

r = np.array([dd["ratio"] for dd in ds_sweep])
M = np.array([dd["modes"] for dd in ds_sweep])
for i, lab in enumerate(("mode 1", "mode 2 (flutters)", "mode 3 (flap)")):
    ax[2].plot(r, M[:, i, 1], "o-", ms=4, label=lab)
ax[2].axhline(0, color="k", lw=.6)
ax[2].set_xlabel("$U/U_f$"); ax[2].set_ylabel(r"$\zeta$ ground truth")
ax[2].set_title("what Phase 2 must recover"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig("phase1_conner.png", dpi=140)
print("\nwrote phase1_conner_linear.npz, phase1_conner_freeplay.npz, phase1_conner.png")
