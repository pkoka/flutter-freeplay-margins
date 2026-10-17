"""Path-following sweep reproducing Conner Fig. 8. Run: python freeplay_sweep.py"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from conner3dof import ConnerSection, ss_flutter
from freeplay import integrate, rms_amplitudes, dominant_freq

sec = ConnerSection()
UF, _ = ss_flutter(sec, np.linspace(1., 40., 3000), True)
gaps = {"Gap 1": 1.15, "Gap 2": 1.83, "Gap 3": 2.12}
ratios = np.arange(0.10, 0.96, 0.025)

fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
out = {}
for name, dg in gaps.items():
    d = np.deg2rad(dg)
    A, B, H, F = [], [], [], []
    z0 = [0, d, 0, 0, 0, 0, 0, 0]
    for r in ratios:
        t, z, ns = integrate(sec, r*UF, d, T=14.0, dt=1.25e-4, z0=z0, record_from=11.0)
        a_, b_, h_ = rms_amplitudes(t, z, sec, d)
        A.append(a_); B.append(b_); H.append(h_)
        F.append(dominant_freq(t, z[:, 1]) if ns > 0 else np.nan)
        z0 = z[-1]                      # path-following: carry the state up
    out[name] = (A, B, H, F)
    for k, (v, lab) in enumerate(zip((A, B, H, F),
                                     (r"$A_{rms}$ pitch", r"$B_{rms}$ flap",
                                      r"$H_{rms}$ plunge", "frequency [Hz]"))):
        axi = ax[k//2, k % 2]
        axi.plot(ratios, v, ".-", ms=4, lw=.9, label=name)
        axi.set_xlabel("$U/U_f$"); axi.set_ylabel(lab); axi.grid(alpha=.3)

for axi in ax.ravel():
    for xv in (0.35, 0.44, 0.50):
        axi.axvline(xv, color="k", ls=":", lw=.8)
ax[0, 0].legend(fontsize=8)
ax[1, 1].set_ylim(0, 14)
plt.tight_layout(); plt.savefig("freeplay_sweep.png", dpi=140)
print(f"U_f = {UF:.3f} m/s;  wrote freeplay_sweep.png")
np.savez("freeplay_sweep.npz", ratios=ratios, UF=UF,
         **{f"{k}_{n}": v for k, vs in out.items()
            for n, v in zip("ABHF", vs)})
