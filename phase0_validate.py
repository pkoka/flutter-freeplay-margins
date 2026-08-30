"""
Phase 0 cross-validation. Run this, look at the three panels, then move on.

Panel 1  V-g and V-f from the k-method (the classical flutter plot)
Panel 2  damping ratio vs U from p-k and from eig(A(U)) overlaid.
         These MUST coincide at the flutter point and MAY differ away from it.
         That is the whole point of the comparison.
Panel 3  root locus of A(U) in the complex plane
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from typical_section import (TypicalSection, theodorsen, theodorsen_jones,
                             k_method, pk_method, state_space,
                             flutter_k, flutter_pk, flutter_ss)

ts = TypicalSection(a=-0.2, e=-0.1, mu=20.0, r2=6/25, sigma=0.4)
b, wt = ts.b, ts.omega_th

k_vals = np.geomspace(3.0, 0.01, 4000)
U_vals = np.linspace(0.05, 4.0, 1500)

Uk, wk = flutter_k(ts, k_vals, theodorsen)
Up, wp = flutter_pk(ts, U_vals, theodorsen)
Upj, wpj = flutter_pk(ts, U_vals, theodorsen_jones)
Us, ws = flutter_ss(ts, U_vals)

print(f"{ts}\n  x_theta = {ts.x_th:.4f}")
print("\n  method                              V* = U/(b w_th)   w_f/w_th")
print(f"  k-method,  exact C(k)                    {Uk/(b*wt):8.4f}   {wk/wt:8.4f}")
print(f"  p-k,       exact C(k)                    {Up/(b*wt):8.4f}   {wp/wt:8.4f}")
print(f"  p-k,       Jones C(k)                    {Upj/(b*wt):8.4f}   {wpj/wt:8.4f}")
print(f"  state space, Jones lags, eig(A(U))       {Us/(b*wt):8.4f}   {ws/wt:8.4f}")
print(f"\n  Jones fit error       : {100*abs(Upj-Up)/Up:6.3f} %")
print(f"  harmonic assumption   : {100*abs(Us-Upj)/Upj:6.3f} %")
print(f"  k-method vs p-k exact : {100*abs(Uk-Up)/Up:6.3f} %")

# ---- data for the plots
Ug, omg, gg = k_method(ts, k_vals, theodorsen)
om_pk, z_pk = pk_method(ts, U_vals, theodorsen_jones)

z_ss = np.full((len(U_vals), 2), np.nan)
w_ss = np.full((len(U_vals), 2), np.nan)
locus = []
for i, U in enumerate(U_vals):
    ev = np.linalg.eigvals(state_space(ts, U))
    locus.append(ev)
    osc = ev[ev.imag > 1e-9]
    osc = osc[np.argsort(osc.imag)]
    for j in range(min(2, len(osc))):
        w_ss[i, j] = abs(osc[j])
        z_ss[i, j] = -osc[j].real/abs(osc[j])
locus = np.array(locus)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

m = (Ug[:, 0] > 0) & (Ug[:, 0] < 4)
for j, c in zip((0, 1), ("C0", "C3")):
    mj = (Ug[:, j] > 0) & (Ug[:, j] < 4)
    ax[0].plot(Ug[mj, j]/(b*wt), gg[mj, j], c, label=f"branch {j+1}")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].axvline(Uk/(b*wt), color="g", ls="--", lw=1,
              label=f"$V^*_F$ = {Uk/(b*wt):.3f}")
ax[0].set_xlim(0, 4); ax[0].set_ylim(-0.6, 0.3)
ax[0].set_xlabel(r"$U/(b\omega_\theta)$"); ax[0].set_ylabel("$g$")
ax[0].set_title("k-method  V-g,  exact C(k)"); ax[0].legend(); ax[0].grid(alpha=.3)

for j, c in zip((0, 1), ("C0", "C3")):
    ax[1].plot(U_vals/(b*wt), z_pk[:, j], c, lw=2, alpha=.5,
               label=f"p-k branch {j+1}")
    ax[1].plot(U_vals/(b*wt), z_ss[:, j], c, ls="--", lw=1.2,
               label=f"eig(A) branch {j+1}")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlim(0, 4); ax[1].set_ylim(-0.15, 0.35)
ax[1].set_xlabel(r"$U/(b\omega_\theta)$"); ax[1].set_ylabel(r"$\zeta$")
ax[1].set_title("p-k vs state space, both Jones C(k)")
ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)

sc = ax[2].scatter(locus.real.ravel(), locus.imag.ravel(),
                   c=np.repeat(U_vals/(b*wt), 6), s=3, cmap="viridis")
ax[2].axvline(0, color="k", lw=0.6)
ax[2].set_xlim(-1.2, 0.4); ax[2].set_ylim(-1.6, 1.6)
ax[2].set_xlabel(r"Re$(s)/\omega_\theta$"); ax[2].set_ylabel(r"Im$(s)/\omega_\theta$")
ax[2].set_title("root locus of $A(U)$"); ax[2].grid(alpha=.3)
plt.colorbar(sc, ax=ax[2], label=r"$U/(b\omega_\theta)$")

plt.tight_layout()
plt.savefig("phase0_validation.png", dpi=140)
print("\n  wrote phase0_validation.png")


# ---------------------------------------------------------------------------
# check 5: does A(U) actually integrate the way its eigenvalues say it should
# ---------------------------------------------------------------------------
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

print("\n  time-domain check (impulse in plunge, peak-envelope fit):")
for Ustar in (0.9*Us, 1.06*Us):
    A = state_space(ts, Ustar)
    sol = solve_ivp(lambda t, z: A @ z, [0, 120], np.array([0.01, 0, 0, 0, 0, 0]),
                    t_eval=np.linspace(0, 120, 40000), rtol=1e-12, atol=1e-14)
    pk, _ = find_peaks(np.abs(sol.y[0]))
    tp, hp = sol.t[pk], np.abs(sol.y[0][pk])
    m = (tp > 40) & (hp > 1e-11)
    rate = np.polyfit(tp[m], np.log(hp[m]), 1)[0]
    print(f"    U* = {Ustar/(b*wt):.3f}   sim rate = {rate:+.5f}   "
          f"max Re(eig) = {np.linalg.eigvals(A).real.max():+.5f}")
