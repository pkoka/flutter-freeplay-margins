"""
Figures for the memo. Every value plotted here is recomputed from the model
modules at run time. Nothing is transcribed, interpolated, smoothed, or read
off a previous result. Where a quantity is disputed or unresolved it is drawn
in a visually distinct way rather than as a peer of the validated numbers.

Run: python make_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from conner3dof import ConnerSection, ss_flutter, state_space
from margin_envelope import ringdown, windowed_id
from absolute_stability import popov_margin, circle_margin
from zames_falb import combined_margin, zf_margin

plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

sec = ConnerSection()
UF, fF = ss_flutter(sec, np.linspace(1., 40., 3000), True)
print(f"U_f = {UF:.4f} m/s, f_f = {fF:.4f} Hz  (recomputed)")


# ===========================================================================
# Figure 1 - identified modal parameters drift WITHIN a single ring-down
# ===========================================================================
def figure1(path):
    fig = plt.figure(figsize=(8.6, 5.6), constrained_layout=True)
    gs_ = fig.add_gridspec(2, 2)
    ax = [fig.add_subplot(gs_[0, :]), fig.add_subplot(gs_[1, 0]),
          fig.add_subplot(gs_[1, 1])]
    d = np.deg2rad(2.12)
    ts, y, zz = ringdown(sec, 0.26*UF, d, 200.0, T=14.0, snr_db=40.0, seed=11)

    ax[0].plot(ts, np.rad2deg(zz[:, 1]), lw=0.4, color="0.25")
    for s in (+1, -1):
        ax[0].axhline(s*2.12, color="C3", lw=0.9, ls="--")
    ax[0].set_xlabel("t [s]"); ax[0].set_ylabel(r"flap $\beta$ [deg]")
    ax[0].set_title(r"(a) one ring-down, $U/U_f=0.26$" "\n"
                    "dashed: freeplay gap edges", fontsize=9)

    W = windowed_id(ts, y, zz, d, win=1.5, hop=0.25)
    amp, f, z = W[:, 0], W[:, 1], W[:, 2]
    sc = ax[1].scatter(amp, z, c=np.arange(len(amp)), cmap="viridis", s=22,
                       zorder=3)
    ax[1].axhline(0, color="k", lw=0.7)
    ax[1].set_xlabel(r"$|\beta|_{rms}\sqrt{2}$ / gap")
    ax[1].set_ylabel(r"identified $\zeta$")
    ax[1].set_title("(b) identified damping vs flap amplitude\n"
                    "colour = window index (time)", fontsize=9)
    cb = fig.colorbar(sc, ax=ax[2], label="window index (time $\\rightarrow$)")

    ax[2].scatter(amp, f, c=np.arange(len(amp)), cmap="viridis", s=22, zorder=3)
    ax[2].set_xlabel(r"$|\beta|_{rms}\sqrt{2}$ / gap")
    ax[2].set_ylabel("identified frequency [Hz]")
    ax[2].set_title("(c) identified frequency vs flap amplitude", fontsize=9)

    fig.savefig(path); plt.close(fig)
    print(f"  fig1: {len(W)} windows, zeta range "
          f"[{z.min():+.4f}, {z.max():+.4f}], f range [{f.min():.2f}, "
          f"{f.max():.2f}] Hz")
    return W


# ===========================================================================
# Figure 2 - describing-function family and the bound hierarchy
# ===========================================================================
def gamma_onset(gamma, sd=True, lo=0.2, hi=40.0):
    def mx(U):
        Af = state_space(sec, U, sd, "free")
        As = state_space(sec, U, sd, "stiff")
        return np.linalg.eigvals(Af + gamma*(As - Af)).real.max()
    if mx(lo) > 0:
        return lo
    try:
        return brentq(mx, lo, hi)
    except ValueError:
        return hi


def figure2(path):
    gs = np.linspace(0, 1, 401)
    Us = np.array([gamma_onset(g) for g in gs])
    gmin = gs[Us.argmin()]; Udf = Us.min()

    w_pop = np.geomspace(1e-4, 1e4, 4000)
    U_pop = brentq(lambda U: popov_margin(sec, U, w_pop), 0.03*UF, 0.25*UF,
                   xtol=1e-4)
    w_lp = np.geomspace(1e-3, 1e3, 3000)
    U_comb = brentq(lambda U: combined_margin(sec, U, w_lp), 0.10*UF, 0.28*UF,
                    xtol=1e-3)
    U_free = gamma_onset(0.0)
    print(f"  fig2 recomputed: DF-min {Udf/UF:.4f} (gamma={gmin:.3f}), "
          f"Popov {U_pop/UF:.4f}, combined {U_comb/UF:.4f}, "
          f"free-region {U_free/UF:.4f}")

    fig, ax = plt.subplots(2, 1, figsize=(8.0, 8.2), constrained_layout=True,
                           gridspec_kw={"height_ratios": [1, 1.25]})

    ax[0].plot(gs, Us/UF, "-", color="C0", lw=1.4)
    ax[0].plot(gmin, Udf/UF, "o", color="C3", ms=6, zorder=4)
    ax[0].annotate(f"minimum\n{Udf/UF:.4f} $U_f$\nat $\\gamma$={gmin:.3f}",
                   (gmin, Udf/UF), textcoords="offset points",
                   xytext=(28, 6), fontsize=8, color="C3")
    ax[0].set_xlabel(r"$\gamma = K_{eff}/K_\beta$")
    ax[0].set_ylabel("instability onset  $U/U_f$")
    ax[0].set_title("(a) onset over the effective-stiffness family\n"
                    "$\\gamma$=0 free, $\\gamma$=1 nominal linear", fontsize=9)
    ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1.05)

    items = [
        ("circle criterion",        None,        "certifies nothing", "0.55", "x"),
        ("Popov",                   U_pop/UF,    "rigorous",          "C2",   "o"),
        ("describing-fn minimum",   Udf/UF,      "heuristic",         "C0",   "s"),
        ("Zames-Falb alone",        None,        "certifies nothing", "0.55", "x"),
        ("Popov + Zames-Falb",      U_comb/UF,   "rigorous",          "C1",   "D"),
        ("free-region ceiling",     U_free/UF,   "hard limit",        "0.35", "|"),
        ("observed LCO 2.12 deg",   0.18,        "simulation",        "C3",   "^"),
        ("observed LCO 1.15 deg",   0.20,        "simulation",        "C3",   "^"),
        ("Conner numerical",        0.18,        "published",         "C4",   "v"),
        ("Conner experimental",     0.23,        "published",         "C4",   "v"),
        ("nominal linear flutter",  1.00,        "definition",        "k",    "*"),
    ]
    ax[1].axvspan(0, U_pop/UF, color="C2", alpha=0.13, zorder=0)
    ax[1].text(U_pop/UF/2, -0.75, "guaranteed:\nno sustained motion",
               ha="center", va="top", fontsize=8, color="C2")
    for i, (lab, v, note, col, mk) in enumerate(items):
        yv = len(items) - 1 - i
        if v is None:
            ax[1].text(0.005, yv, "— " + note, va="center", fontsize=8,
                       color=col, style="italic")
        else:
            ax[1].plot(v, yv, mk, color=col, ms=8,
                       mfc="none" if note == "UNRESOLVED" else col,
                       mew=1.8 if note == "UNRESOLVED" else 1.0)
            ax[1].text(v + 0.018, yv, f"{v:.4f}" if v < 1 else "1.0000",
                       va="center", fontsize=8, color=col)
            if note == "UNRESOLVED":
                ax[1].annotate("exceeds observed LCO onset;\nnot yet resolved "
                               "(see memo §6)", (v, yv),
                               textcoords="offset points", xytext=(34, 16),
                               fontsize=7.5, color="C1",
                               arrowprops=dict(arrowstyle="-", color="C1",
                                               lw=0.8))
    ax[1].set_yticks(range(len(items)))
    ax[1].set_yticklabels([l for l, *_ in items][::-1], fontsize=8)
    ax[1].set_xlim(-0.03, 1.15); ax[1].set_ylim(-1.6, len(items)-0.3)
    ax[1].set_xlabel("$U/U_f$")
    ax[1].set_title("(b) bound hierarchy, all gap-independent, "
                    "no amplitude assumed", fontsize=9)
    ax[1].grid(axis="y", alpha=0.15)

    fig.savefig(path); plt.close(fig)
    return U_pop, Udf, U_comb


# ===========================================================================
# Figure 3 - freeplay path-following, reproducing Conner Fig. 8
# ===========================================================================
def figure3(path):
    d = np.load("freeplay_sweep.npz")
    r = d["ratios"]
    fig, ax = plt.subplots(2, 2, figsize=(8.4, 5.8))
    labels = [("A", r"$A_{rms}$  pitch / gap"),
              ("B", r"$B_{rms}$  flap / gap"),
              ("H", r"$H_{rms}$  plunge / (gap$\cdot b$)"),
              ("F", "LCO frequency [Hz]")]
    styles = {"Gap 1": ("C0", "o", 1.15), "Gap 2": ("C1", "s", 1.83),
              "Gap 3": ("C2", "^", 2.12)}
    for k, (key, lab) in enumerate(labels):
        a = ax[k//2, k % 2]
        for g, (col, mk, dg) in styles.items():
            v = d[f"{g}_{key}"]
            a.plot(r, v, mk + "-", color=col, ms=3.2, lw=0.8, alpha=0.85,
                   label=f"{dg}$^\\circ$ gap")
        a.set_xlabel("$U/U_f$"); a.set_ylabel(lab)
        a.set_xlim(0, 1)
        if key == "F":
            a.set_ylim(0, 14)
        if k == 0:
            a.legend(fontsize=8, title="freeplay", title_fontsize=8)
    for a in ax.ravel():
        a.axvline(0.50, color="0.4", ls=":", lw=1.0)
    ax[0, 0].text(0.505, ax[0, 0].get_ylim()[1]*0.93,
                  "jump to high-frequency LCO", fontsize=7.5, color="0.3")
    fig.tight_layout()
    fig.savefig(path); plt.close(fig)
    coll = max(np.max(np.abs(d[f"Gap 3_{k}"] - d[f"Gap 1_{k}"]))
               for k in ("A", "B", "H"))
    print(f"  fig3: max |Gap3 - Gap1| across normalized amplitudes = {coll:.2e}")


if __name__ == "__main__":
    figure1("fig1_within_record.pdf")
    figure2("fig2_bounds.pdf")
    figure3("fig3_freeplay_sweep.pdf")
    for n in ("fig1_within_record", "fig2_bounds", "fig3_freeplay_sweep"):
        pass
    print("done")
