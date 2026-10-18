# Flutter margin prediction on a structure with control-surface freeplay

Code and data for the memo *Limits of Identification-Based Flutter Margin Prediction on a Structure
with Control-Surface Freeplay* (`flutter_margin_memo.pdf`).

A subcritical flutter test excites a structure at several airspeeds below the flutter boundary,
identifies the modal damping at each one, and extrapolates that damping to zero. This repository
asks what the procedure gives you when the structure has control-surface freeplay, and what you can
guarantee instead when it fails.

The test article is the three-degree-of-freedom typical section from Conner, Tang, Dowell & Virgin
(1997). It was rebuilt from Theodorsen's aerodynamic equations and checked against published values
at every stage.

## Main results

| | result |
|---|---|
| Model validation | flutter speed +0.29%, flutter frequency −0.11%, coupled natural frequencies within 0.13% of Conner et al. |
| Extrapolated flutter speed, same rig and method | ranges over 0.258 to 0.454 U_f with excitation amplitude alone. Two amplitudes report no boundary at all |
| Same, within one decay record | 50 windows give ζ from −0.0075 to +0.0824 and frequency from 4.52 to 26.09 Hz |
| Identification on a limit-cycling structure | returns limit-cycle harmonics as modal parameters. Identical output for two different gap sizes |
| Rigorous amplitude-independent bound | Popov certifies stability below **0.1612 U_f**, against observed limit-cycle onsets of 0.18 and 0.20 U_f |
| Same, about a preloaded operating point (October) | rises to 0.2079, 0.3139, 0.5590, 0.9681 U_f at hinge preloads of 0.25, 0.5, 1, 2 gap widths; no simulated trajectory beats it |
| Closed-loop back-calculation (October) | plant recovered to 2% (pitch) and 0.4% (flap) through and past flutter; flutter speed 3.9% conservative, against 1.7% unconservative for open-loop extrapolation |

The October revision (memo Sections 6–8, 10) states the operating point of the certificate, adds the
preload study, measures model error and test scatter approaching onset, and identifies the plant in
closed loop. It also corrects a sample-rate error and an invalid combined multiplier from August.

## Setup

Requires Python 3.10 or later.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run order

Everything runs from the repository root. The modules import each other by bare name, and
`make_figures.py` writes the figure PDFs into the working directory, so do not reorganize into
subfolders without adding path handling.

**Check the install first.**

```bash
python -c "
import numpy as np
from conner3dof import ConnerSection, ss_flutter
UF, fF = ss_flutter(ConnerSection(), np.linspace(1., 40., 3000), True)
print(f'U_f = {UF:.3f} m/s, f_f = {fF:.3f} Hz')"
```

This should print `U_f = 23.970 m/s, f_f = 6.105 Hz`. If it does not, stop.

**Then, in this order:**

| script | time | what it does |
|---|---|---|
| `phase0_validate.py` | 10 s | 2-DOF section validated against Hodges & Pierce. Four solvers, plus the steady-flow case that matches the printed values. |
| `hodges_kmethod.py` | 5 s | Independent port of Hodges' own aerodynamic formulation. Cross-check on the above. |
| `phase1_conner_check.py` | 2 min | Measurement model: three RVDTs, 800 Hz, shaker. Writes the linear and freeplay datasets. |
| `freeplay_sweep.py` | 10 min | Freeplay path-following at three gap sizes. Reproduces Conner Figure 8. Writes `freeplay_sweep.npz`. |
| `section7_rerun.py` | 8 min | Identification comparison: ERA, PO-MOESP, weak form and strong form at orders 2 and 3. |
| `margin_amplitude.py` | 6 min | Extrapolated flutter speed against excitation amplitude. |
| `margin_envelope.py` | 5 min | Windowed identification inside one record, plus the describing-function family. |
| `absolute_stability.py` | 3 min | Sector decomposition, circle and Popov criteria. |
| `zames_falb.py` | 10 min | Zames-Falb and Popov + Zames-Falb search. |
| `preload.py` | 35 min | Preload: operating point, narrowed-sector Popov, describing function, 3,968 simulations. Writes `preload_results.npz`. |
| `item1_residual_sparsity.py` | 1 min | The sparse-recovery negative result. |
| `stage2_uncertainty.py` | 3 min | Measured model error, robust Popov. Writes a `grid.pkl` cache. |
| `stage2_hinf.py` | 8 min | H∞ designs and the robust-stability test against the measured W. |
| `stage3_scatter.py` | 10 min | Replicate scatter and measurement quality approaching onset. |
| `closed_loop_id.py` | 5 min | Flutter suppression loop and closed-loop back-calculation. |
| `make_figures.py` | 20 min | Regenerates all four figures. Requires `freeplay_sweep.npz` and `preload_results.npz`. |

Times are for a 15-core laptop; the scripts that run many simulations use 12 worker processes.

`make_figures.py` prints every value it plots. Those should read Popov 0.1612, describing-function
minimum 0.1612 at γ = 0.050, Popov + Zames-Falb 0.1612, free-region ceiling 0.2858. Those are the
numbers in memo Section 6.2.

## What each file does

**Plant models**

- `typical_section.py`: 2-DOF pitch-plunge section. k-method, p-k, steady-flow p-method, and a Jones
  two-lag state space.
- `hodges_kmethod.py`: independent implementation of the Hodges aerodynamic formulation.
- `conner3dof.py`: 3-DOF Conner plant. Theodorsen T-constants, p-k solver, 8-state Jones realization.
- `freeplay.py`: piecewise-linear integration with Hénon switching and shaker forcing.

**Measurement**

- `phase1_conner.py`: three RVDTs, 800 Hz sampling, shaker at midspan quarter chord, sensor noise,
  anti-alias filtering and ADC quantization.

**Identification**

- `sysid.py`: ERA, OKID, PO-MOESP, and ERA on free-decay records.
- `weakform.py`: weak-form second-order fitting.
- `weakform_p.py`: weak-form and strong-form fitting at arbitrary order.
- `scoring.py`: mode matching by optimal assignment, and ground truth from `eig A(U)`.

**Results**

- `margin_amplitude.py`: amplitude dependence of the extrapolated margin, across records.
- `margin_envelope.py`: windowed identification within one record, and the describing-function
  envelope.
- `absolute_stability.py`: exact sector decomposition, circle and Popov criteria.
- `zames_falb.py`: Zames-Falb multipliers and the Popov + Zames-Falb search, odd and h ≥ 0 classes.
- `preload.py`: operating point under preload, narrowed-sector Popov, dual-input describing function,
  random-start simulation (memo Section 6.3).
- `section7_rerun.py`: the identification comparison in memo Section 9.

**Measurement near onset and closed-loop identification** (memo Sections 7 and 8)

- `stage2a_harness.py`: deliberately perturbed plant and periodic-multisine test harness.
- `stage2_uncertainty.py`: measured W, robust Popov, antiresonance floor and fitted weight.
- `hinf.py`: self-contained two-Riccati H∞ synthesis.
- `stage2_hinf.py`: H∞ designs, the robust-stability test against the measured W, and the
  freeplay check of the design that passes it.
- `stage3_scatter.py`: replicate scatter and how measurement quality degrades approaching onset.
- `closed_loop_id.py`: digital flutter suppression on the perturbed plant, and four ways of
  back-calculating the open-loop plant from closed-loop sweeps.
- `item1_residual_sparsity.py`: the sparse-recovery negative result.

**Figures**

- `make_figures.py`: regenerates all four figures. Every plotted value is recomputed from the model
  at run time rather than read from stored results, except the simulated preload amplitudes in
  figure 4, which come from `preload.py`.

## Reproducibility

Every number in memo Section 3 is reproduced by running these files. `make_figures.py` recomputes
each plotted quantity from the model instead of reading stored results, and its output is
reproducible across runs.

`.npz` and `.png` outputs are gitignored because every one of them is regenerated by the scripts.
The three figure PDFs are tracked, because the memo will not compile without them.

Note that matplotlib writes a creation timestamp into PDF metadata, so regenerating the figures can
show them as modified in git even when the plotted content is byte-identical.

## Rebuilding the memo

The committed `flutter_margin_memo.pdf` was built from `flutter_margin_memo.md` and
`memo_meta.yaml`. Rebuilding needs pandoc and a TeX installation:

```bash
brew install pandoc
brew install --cask basictex          # then restart your shell
sudo tlmgr update --self && sudo tlmgr install lmodern sectsty

pandoc memo_meta.yaml flutter_margin_memo.md -o flutter_margin_memo.pdf \
  --pdf-engine=xelatex --toc --toc-depth=2 -V toc-title="Contents"
```

## Known limitations

Read Section 11 of the memo before using any of this. The short version:

- The multiplier searches enforce their constraints on a frequency grid, so they can miss violations
  between grid points. Real rigor needs a KYP/LMI formulation.
- The describing function misses limit-cycle amplitude by up to about 20% under aerodynamic preload.
- Closed-loop flutter-mode damping is biased near the flutter speed, by up to 0.039 at 1.05 U_f.
- The "true" plant in memo Sections 7 and 8 is the model with 4% parameter scatter, a stand-in for
  real model error.
- Model order was chosen knowing the true system dimension, and the weak form needed retuning per
  noise level. On real data neither is known.
- The aerodynamics are Theodorsen incompressible two-dimensional with a two-lag Jones approximation,
  contributing about 0.75% error on flutter speed.
- This work reproduces Conner's numerical model, not their experiment.

Memo Section 10 lists twenty-seven errors found and corrected during this work, along with how each
was caught.

## License

MIT, see `LICENSE`. This covers the code. The memo text and figures are the author's own work and
may be cited normally.

## References

Conner, M. D., Tang, D. M., Dowell, E. H. & Virgin, L. N. (1997). Nonlinear behavior of a typical
airfoil section with control surface freeplay: a numerical and experimental study. *Journal of Fluids
and Structures* **11**, 89–109.

Edwards, J. W., Ashley, H. & Breakwell, J. V. (1979). Unsteady aerodynamic modeling for arbitrary
motions. *AIAA Journal* **17**(4), 365–374.

Hodges, D. H. & Pierce, G. A. (2002). *Introduction to Structural Dynamics and Aeroelasticity*.
Cambridge University Press.

Jones, R. T. (1938). Operational treatment of the nonuniform lift theory in airplane dynamics.
NACA TN 667.

Perry, B. III (2015). Recomputation of numerical results contained in NACA Report No. 496.
NASA/TP-2015-218765.

Perry, B. III (2017). Comparison of Theodorsen's unsteady aerodynamic forces with doublet lattice
generalized aerodynamic forces. NASA/TM-2017-219667.

Theodorsen, T. (1935). General theory of aerodynamic instability and the mechanism of flutter.
NACA Report 496.

## Author

Pratham Koka
