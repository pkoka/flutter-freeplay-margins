## Abstract

A subcritical flutter test works like this. You excite the structure at several airspeeds below the
flutter boundary, identify the modal damping at each one, and extrapolate the damping to zero. The
speed where it reaches zero is your flutter speed. This memo asks what that procedure gives you when
the structure has control-surface freeplay, and what you can guarantee instead when it fails.

The test article is the three-degree-of-freedom typical section from Conner, Tang, Dowell & Virgin
(1997). I rebuilt it from Theodorsen's aerodynamic equations and checked it against published values
at every stage: the worked flutter example in Hodges & Pierce, Theodorsen's 1935 constant tables, and
the measured Duke rig. The model matches Conner's flutter speed to 0.29% and his coupled natural
frequencies to 0.13%.

On that model, the extrapolated flutter speed changes by a factor of 2.1 when the only thing that
changes is how hard you hit the structure. At large amplitude the test reports no flutter boundary at
all. The same variation shows up inside a single decay record. Nothing in the identification
indicates a problem.

There is still something you can guarantee. The freeplay restoring moment is sector-bounded, so the
plant splits exactly into a linear system with a nonlinearity in feedback. A Popov multiplier then
certifies a speed below which no disturbance of any size can produce sustained motion. That bound
needs no knowledge of the amplitude, the waveform, or the gap size.

### Principal results

| quantity | result |
|---|---|
| Extrapolated flutter speed, same rig and method | ranges over 0.259 to 0.544 $U_f$ with excitation amplitude alone; two amplitudes report no boundary |
| Same, within one decay record | 53 windows give ζ from −0.0041 to +0.0551 and frequency from 4.35 to 25.37 Hz |
| Identification on a limit-cycling structure | returns limit-cycle harmonics as modal parameters, and gives identical output for two different gap sizes |
| Rigorous amplitude-independent bound | Popov certifies stability below **0.1612 $U_f$**, against observed onsets of 0.18 and 0.20 $U_f$ |
| Model validation | flutter speed +0.29%, flutter frequency −0.11%, coupled frequencies within 0.13% |

---

## 1. The question

A subcritical flutter test assumes the structure is linear and time-invariant. That assumption is
what lets damping measured at 60% of the flutter speed tell you anything about behaviour at 95%.
Control-surface freeplay breaks it. The flap has no restoring stiffness inside the gap and full
stiffness outside it, so the effective stiffness the test sees depends on how hard you hit the
structure.

So the question is not whether identification is accurate. On a linear plant it is very accurate, as
Section 7 shows. The question is whether the thing it identifies is a property of the structure at
all. That splits into three parts:

1. Does the identified flutter margin depend on the test amplitude, and by how much?
2. Does anything in the identification tell you there is a problem?
3. If the point estimates cannot be trusted, what can you still guarantee?

Sections 4, 5 and 6 answer these in order. Section 3 comes first because none of those answers mean
anything unless the model is right.

## 2. Scope and contribution

**Reproduced from the literature.** This is here to show the tooling is correct before it is used to
claim anything: the Hodges & Pierce worked flutter example in both steady-flow and unsteady form,
Theodorsen's tabulated constants, the coupled natural frequencies and linear flutter point of the
Conner rig, and the freeplay limit-cycle amplitudes and regime boundaries from Conner's Figure 8.

**Built.** An eight-state model of the Conner three-degree-of-freedom section, assembled from
Theodorsen's aerodynamic equations with a two-lag Jones approximation. A piecewise-linear integrator
with Hénon switching and shaker forcing. A measurement model that matches the rig's actual
instrumentation. Implementations of ERA/OKID, PO-MOESP and weak-form second-order identification.

**Measured.** How much the extrapolated flutter margin depends on excitation amplitude, both across
records and within one record. How output-only identification behaves on a limit-cycling structure.
How the three identification methods compare against known ground truth across noise levels.

**Applied.** An exact sector decomposition of the freeplay plant, then circle, Popov, Zames–Falb and
combined-multiplier bounds on the speed below which no disturbance can produce sustained motion.

Absolute stability and describing-function methods are well established in aeroelasticity with
structural nonlinearity, and I am not claiming anything new against that literature. What this memo
adds is the quantified demonstration on one validated rig, plus the verification record in Section 8.

## 3. Model and validation

### 3.1 Two-degree-of-freedom section

I built a pitch-plunge typical section with Theodorsen unsteady aerodynamics and checked it against
the configuration in Hodges & Pierce Chapter 5 (a = −1/5, e = −1/10, µ = 20, r² = 6/25, σ = 2/5). The
values printed there use steady-flow aerodynamics.

| quantity | computed | printed |
|---|---|---|
| $V_F$ = $U_F/(b\omega_\theta)$ | 1.8425 | 1.843 |
| $\Omega_F/\omega_\theta$ | 0.5568 | 0.5568 |
| $V_D$, divergence | 2.8284 | 2.828 |

The divergence value also satisfies the closed form $V_D$ = r·√(µ/(1+2a)) exactly. That is a separate
check on a, µ and r².

With Theodorsen unsteady aerodynamics the same configuration gives $V_F$ = 2.1839,
$\Omega_F/\omega_\theta$ = 0.6490. That number is not comparable to 1.843. Same structure, different
aerodynamic model. The 18% gap between them is what neglecting unsteadiness costs.

Four solvers agree at the flutter point. The k-method and p-k with exact C(k) both give 2.1839. The
p-k method and the time-domain state space with the R. T. Jones two-lag Wagner approximation both
give 2.1704. Two of those agreements are guaranteed by theory rather than earned. At neutral
stability the motion is harmonic, so the harmonic assumption is exact, which makes them algebra
checks. The number that means something is the **0.62% Jones approximation error**. I also
cross-checked against Hodges' own aerodynamic formulation, which uses nondimensional coefficients
instead of dimensional lift and moment. It agrees to four decimals on two different parameter sets.

### 3.2 Theodorsen constants

The three-degree-of-freedom model needs Theodorsen's T-constants. I computed them in closed form and
compared them against Theodorsen's own tabulated values in NACA Report 496, Table I, at c = 0.5.

The comparison is against the last printed decimal place rather than as a percentage. Theodorsen's
table is printed to fixed decimals, so a percentage makes the smallest constants look worst. "ulp"
below means units in the last printed place.

| constant | computed | NACA 496 | ulp | reading |
|---|---|---|---|---|
| T₁ | −0.125920 | −0.1259 | 0.20 | rounding |
| T₄ | −0.614185 | −0.6142 | 0.15 | rounding |
| T₇ | +0.013250 | +0.0132 | 0.50 | rounding |
| T₁₀ | +1.913223 | +1.9132 | 0.23 | rounding |
| T₁₁ | +1.299038 | +1.2990 | 0.38 | rounding |
| T₅ | −0.939723 | −0.9398 | 0.77 | about 1 ulp |
| T₁₂ | +0.070668 | +0.07066 | 0.84 | about 1 ulp |
| **T₃** | −0.053203 | −0.05313 | **7.26** | outside rounding |
| **T₈** | +0.090586 | +0.0903 | **2.86** | outside rounding |

Seven of the nine agree to within one unit in the last place, which is normal hand-rounding. Two do
not. T₃ is off by 7.3 units and T₈ by 2.9 units. Rounding cannot account for either.

I am reporting this as an open minor discrepancy, not as a claim that the table is wrong. Perry's
recomputation of NACA 496 (TP-2015-218765) covers Table I and reports excellent agreement with the
original, and his list of identified errors in NACA 496 does not include Table I. But his comparison
is graphical. He plots the tabulated values as symbols against a computed curve, which cannot resolve
differences of order 10⁻⁴. So the two results do not actually contradict each other, and the
discrepancy in T₃ and T₈ is still unexplained. Neither constant affects any result here. T₈ does not
appear in the aerodynamic matrices at all, and substituting Theodorsen's tabulated T₃ moves the
computed flutter speed of the Conner rig by less than 10⁻⁵ percent. It is recorded because I found
it, not because it matters.

One obstacle is worth writing down. NASA/TM-2017-219667 recasts Theodorsen's equations into the
matrix form I used, and it contains the symbol *p* in three entries of its noncirculatory damping
matrix without ever defining it. The symbol is missing from its nomenclature. The definition is in
the original NACA 496, Appendix I: p = −(1/3)(1−c²)^(3/2). The identity
$T_9$ = ½[(1/3)(√(1−c²))³ + aT₄] = ½(−p + aT₄) confirms it, and the closed forms reproduce that
exactly. Perry notes separately that Table I lists only 11 of the 14 constants. The missing ones are
$T_9$, T₁₃ and T₁₄, which depend on the elastic axis location as well as the hinge location. That is
why they have to be reconstructed rather than looked up.

### 3.3 The Conner rig

I then rebuilt the plant as the Duke apparatus from Conner et al. (1997): NACA 0012, 0.254 m chord,
0.52 m span, elastic axis at the quarter chord, 25% chord flap.

**The parameter table checks against itself.** Appendix A of that paper gives both raw and derived
quantities, so they can be compared.

| quantity | recomputed | printed |
|---|---|---|
| $x_\alpha$ = $S_\alpha$/(Mb) | 0.43398 | 0.434 |
| $x_\beta$ = $S_\beta$/(Mb) | 0.01996 | 0.01996 |
| $r_\alpha$ = √($I_\alpha$/Mb²) | 0.73214 | 0.7321 |
| $r_\beta$ = √($I_\beta$/Mb²) | 0.11397 | 0.11397 |
| κ = πρb²/M | 0.03984 | 0.03984 |

The κ row pins ρ = 1.225 kg/m³ exactly, which the paper never states. Using ρ = 1.20 would give
0.03903 instead.

Table 1 of the paper is also over-determined. It prints flutter speed, flutter frequency and reduced
velocity for both the numerical and experimental cases, so the semichord can be solved for twice:
12.701 cm and 12.719 cm. The rig has a 10.000 inch chord exactly, and the paper's "19 cm main wing
chord" is a rounding of 7.5 in = 19.05 cm. This also fixes the reduced-velocity normalization as
$U/(b\omega_\alpha)$ and gives a = −0.500, c = +0.500.

**One line in that table decides the model.** Appendix A lists mass per unit length of wing-aileron
as 1.558 kg/m, and separately lists mass of support blocks as 0.47485 × 2 kg. The support blocks only
move in plunge, so their mass goes into the h-h entry of the mass matrix and nowhere else.

| | plunge | pitch | flap |
|---|---|---|---|
| wing + aileron only | 5.828 Hz | 11.671 | 19.658 |
| **+ support blocks in h-h** | **4.449** | **9.213** | **19.431** |
| Conner Table 1 | 4.455 | 9.218 | 19.442 |

That is 0.13%, 0.05% and 0.06%. Read the other way, the plunge frequency is wrong by 31%, and the
flutter speed comes out at 18.018 m/s against 23.970. **A 24.8% error from one ambiguous line in an
appendix.**

**Flutter validation.** Everything below includes Conner's measured structural damping, as his model
does. $U_f$ = 23.970 m/s is the normalizer used throughout this memo.

| method | $U_F$ | f_F |
|---|---|---|
| p-k, exact C(k) | 24.150 m/s (+1.05%) | 6.090 Hz (−0.36%) |
| **8-state Jones lag, eig A(U)** | **23.970 m/s (+0.29%)** | **6.105 Hz (−0.11%)** |
| Conner, numerical | 23.9 m/s | 6.112 Hz |

The 0.74% gap between the two solvers is the Jones fit error, which matches the 0.62% seen on the
two-degree-of-freedom section. Without structural damping the same two solvers give 23.915 m/s and
23.736 m/s, so how you model damping is worth about 1% on flutter speed.

A time-domain check confirms A(U) is a correct realization and not just a matrix with the right
eigenvalues. Integrating it and fitting the response envelope gives a growth rate of −1.5002 against
a predicted −1.5001 below the boundary, and +1.6888 against +1.6888 above it.

### 3.4 Freeplay

Conner's piecewise-linear model is three linear regions that share one A matrix. The free region
differs only in $K_\beta$, and the two stiff regions differ only in a constant offset vector. I
verified this directly by comparing the coefficient matrices: the only entries that differ are the
three acceleration rows in the β column. Switching points are found with Hénon's method, which
integrates in space onto the switching surface instead of interpolating or bisecting.

The model reproduces the structure of Conner's Figure 8. Zero response, then abrupt onset, then a
low-frequency regime near 5 Hz, then an abrupt jump near 0.50 $U_f$ to a high-frequency regime near
11.3 Hz where the plunge amplitude collapses, then growth toward the linear boundary. Amplitudes at
0.49 $U_f$ and 0.73 $U_f$ match Conner's Figures 10 and 11 across all three degrees of freedom and
over two orders of magnitude in plunge.

The three gap sizes collapse onto a single curve. The largest difference between the 1.15° and 2.12°
gaps across all normalized amplitudes is 2.3 × 10⁻¹². That collapse is Conner's analytical scaling
claim, and I did not impose it.

![**Figure 1.** Freeplay path-following at three gap sizes, with amplitudes normalized by gap size.
The three curves collapse onto one. The dotted line marks the jump to the high-frequency limit cycle
near 0.50 $U_f$.](fig3_freeplay_sweep.pdf){width=95%}

## 4. The identified margin depends on the test amplitude

### 4.1 Across records

Five subcritical speeds from 0.10 to 0.26 $U_f$, shaker raps of varying size, PO-MOESP
identification, then a linear extrapolation of damping to zero.

| max \|β\|/gap | extrapolated $U_F$ / true $U_F$ |
|---|---|
| 0.10 | 0.381 |
| 0.96 | 0.381 |
| 1.51 | 0.259 |
| 2.78 | 0.544 |
| 6.00 | ∞, damping rises with airspeed |
| 19.32 | ∞ |

The first two rows are bit-identical. Below one gap width the flap never leaves the deadband, so
those two tests are the same experiment. Once the flap reaches the stops the identified damping jumps
by a factor of five and the shape of the curve inverts.

Every finite prediction here is conservative, by 1.8× to 3.9×. But that is a property of this
configuration and I have not shown it holds in general. **A procedure whose answer moves by a factor
of 2.1 with a parameter nobody records is not usable for certification, whichever way it errs.**

### 4.2 Within a single record

A decaying record on a freeplay structure is not a stationary system. Early in the decay the flap
hits the stops. Late in the decay it rattles around inside the gap. I slid a 1.5 s window in 0.25 s
steps down one ring-down at 0.26 $U_f$, gap 2.12°, 40 dB SNR. That gives 53 estimates from one
record.

| quantity | range across windows |
|---|---|
| flap amplitude | 2.176 down to 1.056 gap widths |
| identified ζ | −0.0041 to +0.0551 |
| identified frequency | 4.35 to 25.37 Hz |

![**Figure 2.** One ring-down at $U/U_f$ = 0.26, gap 2.12°, 40 dB SNR. (a) is the flap response with
the gap edges marked. (b) and (c) are the identified damping and frequency from 53 overlapping
windows of that same record, coloured by window index. The record does not decay to zero. It settles
onto a limit cycle.](fig1_within_record.pdf){width=95%}

One structure, one airspeed, one record. Two engineers looking at different windows of the same data
get modal frequencies that differ by a factor of 5.8, and damping ratios that disagree on sign.

The record also does not decay to zero. It settles at |β|/gap = 1.072 with ζ near zero. That is a
limit cycle at 0.26 $U_f$, triggered by a large rap, below the fully-free instability speed of
0.286 $U_f$. The bifurcation is subcritical.

## 5. Identification does not reveal the problem

Run output-only identification on a self-excited limit cycle and what comes back is the harmonic
content of the cycle.

| gap | $U/U_f$ | identified |
|---|---|---|
| 1.83° | 0.73 | 6.89, 11.01, 12.94, 33.04 Hz |
| 2.12° | 0.73 | 6.89, 11.01, 12.94, 33.04 Hz, identical |

33.04/11.01 = 3.001, an exact third harmonic. The other lines are sidebands. Every damping ratio
satisfies |ζ| < 0.003, which is what a periodic signal looks like to a linear realization algorithm.
The truth at that speed is 5.14/8.63/19.39 Hz in the stiff region, or 5.50/6.64/8.31 Hz in the free
region. The identification matches neither, returns bit-identical answers for two different gap
sizes, and gives no warning.

This is the part that matters operationally. Subspace identification returns clean, well-conditioned
poles with plausible damping at every amplitude. Nothing in the residual, the pole locations, or the
extrapolation flags a problem. The failure is silent.

## 6. What can still be guaranteed

### 6.1 Exact sector decomposition

The only thing that differs between the stiff and free regions is the flap restoring moment. So the
plant splits exactly into a linear system with a memoryless nonlinearity in feedback:

> ż = $A_{\mathrm{free}}$·z + $B_\varphi$·(−φ(β)),  β = $C_\beta$·z
>
> $A_{\mathrm{stiff}}$ = $A_{\mathrm{free}}$ − $B_\varphi$·$K_\beta$·$C_\beta$  (verified to 3.6 × 10⁻¹²)

The freeplay moment satisfies 0 ≤ φ(β)/β ≤ $K_\beta$ for every β, because the deadband is flat in the
middle and has slope $K_\beta$ outside it. This is exact. No averaging, no assumed waveform, no
amplitude parameter. That is what puts absolute stability theory in play.

### 6.2 Bound hierarchy

| criterion | certified speed | status |
|---|---|---|
| circle | nothing at any speed | rigorous, but allows time-varying φ |
| **Popov** | **0.1612 $U_f$** | **rigorous** |
| describing-function minimum | 0.1612 $U_f$ | heuristic, minimum at γ = 0.050 |
| Zames–Falb alone | nothing | see below |
| Popov × Zames–Falb | 0.2371 $U_f$ | unresolved, see Section 9 |
| free-region limit | 0.2858 $U_f$ | hard ceiling, $A_{\mathrm{free}}$ loses stability |
| observed LCO, 2.12° / 1.15° | 0.18 / 0.20 $U_f$ | simulation |
| Conner, numerical / experimental | 0.18 / 0.23 $U_f$ | published |

![**Figure 3.** (a) Instability onset over the effective-stiffness family
$K_{\mathrm{eff}}$ = γ$K_\beta$, showing the interior minimum at γ = 0.050. (b) The bound hierarchy on
a single axis. The shaded region is certified by the Popov multiplier. The combined Popov × Zames–Falb
value is drawn as an open marker because it exceeds the observed limit-cycle onsets and is
unresolved.](fig2_bounds.pdf){width=88%}

**Below 0.1612 $U_f$, no disturbance of any size can produce sustained motion, for any nonlinearity
in the sector.** The bound needs no knowledge of amplitude, waveform or gap size. It sits below both
limit-cycle onsets I observed and below both onsets Conner et al. report.

**Why circle fails and Popov works.** The circle criterion allows the nonlinearity to be
time-varying. Freeplay is memoryless and time-invariant. Using that extra piece of information is the
entire difference between certifying nothing and certifying to within 10% of the observed onset.

**Why Zames–Falb alone fails.** Zames–Falb multipliers use monotonicity, which freeplay has. But with
‖h‖₁ ≤ 1 the multiplier satisfies |M(jω)| ≤ 2, so it can only supply bounded phase advance. The Popov
multiplier 1 + jωη is unbounded. On this loop the binding requirement turns out to be phase, not
monotonicity.

**The agreement with the describing function is not a theorem.** Popov and the describing-function
minimum agree to four digits on this configuration. Take the structural damping out and they diverge
completely. The describing-function minimum drops to 0.00209 $U_f$ while Popov certifies nothing
above about 0.005 $U_f$. The agreement is specific to this configuration.

## 7. Which identification method

All comparisons use the rig's own instrumentation: three RVDTs measuring pitch, flap and plunge
displacement, 800 Hz sampling, 4000-point records, and a shaker at midspan quarter chord with
measured force. Ground truth is eig A(U) of the plant that generated the data. The metric is mean
relative error in modal damping over the three modes. Modes are matched to truth by optimal
assignment on relative frequency, after discarding candidates with ζ > 0.5, which the true modes of
this plant never approach. Every method is scored the same way.

### 7.1 What model order is even admissible

The plant has 8 states and 3 measured outputs, so a model written in the measured outputs alone
cannot represent it. This can be made exact rather than qualitative. An order-d delay model must
satisfy X·O_d = C_d A_d^d, where O_d is the dl × n observability stack. Once O_d has rank n, the
solution set is affine with dimension l(dl − n), and that is the entire freedom any estimator has to
set coefficients to zero. Everything else is pinned by the plant.

| d | coefficients | freedom | minimum nonzero | max achievable sparsity |
|---|---|---|---|---|
| 2 | 18 | no exact model exists | | |
| 3 | 27 | 3 | 24 | 89% |
| 4 | 36 | 12 | 24 | 67% |
| 6 | 54 | 30 | 24 | 44% |

At d = 2 there is no exact model at all, with a residual of 1.3 × 10⁻³. That is the observability
index: 2 × 3 = 6 < 8. At d = 3 an L1 search over the full family of exact models finds the sparsest
one has exactly 24 of 27 nonzero coefficients, verified to a residual of 5.7 × 10⁻¹³. The three zeros
sit at the 10⁻¹⁶ level and the next smallest coefficient is 8.4 × 10⁻⁵.

Two consequences. Order 3 is the lowest order that can represent this plant, in delay coordinates or
as a matrix ODE. And sparse regression has nothing to find: the true representation is 89% dense at
the minimum viable order and no sparser in absolute terms at any higher order.

### 7.2 Results

| method | 60 dB | 40 dB | 30 dB | 20 dB | all |
|---|---|---|---|---|---|
| **ERA, free decay** | **0.04** | **0.08** | **0.16** | **0.61** | **0.22** |
| PO-MOESP order 8 | 0.04 | 0.60 | 0.72 | 2.17 | 0.88 |
| PO-MOESP order 6 | 0.58 | 0.49 | 0.65 | 3.52 | 1.31 |
| weak form, order 2 | 2.36 | 2.53 | 3.34 | 3.49 | 2.93 |
| strong form, order 2 | 4.42 | 3.91 | 1.34 | 16.09 | 6.44 |
| weak form, order 3 | 0.66 | 4.47 | 13.44 | 50.73 | 17.32 |
| strong form, order 3 | 50.88 | 43.17 | 71.37 | 39.04 | 51.12 |

**The weak form does what it claims, and the effect is largest where it matters.** At order 3 the
strong form has to differentiate noisy data three times, and it never recovers: 50.88% error even at
60 dB against 0.66% for the weak form, a factor of 77. At order 2 the gap is only 3.22% against
2.74% over 30 to 60 dB. So the weak form's advantage grows with model order, which is the regime the
plant actually requires.

**But at admissible order, regression loses to projection.** Weak form order 3 goes 0.66, 4.47,
13.44, 50.73 as noise increases. PO-MOESP order 8 goes 0.04, 0.60, 0.72, 2.17. ERA on free-decay
records, which is the actual flutter-test protocol, holds at 0.61% even at 20 dB. Subspace methods
use far more delays than order 3 and still win, because they project the output noise out instead of
regressing noisy outputs against noisy outputs.

The apparent tie between the weak form and subspace identification at 20 dB, 3.49% against 3.52%,
occurs only at order 2. Order 2 cannot represent this plant. At order 3 the same comparison is
50.73% against 2.17%.

**Observability.** With three displacement channels, all six structural state directions sit above a
40 dB noise floor at every test speed, with normalized Hankel singular values of 5 × 10⁻² or higher.
The two aerodynamic lag states sit one to two orders below that and cannot be identified separately.
Adding a fourth sensor does not change this.

## 8. Verification record

I found and fixed fourteen errors while doing this work. They are listed because the headline result is
that a standard procedure fails silently, so the discipline used to catch silent failures is part of
the result.

| # | error | how it was caught |
|---|---|---|
| 1 | Wagner's φ used for the gust input channel | φ(0) = 0.5 but Küssner's ψ(0) = 0, since a gust has to penetrate the chord first |
| 2 | "Dryden" turbulence was a first-order lag with the label attached | Compared against the MIL-F-8785C spectrum, which is second-order with a √3 zero |
| 3 | SNR setting delivered 44 dB when 30 was requested | Broadband noise added at the analog rate is mostly removed by the anti-alias filter |
| 4 | Two rational C(k) fits claimed as independent | Expanding them algebraically showed they are the same fit |
| 5 | Freeplay switching target taken from the sign of β after the step | Leaving the stiff region downward aimed at the opposite stop, giving 426° of flap rotation on a 25 cm chord |
| 6 | Weak-form test function had about 1 Hz bandwidth, modes at 5 to 19 Hz | Damping error rose steadily as the bandwidth fell below the top mode |
| 7 | ERA paired with OKID on swept-sine data | ERA put all its states on 180 to 290 Hz noise poles; on free decay it gives 0.14% |
| 8 | Popov multiplier search capped at η = 100 | The optimum sat at the search endpoint; the real optimum spans many decades |
| 9 | Frequency grid too coarse for the multiplier LP | The bound moved from 0.2432 to 0.2370 on refinement |
| 10 | Nearly compared an unsteady flutter speed against a steady-flow published value | The figure captions in Hodges & Pierce say steady-flow theory |
| 11 | Identification methods tuned on noise-free data, then scored under noise | Regularized methods collapse, and the protocol picks the wrong hyperparameters |
| 12 | Windowed identification used an unseeded dither | Two consecutive runs gave different ζ ranges |
| 13 | Modes matched to truth greedily by nearest frequency | A spurious pole at ζ = 0.985 was assigned to a real mode, scoring PO-MOESP at 787% in a run where it recovered every mode to four digits |
| 14 | Section 7 originally compared methods at order 2 | Order 2 admits no exact model of this plant, so the comparison was run in a model class that cannot represent it |

Items 4, 7, 10, 11, 13 and 14 would each have produced a published claim that is false rather than
just imprecise.

## 9. Open items

**The combined-multiplier bound is unresolved, so I am not offering it as a result.** Popov ×
Zames–Falb gives 0.2371 $U_f$, which is above the earliest limit cycles I observed at 0.18 and
0.20 $U_f$. A valid global stability certificate cannot be beaten by a real instability, so one of
the two is wrong. The most likely explanation is that those limit cycles sit on the subcritical
branch and are only reached by large finite disturbances, while absolute stability certifies global
asymptotic stability. But global stability and finite-disturbance limit cycles cannot both hold, so
that explanation needs checking rather than assuming. The other candidate is that I misapplied the
‖h‖₁ ≤ Σ|cᵢ| sufficiency condition in the odd-nonlinearity case. Only the Popov bound of
0.1612 $U_f$ is quoted above.

**The certificate is gridded, not analytic.** The multiplier search enforces its constraint on a
frequency grid, so it can miss violations between grid points. The value is converged, going 0.2432,
0.2371, 0.2371, 0.2370 over 600 to 6000 points. Real rigor would need a KYP/LMI formulation.

**The bound is valid but not sharp.** The sector argument is gap-independent by construction, but the
observed onset is not. It is 0.18 $U_f$ for the 2.12° gap against 0.20 $U_f$ for 1.15°. That leftover
dependence is physics the quasi-static sector argument throws away.

**Order and bandwidth selection are not costed.** I picked the model order knowing the true system
dimension, and the weak form needed retuning for each noise level. On real data you know neither.
That uncertainty is not included in any number above.

**Model fidelity.** Theodorsen incompressible two-dimensional aerodynamics with a two-lag Jones
approximation, which contributes about 0.75% error on flutter speed. Conner et al. attribute the 15%
gap between their numerical and experimental flutter speeds to aerodynamics outside the Theodorsen
model and to three-dimensional tunnel effects, with a measured 12.3% drift in pitch torsional
stiffness over the test campaign as part of the cause. This work reproduces their numerical model,
not their experiment.

**Not attempted.** Control design. The nonperiodic regime between 0.35 and 0.44 $U_f$, which needs a
Poincaré map rather than RMS amplitude. Any experimental data of my own.

## Appendix: reproducibility

| file | contents |
|---|---|
| `typical_section.py` | 2-DOF plant; k-method, p-k, steady-flow p-method, Jones state space |
| `hodges_kmethod.py` | Independent port of the Hodges aerodynamic formulation |
| `conner3dof.py` | 3-DOF Conner plant; Theodorsen T-constants, p-k, 8-state Jones realization |
| `freeplay.py` | Piecewise-linear integration with Hénon switching and shaker input |
| `phase1_conner.py` | Measurement model: three RVDTs, 800 Hz, shaker, noise, anti-alias, ADC |
| `sysid.py` | ERA, OKID, PO-MOESP, free-decay ERA |
| `weakform.py` | Weak-form second-order identification |
| `margin_amplitude.py` | Amplitude dependence across records |
| `margin_envelope.py` | Windowed identification within one record; describing-function family |
| `absolute_stability.py` | Sector decomposition; circle and Popov criteria |
| `zames_falb.py` | Zames–Falb and combined multiplier search |
| `make_figures.py` | Generates all three figures; every plotted value recomputed at run time |

Every validation number in Section 3 can be reproduced by running these files. `make_figures.py`
recomputes each plotted quantity from the model instead of reading stored results, and its output is
reproducible across runs.

## References

Conner, M. D., Tang, D. M., Dowell, E. H. & Virgin, L. N. (1997). Nonlinear behavior of a typical
airfoil section with control surface freeplay: a numerical and experimental study. *Journal of
Fluids and Structures* **11**, 89–109.

Edwards, J. W., Ashley, H. & Breakwell, J. V. (1979). Unsteady aerodynamic modeling for arbitrary
motions. *AIAA Journal* **17**(4), 365–374. (Cited by Conner et al. as the source of the state-space
submatrix forms. The realization used here was assembled independently from Theodorsen (1935) via
Perry (2017), so this reference is given for provenance rather than as a source used directly.)

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
