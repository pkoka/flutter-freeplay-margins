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

On that model, the extrapolated flutter speed changes by a factor of 1.8 when the only thing that
changes is how hard you hit the structure. At large amplitude the test reports no flutter boundary at
all. The same variation shows up inside a single decay record. Nothing in the identification
indicates a problem.

There is still something you can guarantee. The freeplay restoring moment is sector-bounded, so the
plant splits exactly into a linear system with a nonlinearity in feedback. A Popov multiplier then
certifies a speed below which no disturbance of any size can produce sustained motion. That bound
needs no knowledge of the amplitude, the waveform, or the gap size.

**October revision.** The August version certified stability of the origin of a symmetric deadband.
A control surface under load is pushed against one side of its gap, so the origin is not an
equilibrium and that certificate did not say what it was about. Section 6 now states the operating
point. Doing so does more than fix the statement: a flap resting on a stop sees a narrower sector, and
the certified speed rises from 0.16 to 0.97 $U_f$ as the hinge preload grows. Section 7 measures how
model error and test scatter behave approaching onset, and Section 8 identifies the plant through and
past the flutter speed by back-calculating it from closed-loop sweeps, which is how the field answers
the open-loop identification problem. Two August results were wrong and are corrected, a sample-rate
error and an invalid combined multiplier, and so was a claim made after August that measurement
quality falls off a cliff near onset (Section 10).

### Principal results

| quantity | result |
|---|---|
| Extrapolated flutter speed, same rig and method | ranges over 0.258 to 0.454 $U_f$ with excitation amplitude alone; two amplitudes report no boundary |
| Same, within one decay record | 50 windows give ζ from −0.0075 to +0.0824 and frequency from 4.52 to 26.09 Hz |
| Identification on a limit-cycling structure | returns limit-cycle harmonics as modal parameters, and gives identical output for two different gap sizes |
| Rigorous amplitude-independent bound | Popov certifies stability below **0.1612 $U_f$**, against observed onsets of 0.18 and 0.20 $U_f$ |
| Same, about a preloaded operating point | certified speed rises to 0.2079, 0.3139, 0.5590 and 0.9681 $U_f$ at hinge preloads of 0.25, 0.5, 1 and 2 gap widths; no simulated trajectory beats it |
| Describing function under preload | onset matches simulation where the flap floats in the gap; amplitude 4–20% low; preload above half a gap removes the limit cycle in both |
| Closed-loop back-calculation | recovers the plant to 2% (pitch) and 0.4% (flap) through and past flutter, where open-loop testing is impossible; flutter speed 3.9% conservative |
| Indirect method with the documented controller | 34–52% plant error from a 5% gain error and an unmodelled filter, amplified by the loop gain |
| Model validation | flutter speed +0.29%, flutter frequency −0.11%, coupled frequencies within 0.13% |

---

## 1. The question

A subcritical flutter test assumes the structure is linear and time-invariant. That assumption is
what lets damping measured at 60% of the flutter speed tell you anything about behaviour at 95%.
Control-surface freeplay breaks it. The flap has no restoring stiffness inside the gap and full
stiffness outside it, so the effective stiffness the test sees depends on how hard you hit the
structure.

So the question is not whether identification is accurate. On a linear plant it is very accurate, as
Section 9 shows. The question is whether the thing it identifies is a property of the structure at
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
combined-multiplier bounds on the speed below which no disturbance can produce sustained motion,
about the origin and, since October, about a preloaded operating point. A dual-input describing
function checked against simulation. A digital flutter-suppression loop on a deliberately mismodelled
plant, used to back-calculate the open-loop plant from closed-loop sweeps.

Absolute stability and describing-function methods are well established in aeroelasticity with
structural nonlinearity, and I am not claiming anything new against that literature. What this memo
adds is the quantified demonstration on one validated rig, plus the verification record in Section 10.

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

The three gap sizes collapse onto a single curve. Below the jump at 0.50 $U_f$ the largest difference
between the 1.15° and 2.12° gaps across all normalized amplitudes is under 10⁻⁹. That collapse is
Conner's analytical scaling claim, and I did not impose it. At the jump itself the curves differ by up
to 6% in flap amplitude, and by 0.01 to 0.4% above it: path-following carries the state from one speed
to the next, and at a branch switch rounding differences between the two runs grow into a different
transient. The August run, at a different time step, happened to agree to 2.3 × 10⁻¹² everywhere.

![**Figure 1.** Freeplay path-following at three gap sizes, with amplitudes normalized by gap size.
The three curves collapse onto one. The dotted line marks the jump to the high-frequency limit cycle
near 0.50 $U_f$.](fig3_freeplay_sweep.pdf){width=95%}

## 4. The identified margin depends on the test amplitude

### 4.1 Across records

Five subcritical speeds from 0.10 to 0.26 $U_f$, shaker raps of varying size, PO-MOESP
identification, then a linear extrapolation of damping to zero.

| max \|β\|/gap | extrapolated $U_F$ / true $U_F$ |
|---|---|
| 0.10 | 0.383 |
| 0.96 | 0.382 |
| 1.51 | 0.258 |
| 2.78 | 0.454 |
| 6.00 | ∞, damping rises with airspeed |
| 19.32 | ∞ |

The first two rows agree to 0.001. Below one gap width the flap never leaves the deadband, so those
two tests are the same experiment. Once the flap reaches the stops the identified damping rises by a
factor of six and the shape of the curve inverts.

Every finite prediction here is conservative, by 2.2× to 3.9×. But that is a property of this
configuration and I have not shown it holds in general. **A procedure whose answer moves by a factor
of 1.8 with a parameter nobody records is not usable for certification, whichever way it errs.**

The August version of this table read 0.381, 0.381, 0.259 and 0.544. The simulation was sampled at
833.33 Hz rather than the rig's 800 Hz (Section 10, item 15), which moved every identified frequency
by 4%; the 2.78 row moved most.

### 4.2 Within a single record

A decaying record on a freeplay structure is not a stationary system. Early in the decay the flap
hits the stops. Late in the decay it rattles around inside the gap. I slid a 1.5 s window in 0.25 s
steps down one ring-down at 0.26 $U_f$, gap 2.12°, 40 dB SNR. That gives 50 estimates from one
record.

| quantity | range across windows |
|---|---|
| flap amplitude | 2.154 down to 1.060 gap widths |
| identified ζ | −0.0075 to +0.0824 |
| identified frequency | 4.52 to 26.09 Hz |

![**Figure 2.** One ring-down at $U/U_f$ = 0.26, gap 2.12°, 40 dB SNR. (a) is the flap response with
the gap edges marked. (b) and (c) are the identified damping and frequency from 50 overlapping
windows of that same record, coloured by window index. The record does not decay to zero. It settles
onto a limit cycle.](fig1_within_record.pdf){width=95%}

One structure, one airspeed, one record. Two engineers looking at different windows of the same data
get modal frequencies that differ by a factor of 5.8, and damping ratios that disagree on sign.

The record also does not decay to zero. It settles at |β|/gap = 1.074 with ζ near zero. That is a
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

### 6.1 Exact sector decomposition about an operating point

The only thing that differs between the stiff and free regions is the flap restoring moment. So the
plant splits exactly into a linear system with a memoryless nonlinearity in feedback:

> ż = $A_{\mathrm{free}}$·z + $B_\varphi$·(−φ(β)) + $f_0$,  β = $C_\beta$·z
>
> $A_{\mathrm{stiff}}$ = $A_{\mathrm{free}}$ − $B_\varphi$·$K_\beta$·$C_\beta$  (verified to 3.6 × 10⁻¹²)

where $f_0$ is any constant load. The August version set $f_0$ = 0 and certified the origin. A
control surface is under load in practice, so it is pushed against one side of the gap and sees the
freeplay only when it swings back through it. Then the origin is not an equilibrium, and a
certificate about the origin says nothing.

The equilibrium $z^*$ solves 0 = $A_r z + a_r + f_0$ in a region r that contains β\*. I modelled two
loads: a bench hinge moment $M_0$, written as $p = M_0/(K_\beta\delta)$, the number of gap widths it
would push the flap past the stop in still air; and an aerodynamic load from a mounting angle
$\alpha_0$, which enters every aerodynamic term that pitch does but not the pitch spring, and scales
with $U^2$. Exactly one region-consistent equilibrium exists at every speed and load tested, with a
residual of 4.5 × 10⁻¹³.

Shifting to $z = z^* + x$ gives

> ẋ = $A_{\mathrm{free}}$·x − $B_\varphi$·φ_s($C_\beta$ x),  φ_s(s) = φ(β\* + s) − φ(β\*).

φ is monotone with slope in [0, $K_\beta$], so φ_s is too, and φ_s(0) = 0. The linear part is
unchanged. **So the sector argument survives, but what it certifies is the operating point $z^*$,
not the origin.**

Two things change. First, φ_s is not odd unless β\* = 0, so Zames–Falb multipliers must keep
h(t) ≥ 0 (Turner 2021, Theorems 4 and 5). Second, and more useful, if the flap rests on a stop
(|β\*| > δ) the sector narrows. The smallest ratio φ_s(s)/s occurs when the swing just reaches the
far edge of the gap:

> $k_1 = K_\beta\,(|\beta^*| - \delta)/(|\beta^*| + \delta)$.

A sector [$k_1$, $K_\beta$] with $k_1$ > 0 lets a loop transformation move $k_1$ of the flap stiffness
into the linear part, and the Popov test then runs on a plant that is much further from instability.

### 6.2 Bound hierarchy at zero preload

| criterion | certified speed | status |
|---|---|---|
| circle | nothing at any speed | rigorous, but allows time-varying φ |
| **Popov** | **0.1612 $U_f$** | **rigorous** |
| describing-function minimum | 0.1612 $U_f$ | heuristic, minimum at γ = 0.050 |
| Zames–Falb alone | nothing | see below |
| Popov + Zames–Falb | 0.1612 $U_f$ | rigorous; no gain over Popov |
| free-region limit | 0.2858 $U_f$ | hard ceiling, $A_{\mathrm{free}}$ loses stability |
| observed LCO, 2.12° / 1.15° | 0.18 / 0.20 $U_f$ | simulation |
| Conner, numerical / experimental | 0.18 / 0.23 $U_f$ | published |

![**Figure 3.** (a) Instability onset over the effective-stiffness family
$K_{\mathrm{eff}}$ = γ$K_\beta$, showing the interior minimum at γ = 0.050. (b) The bound hierarchy on
a single axis. The shaded region is certified by the Popov multiplier.](fig2_bounds.pdf){width=88%}

**Below 0.1612 $U_f$, no disturbance of any size can produce sustained motion, for any nonlinearity
in the sector.** The bound needs no knowledge of amplitude, waveform or gap size. It sits below both
limit-cycle onsets I observed and below both onsets Conner et al. report. The 0.18 $U_f$ onset was
re-checked this revision with 48 random starts carrying both displacement and velocity: 34 end on a
cycle of 1.22 gap widths at 0.18 $U_f$, none at 0.16 $U_f$.

**Why circle fails and Popov works.** The circle criterion allows the nonlinearity to be
time-varying. Freeplay is memoryless and time-invariant. Using that extra piece of information is the
entire difference between certifying nothing and certifying to within 10% of the observed onset.

**Why Zames–Falb adds nothing.** Zames–Falb multipliers use monotonicity, which freeplay has. But
with ‖h‖₁ ≤ 1 the multiplier satisfies |M(jω)| ≤ 2, so it can only supply bounded phase advance. The
Popov multiplier 1 + jωη is unbounded, and on this loop the binding requirement turns out to be
phase. Adding a Zames–Falb term to the Popov multiplier, the combination used in the literature
(Carrasco et al. 2014), certifies 0.1612 $U_f$ in both the odd and the h ≥ 0 classes.

**The August combined bound was invalid.** The August version multiplied the two multipliers,
(1 + ηs)(1 − Ĥ(s)), and got 0.2371 $U_f$, above the 0.18 $U_f$ limit cycle. I found no theorem that
covers that product, and the plant falsifies it: a valid global certificate cannot be beaten by a
sustained cycle. Restricting to h ≥ 0 does not rescue it (0.2356 $U_f$), so the product, not the odd
class, was the error. This closes the first open item of the August memo.

**The agreement with the describing function is not a theorem.** Popov and the describing-function
minimum agree to four digits on this configuration. Take the structural damping out and they diverge
completely. The describing-function minimum drops to 0.00209 $U_f$ while Popov certifies nothing
above about 0.005 $U_f$. The agreement is specific to this configuration.

### 6.3 Under preload

| preload | β\*/δ at 0.3 $U_f$ | $k_1/K_\beta$ at the certified speed | Popov on [$k_1$, $K_\beta$] | describing-function onset | simulated onset |
|---|---|---|---|---|---|
| none | 0 | 0 | 0.1612 | 0.18 | 0.18 |
| hinge, p = 0.25 | +1.21 | 0.10 | 0.2079 | 0.50 | 0.68 |
| hinge, p = 0.5 | +1.46 | 0.18 | 0.3139 | none to 0.70 | none to 0.70 |
| hinge, p = 1 | +1.94 | 0.29 | 0.5590 | none to 0.70 | none to 0.70 |
| hinge, p = 2 | +2.91 | 0.41 | 0.9681 | none to 0.70 | none to 0.70 |
| aero, α₀ = 2° | −0.58 | 0 | 0.1612 | 0.18 | 0.18 |
| aero, α₀ = 4° | −1.00 | 0.0006 | 0.1612 | 0.18 | 0.18 |
| aero, α₀ = 6° | −1.02 | 0.003 | 0.1612 | 0.18 | 0.18 |

Onsets are $U/U_f$ on a 0.02 grid, simulated from 16 random starts per speed with displacements up to
12 gap widths and velocities. Every one of the 1,392 trajectories below a certified speed returned to
$z^*$ to within 2 × 10⁻⁶ δ after 40 s. Figure 4 shows the same data.

![**Figure 4.** (a) Onset of sustained motion against hinge preload. The line is the Popov certificate
on the preload-narrowed sector; the shaded region below it is certified. Hollow triangles mark preloads
where no cycle was found up to the 0.70 $U_f$ search ceiling. (b) Limit-cycle amplitude against speed:
filled markers are simulation, hollow markers are every dual-input describing-function solution,
stable or not.](fig4_preload.pdf){width=100%}

Four things follow.

**A bench preload is stabilising, and the certificate shows it.** Half a gap of preload removes the
limit cycle up to at least 0.70 $U_f$ in both simulation and describing function, and two gaps
certify the operating point almost to the linear flutter speed. The certified speed rises with
preload because a flap held on a stop sees a stiffness of at least $k_1$, not zero.

**An aerodynamic preload is not.** On this rig the aerodynamic hinge moment and the aerodynamic hinge
stiffness both scale with $U^2$, so a mounting angle floats the flap at a fixed fraction of the gap,
0.28 gap widths per degree, at every speed. Below α₀ ≈ 3.5° it never reaches a stop. Above that it
only touches one, and $k_1$ stays near zero. Onset stays at 0.18 $U_f$ and cycle amplitudes grow by
up to 18%. Whether preload helps depends on whether the load source has its own stiffness.

**The describing function gets the onset and misses the amplitude.** With zero preload its stable
branch is 4–6% below the simulated amplitude. Under aerodynamic preload the gap grows to about 20%
(1.21 against 1.50 gap widths at 0.30 $U_f$, α₀ = 4°). It also predicts onset at 0.50 $U_f$ for
p = 0.25, where simulation finds none below 0.68.

**The describing function needs the preload, but its error is not the preload.** Scaling the
aerodynamic preload it is given by ±25% moves its amplitude by under 1%, against the 19% it already
misses by. Bench preload is the opposite case: whether the flap is on the stop at all decides whether
there is a cycle, so it has to be measured, which is the bench-versus-sweep question I put to Ting.

## 7. Measuring model error near onset

The memo's certificate needs a sector, which the plant supplies exactly. A robust controller needs
something the plant does not supply: a bound on how wrong the linear model is. This section measures
that bound the way a tunnel campaign would, from many forced runs on a plant that is deliberately not
the model, and asks where it can be used. The runs perturb stiffnesses, inertias and damping ratios by
4% (flutter speed +7.75%), excite with a periodic multisine through the shaker, and estimate the
frequency response at the excited lines from integer periods.

**The measurement resolves model error, with a known inflation.** With the flap inside the gap the
true plant is exactly the free-region linear system, so the measured mismatch can be checked against
the true one. At 40 dB it over-reads by 2.0 to 3.2×, the safe direction; in the September runs the
excess fell toward zero as the noise was reduced, which identifies it as the noise floor. Twelve
replicates per condition put the scatter at 3 to 5% of the measured value, so signal-to-floor is 6.6
to 16.

**What the August certificate can absorb: about 2%.** Robustifying Popov with a flat multiplicative
weight W leaves 0.1565 $U_f$ at W = 0.02 and certifies nothing at 0.03, because at the binding
frequency |G| = 32.8 against 1/k = 0.257. The true model error in this loop is 3 to 4% (median over
the band), and the measured one 7 to 14%. **Any honestly measured W destroys the certificate,
while the unrobustified certificate survives 32% parameter error**, because Popov at 0.161 $U_f$
already sits far below the perturbed plant's onset. The measured W and this certificate are not in
series. W belongs to a robust controller design.

**In a robust controller design, the unstructured weight cannot be met.** I measured W with an
antiresonance floor (relative error is meaningless where |G| → 0) and tested robust stability,
‖W T_i‖∞ < 1, on H∞ designs over a grid of actuator authority and performance weight:

| $U/U_f$ | default design | best design found | stabilises the true plant |
|---|---|---|---|
| 0.30 | 2.32 | 0.46 | yes |
| 0.60 | 2.49 | 1.32 | yes |
| 1.00 | 3.28 | 2.43 | yes |
| 1.20 | 4.66 | 2.25 | yes |

The violation sits at 10.4 Hz, the free-region mode, which the 4% perturbation shifts. Near a
lightly damped mode a frequency error of a few percent is a relative error above 1 (here 2.2). That is
real model error, not an artifact, and a low-order weight cannot follow it. Every design stabilises
the true plant, so the test is conservative rather than violated. The error the data show is the
frequency shift of one mode; it needs a structured, real-parametric description, not a larger
unstructured weight.

**Two further results limit what a measured weight can do.** W measured through the shaker covers
the true error in the hinge-moment channel, the channel the test is about, at only 28% of frequencies;
it has to be measured through the control surface. And at 0.30 $U_f$ the design that meets the
measured W sustains a 1.4-gap limit cycle on the true freeplay plant from two of four starts, while
the default design, which fails the test, settles from all four. A linear robustness test says nothing
about the freeplay cycle. Only the absolute-stability side of this memo does.

**Measurement quality approaching onset depends on dwell time.** Holding the flap inside the gap and
approaching the free-region instability, I measured signal-to-floor against the damping ratio of the
plant being tested. With a fixed protocol that discards one 16 s period, it falls from 1.6–5.5 to
0.3–1.1 between ζ = 0.004 and ζ = 0.0024, which is where the flutter mode's decay time overtakes the
discarded period. If instead each point waits five decay times before averaging, it declines
gradually: 3 or more in every channel down to ζ ≈ 0.0016, 2.6 in the flap channel at ζ = 0.001 and
1.7 at ζ = 0.00045, with pitch above 5 throughout. The cost is the dwell, which grows as 1/(ζω) to
three minutes per point at ζ = 0.001 and seven at 0.00045, and a bias that grows from about 3× to
5–11×. An open-loop test near onset is not blind, but it is slow and increasingly biased, which is
the practical form of the instability input–output problem. (Ting was told in September that the
quality "falls off a cliff" between ζ ≈ 0.003 and 0.001. Both the cliff and the ζ values were wrong:
the cliff was the fixed protocol, and the ζ values were those of the nominal model rather than the
plant being measured, which is 1.4 to 2 times better damped at the same speeds.)

## 8. Identifying the plant through flutter, in closed loop

Past the flutter speed an open-loop test is impossible, and near it the previous section shows it is
slow. The field's answer is to run the flutter suppression system, sweep, identify the closed loop,
and back-calculate the open-loop plant. I did that in simulation, where the truth is known:

- **Truth.** The 4%-perturbed plant (flutter at 1.0775 $U_f$ of the nominal model), a flap hinge-moment
  actuator with a 30 Hz lag, unmeasured turbulence, RVDT noise and a 16-bit ADC.
- **Controller.** A discrete LQG at 800 Hz designed at each speed on the nominal model, with the
  designer's actuator estimate of 40 Hz. The as-built controller adds a 5% gain error and a 100 Hz
  output filter that the design document does not show.
- **Test.** A periodic multisine added to the controller output, u = K y + r; u and y logged.
- **Estimators.** The direct method from u and y alone; the indirect method, T_yr/(1 + K T_yr), with
  the as-built controller and with the documented one; and the joint input–output method,
  T_yr/T_ur, which uses r as the instrument and needs no controller model.

The as-built loop is stable on the true plant from 0.50 to 1.20 $U_f$ in 0.01 steps. One result
appeared on the way. Near 1.05 $U_f$ the nominal model is past its flutter speed and the true plant is
not, and the model-based controller fights an instability that is not there: it lowers the true
flutter-mode damping from 0.016 to 0.008. A measured plant is what would let the designer see that.

| | below flutter | above flutter |
|---|---|---|
| joint input–output, pitch / flap / plunge | 2% / 0.4% / 5–10% | 2% / 0.3% / 3–4% |
| indirect, as-built controller | 2–3% / 0.6–1.6% / 5–10% | 2% / 0.5–0.6% / 3–4% |
| indirect, documented controller | 40–43% / 34–42% / 39–43% | 42–49% / 42–48% / 44–52% |
| open loop, settled | 0.5% / 0.1% / 1–3% | not possible |

(Median error over the excited lines, 0.5–30 Hz, six replicates.)

**On Ting's first question.** If the controller is known exactly, the indirect method works as well
as the joint input–output method. If the analyst uses the documented controller, a 5% gain error
and an unmodelled filter become a 34–52% plant error, because the controller error is multiplied by
the loop gain, |KG| ≈ 15 at the flutter frequency. The joint input–output method gives the same answer
as an exactly known controller without needing one. In simulation, then, dividing out a measured
controller response is not optional unless the controller is characterised to a precision the loop
gain will tolerate.

**The flutter mode is the hard part.** At the flutter frequency the loop drives the plant input down
to |T_ur| ≈ 0.07 of the excitation, so the plant is least observed exactly where the information is.
An unweighted rational fit chases that noise and returns spurious, nearly undamped poles. Weighting
each line by its period-to-period scatter and fitting the replicate average fixes the pole locations
but not all of the bias:

| $U/U_f$ | 0.80 | 1.00 | 1.05 | 1.10 | 1.20 |
|---|---|---|---|---|---|
| true flutter-mode ζ | +0.064 | +0.039 | +0.016 | −0.015 | −0.080 |
| joint input–output | +0.072 | +0.056 | −0.023 | −0.027 | −0.085 |
| open loop, settled | +0.067 | +0.044 | +0.021 | — | — |

The back-calculated damping crosses zero at 1.036 $U_f$ against a true 1.0775, 3.9% early. The
conventional answer, extrapolating the last two settled open-loop points, gives 1.096: 1.7% late,
the unsafe direction, and it cannot check itself past the boundary. The worst closed-loop point is
1.05 $U_f$, where the controller itself is barely damped (0.008).

**With freeplay, the plant you identify is set by the amplitude.** Preloaded one gap width onto its
stop, the flap stays there and the estimate matches the stiff plant to 0.4–2.3%. With a symmetric gap
the flap rattles across it at about 1.3 gap widths and the estimate matches the free-region plant to
3–4%, and the stiff plant not at all, because at that amplitude the effective stiffness is close to
zero. A closed-loop test on a freeplay structure returns whichever linear plant its amplitude selects.

## 9. Which identification method

All comparisons use the rig's own instrumentation: three RVDTs measuring pitch, flap and plunge
displacement, 800 Hz sampling, 4000-point records, and a shaker at midspan quarter chord with
measured force. Ground truth is eig A(U) of the plant that generated the data. The metric is mean
relative error in modal damping over the three modes. Modes are matched to truth by optimal
assignment on relative frequency, after discarding candidates with ζ > 0.5, which the true modes of
this plant never approach. Every method is scored the same way.

### 9.1 What model order is even admissible

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

### 9.2 Results

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

## 10. Verification record

I found and fixed twenty-seven errors while doing this work, fourteen before August and thirteen in
the October revision. They are listed because the headline result is
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
| 14 | Section 9 originally compared methods at order 2 | Order 2 admits no exact model of this plant, so the comparison was run in a model class that cannot represent it |

| 15 | Simulation sampled at 833.33 Hz, not 800: dt = 2 × 10⁻⁴ does not divide 1/800 | Measured model error was inflated at every lightly damped resonance; the identified frequencies were all 4% high. Moved the Section 4.1 margin at 2.78 gap widths from 0.544 to 0.454 $U_f$ |
| 16 | Welch averaging applied to a periodic multisine | Smeared the excitation lines and inflated measured mismatch 3 to 4×; replaced by integer periods and DFT at the excited lines |
| 17 | Measured W compared against the stiff nominal | Counted the freeplay twice, once in the sector and once in W: 1.57 against a true 0.04 in the flap channel |
| 18 | Popov × Zames–Falb product used as a multiplier | Certified 0.2371 $U_f$ against a sustained cycle at 0.18 $U_f$; no theorem covers the product. The sum certifies 0.1612 |
| 19 | Near-onset damping quoted from the nominal model, not the plant being measured | The two differ by 1.4 to 2× near onset; every ζ in the September cliff was mislabelled |
| 20 | The near-onset "cliff" was the test protocol | Waiting five decay times before averaging restores resolution down to ζ ≈ 0.001 |
| 21 | A near-onset "collapse" row was freeplay | The flap had left the gap: 950 switching events in a run labelled linear |
| 22 | A closed loop called bistable | The H∞ design had a tenth of the actuator authority it needed; with enough, every start converged |
| 23 | Failure of the robust-stability test blamed on antiresonance bins | With a floor on those bins it still fails; the cause is a 10.4 Hz mode shifted by the perturbation |
| 24 | Shaker-channel W said to cover the hinge-channel error at 98–100% of frequencies | Compared correctly, 28% |
| 25 | Limit-cycle search used flap-position kicks only | Missed the subcritical cycle, which needs initial velocity, and put the zero-preload onset at 0.30 $U_f$ instead of 0.18 |
| 26 | Closed-loop modal fit unweighted | Returned spurious, nearly undamped poles where the loop starves the plant input |
| 27 | Flutter-suppression loop checked on a coarse speed grid | The as-built controller was unstable on the true plant at 1.05 $U_f$, between grid points |

Items 4, 7, 10, 11, 13, 14, 18, 19, 20 and 24 would each have produced a published claim that is false
rather than just imprecise. Items 19 and 20 were sent to a reader before they were caught.

## 11. Open items

**The certificate is gridded, not analytic.** The multiplier searches enforce their constraints on a
frequency grid, so they can miss violations between grid points. Real rigour would need a KYP/LMI
formulation.

**The bound is valid but not sharp.** At zero preload Popov certifies 0.1612 $U_f$ against a
simulated onset of 0.18. The August memo read the onsets of 0.18 $U_f$ (2.12° gap) and 0.20 $U_f$
(1.15°) as gap dependence the sector argument throws away. That was wrong: the system is exactly
homogeneous in the gap (Section 3.4), so with disturbances scaled to the gap the onset cannot depend
on it. The 0.20 came from disturbances that were not scaled.

**The describing function misses amplitude under aerodynamic preload by about 20%,** and predicts a
cycle at p = 0.25 that simulation does not find below 0.68 $U_f$. The preload search stops at
0.70 $U_f$, so the p = 0.5 and p = 1 cases are certified to 0.31 and 0.56 $U_f$ and cycle-free to
0.70, but not examined between 0.70 and the linear flutter speed.

**Closed-loop damping is biased near the flutter speed.** The plant frequency response comes back to
a few percent, but the extracted flutter-mode damping is off by up to 0.039 at 1.05 $U_f$, where the
flutter suppression loop is itself barely damped. A higher-authority excitation, or a fit that uses
the closed-loop structure instead of the back-calculated plant, are the obvious next steps.

**The uncertainty a robust controller needs is structured.** The measured error is the frequency
shift of one mode, which no unstructured weight can bound without forbidding control near that mode.
A real-parametric (μ) description, measured through the control surface rather than the shaker, is
untested.

**Order and bandwidth selection are not costed.** I picked the identification model order knowing the
true system dimension. On real data you know neither.

**Model fidelity.** Theodorsen incompressible two-dimensional aerodynamics with a two-lag Jones
approximation, which contributes about 0.75% error on flutter speed. Conner et al. attribute the 15%
gap between their numerical and experimental flutter speeds to aerodynamics outside the Theodorsen
model and to three-dimensional tunnel effects, with a measured 12.3% drift in pitch torsional
stiffness over the test campaign as part of the cause. This work reproduces their numerical model,
not their experiment. The "true" plant of Sections 7 and 8 is the same model with 4% parameter
scatter, which is a stand-in for real model error, not a measurement of it.

**Not attempted.** The nonperiodic regime between 0.35 and 0.44 $U_f$, which needs a Poincaré map
rather than RMS amplitude. Any experimental data of my own.

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
| `zames_falb.py` | Zames–Falb and Popov + Zames–Falb search; odd and h ≥ 0 classes |
| `preload.py` | Operating point under preload, narrowed-sector Popov, dual-input describing function, random-start simulation (Section 6) |
| `stage2a_harness.py` | Perturbed nominal model and periodic-multisine test harness (Section 7) |
| `stage2_uncertainty.py` | Measured W, robust Popov, antiresonance floor, fitted weight (Section 7) |
| `hinf.py`, `stage2_hinf.py` | Two-Riccati H∞ synthesis; robust-stability test against the measured W (Section 7) |
| `stage3_scatter.py` | Replicate scatter and measurement quality approaching onset (Section 7) |
| `item1_residual_sparsity.py` | The sparse-recovery negative result referred to in correspondence |
| `closed_loop_id.py` | Flutter suppression loop and closed-loop back-calculation (Section 8) |
| `make_figures.py` | Generates all four figures; every plotted value recomputed at run time, except the simulated preload amplitudes, read from `preload.py`'s output |

Every validation number in Section 3 can be reproduced by running these files. `make_figures.py`
recomputes each plotted quantity from the model instead of reading stored results, and its output is
reproducible across runs.

## References

Carrasco, J., Maya-Gonzalez, M., Lanzon, A. & Heath, W. P. (2014). LMI searches for anticausal and
noncausal rational Zames–Falb multipliers. *Systems & Control Letters* **70**, 17–22.

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

Turner, M. C. (2021). Zames–Falb multipliers: don't panic. arXiv:2106.15913.
