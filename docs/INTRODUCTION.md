# Quantum–chaos spectral transitions in functional brain networks

**Do task-engaged brain networks sit at a different point of the Poisson–GOE
spectral transition than rest — and is the position a fingerprint of the
individual?**

## Background

### Random matrix theory (RMT)

Random matrix theory is the statistical study of eigenvalues of matrices whose
entries are drawn at random (Mehta, 2004). Its founding observation — Wigner's,
in 1950s nuclear physics — is that the *fluctuation statistics* of complicated
interacting systems are universal: they depend only on the symmetries of the
system, not on its microscopic details. The most important ensemble for us is
the **GOE** (Gaussian Orthogonal Ensemble): real symmetric matrices with
independent Gaussian entries. The GOE is the canonical model of *quantum chaos*:
the eigenvalues of a classically chaotic quantum system (a "quantum billiard"
whose classical counterpart is chaotic) are conjectured — the
**Bohigas–Giannoni–Schmit (BGS) conjecture** (Bohigas, Giannoni & Schmit,
1984), supported by overwhelming numerical and experimental evidence — to have
the same statistics as GOE eigenvalues.

### Level repulsion and the spacing ratio

The single most famous GOE phenomenon is **level repulsion**: neighbouring
eigenvalues avoid each other. If s is the gap between two adjacent eigenvalues,
the probability of finding a tiny gap vanishes, p(s) → 0 as s → 0 — in stark
contrast to uncorrelated ("Poisson") levels, for which p(s) = exp(−s) is
*maximal* at s = 0. The integrable→chaotic transition is thus visible as a
continuous deformation of the spacing distribution from Poisson to Wigner-like
form.

Comparing spacing distributions across systems requires *unfolding* — rescaling
each spectrum so its mean level density is uniform — a step that is delicate
and error-prone. **Oganesyan and Huse (2007)** introduced a remarkably robust
alternative, the **ratio of adjacent spacings**: for each consecutive pair of
gaps (s_i, s_{i+1}),

```
r_i = min(s_i, s_{i+1}) / max(s_i, s_{i+1})   ∈ [0, 1].
```

Because the local density cancels in the ratio, ⟨r⟩ can be computed **without
any unfolding**, making it ideal for short, noisy, biological spectra. The
universal anchor values (Oganesyan & Huse, 2007; Atas et al., 2013) are

```
Poisson (uncorrelated levels):   ⟨r⟩ = 2 ln 2 − 1 ≈ 0.38629
GOE     (chaotic, repelling):    ⟨r⟩ ≈ 0.53590
```

with a continuous family in between for mixed/transition systems. For finite
samples we use the analytical Poisson value and the large-N GOE value, plus
bootstrap confidence intervals — we never treat the anchors as exact finite-N
targets.

### The Wigner surmise

Why do levels repel? Wigner's one-line derivation: in a GOE 2×2 matrix the
level-splitting depends on two independent Gaussian "coupling coordinates",
so the probability of a small splitting s is proportional to the available
phase space ∝ s — repulsion — times a Gaussian fall-off:

```
p(s) = (π/2) s exp(−π s² / 4)        (GOE, β = 1)
p(s) = exp(−s)                        (Poisson, no repulsion)
```

The general Wigner surmise p(s) = a s^β exp(−b s²) with repulsion exponent β
interpolates: β = 0 (Poisson), β = 1 (GOE), β = 2 (GUE, complex Hermitian),
β = 4 (GSE, symplectic). Fitting an effective β gives a second, complementary
transition coordinate. The **spectral form factor** K(τ) — the Fourier
transform of the two-point level correlation — is a third: K(τ) ≈ 1 at all τ
for Poisson, while GOE shows the celebrated "ramp" (linear growth) and
"plateau". Our library implements all three statistics.

### Why brain functional connectivity?

The functional-connectivity (FC) matrix of fMRI signals — correlations between
parcellated regional BOLD time series — is a real symmetric matrix, and its
eigenvalues are natural candidates for RMT analysis. Two established findings
motivate the transition hypothesis:

1. **FC eigenvalues already repel.** Wishart-type noise covariance matrices
   show GOE-like bulk statistics (a fact we verify in our test suite). The
   question is therefore quantitative: *how far* along the Poisson→GOE axis
   does the brain sit, and does the position *move* with cognitive state?
2. **Task modulates the correlation hierarchy.** During attention-demanding
   tasks, within-network correlations strengthen and global modularity shifts
   (Cole et al., 2014). Repulsion between FC eigenvalues is controlled by the
   strength and patterning of correlations — stronger, more structured coupling
   pushes the spectrum toward the GOE end, exactly as stronger hopping
   integrals drive the Poisson→GOE transition in condensed-matter systems.

## Core hypothesis

**H1 (state transition).** Task-engaged FC spectra sit closer to the GOE end
(higher ⟨r⟩, larger fitted β) than rest spectra, after controlling for temporal
autocorrelation with phase-randomised nulls.

**H2 (fingerprint).** The *position* ⟨r⟩_subject is a reliable individual
fingerprint across scan sessions (test–retest ICC), as other FC-derived
"connectome fingerprinting" measures are known to be (Finn et al., 2015).

**H3 (behaviour).** Within task, the transition coordinate predicts individual
differences in cognitive performance: subjects whose spectra repel more
strongly during demanding tasks perform better (or, in the alternative
"dynamical range" version, subjects whose task-minus-rest *shift* is largest).

**Falsifier.** If task and rest spectra are indistinguishable from their
respective matched nulls, or if the apparent shift is fully explained by
degrees-of-freedom loss (effective sample size) in task windows, the framework
is falsified at this spatial/temporal resolution.

## What is in this repository

`src/spectral_chaos/spectral.py` — the complete statistic toolkit with
Google-style docstrings:

| Function | What it computes |
|---|---|
| `correlation_matrix` | Fisher-z-ready FC matrix from regional time series |
| `sorted_eigenvalues` | eigenvalue spectrum of any symmetric matrix |
| `unfold_spectrum` | polynomial-staircase unfolding (for SFF and surmise fits) |
| `spacing_ratios` | Oganesyan–Huse ratios from a raw (not unfolded) spectrum |
| `mean_spacing_ratio` | ⟨r⟩ with a trimmed-mean option |
| `spectral_form_factor` | K(τ) from an unfolded spectrum |
| `wigner_surmise_pdf` | p(s; β) for β = 0, 1, 2, 4 |
| `goe_matrix`, `poisson_spectrum` | exact synthetic anchors |

`tests/` anchors every statistic against its theoretical value:
Poisson spectra give ⟨r⟩ ≈ 0.38629; GOE spectra (pooled over realisations, the
correct way to average ratios) give ⟨r⟩ ≈ 0.536; the SFF is flat for Poisson
and ramped for GOE; and — instructively — white-noise Wishart matrices are
*not* Poisson, which is why null models matter.

## Key terms

- **BOLD / fMRI**: blood-oxygen-level-dependent functional MRI; an indirect,
  slow (≈seconds) haemodynamic proxy for neural activity, measured
  simultaneously across the whole brain.
- **Functional connectivity (FC)**: statistical dependence between regional
  BOLD signals — operationally, the Pearson correlation matrix of parcellated
  time series. Used *when*: one wants a whole-brain snapshot of coupling;
  limitations: correlation ≠ causation, and FC mixes neural and
  vascular/motion structure.
- **Parcellation**: a partition of the cortex into regions ("nodes"). We use
  standard multi-modal parcellations (e.g. Glasser et al., 2016) with 100–400
  parcels — large enough for spectral statistics, small enough for stable
  correlations.
- **Eigenvalue / spectrum**: for a symmetric matrix M, the numbers λ with
  Mv = λv. The spectrum of FC summarises its correlation hierarchy: one big
  global mode, then network-scale modes, then a noisy bulk.
- **Level repulsion**: the tendency of neighbouring eigenvalues of interacting
  systems to avoid each other (p(s) → 0 at s → 0). *Why it matters*: it is the
  smoking gun of collective, strongly-coupled dynamics.
- **Poisson statistics**: the spacing statistics of *uncorrelated* levels
  (p(s) = e^{−s}); the universal signature of integrable/mixed/uncoupled
  systems.
- **GOE (Gaussian Orthogonal Ensemble)**: real symmetric random matrices with
  Gaussian entries; the universality class of time-reversal-symmetric quantum
  chaos. *When to use*: as the strong-coupling anchor for real symmetric FC.
- **BGS conjecture**: chaotic quantum systems have RMT level statistics. *Why
  here*: it is the license to read FC spectra with RMT eyes.
- **Spacing ratio ⟨r⟩**: unfolded-free measure of repulsion (0.386 Poisson,
  0.536 GOE). *Why*: robust for ~100–400 levels; no unfolding needed.
- **Unfolding**: rescaling a spectrum to unit mean spacing, needed for the
  SFF and spacing-distribution fits (not for ⟨r⟩).
- **Spectral form factor (SFF)**: Fourier transform of two-point level
  correlations; distinguishes Poisson (flat) from RMT (ramp + plateau) and
  probes longer-range spectral rigidity than ⟨r⟩.
- **Null model (phase randomisation)**: surrogate time series with identical
  power spectra but destroyed cross-dependencies. *Why*: to prove an effect is
  in the *coupling*, not in autocorrelation or finite-sample artefacts.
- **HCP (Human Connectome Project)**: the reference open fMRI dataset
  (rest + 7 tasks, ~1 200 young adults, densely sampled) — our test bed.
- **Test–retest ICC**: intra-class correlation between repeated measurements
  of the same subject; quantifies fingerprint reliability.

## Step-by-step: from raw spectra to the transition coordinate

1. Compute FC: C = corr(X), a p × p symmetric matrix.
2. Eigen-decompose: {λ_1 ≤ … ≤ λ_p}. Drop the trivial top global mode only if
   the analysis targets the *bulk*; always report both choices.
3. Form consecutive gaps s_i = λ_{i+1} − λ_i in the bulk (trim the edges where
   density varies fastest).
4. Form ratios r_i = min(s_i, s_{i+1})/max(s_i, s_{i+1}).
5. ⟨r⟩ = mean(r_i) with bootstrap CI over eigenvalues **and** over subjects.
6. Compare across states (rest vs task) against phase-randomised nulls.
7. If deeper confirmation is needed: unfold (polynomial staircase), fit the
   Wigner-surmise β, and compute the SFF ramp — three coordinates must agree.

## Related work and differentiation

RMT has been applied to brain data before — mostly as a *noise-cleaning* tool
for correlation matrices (Marchenko–Pastur edge filtering) or as a bulk
curiosity. This project instead treats the **spectral transition coordinate
itself as the observable**: a one-dimensional, theory-anchored, null-controlled
summary of how "collectively chaotic" the functional network is, tracked across
cognitive states and individuals. It is the brain analogue of the
integrable-to-chaos transition that defines much of modern quantum-chaos and
many-body-localisation physics (Oganesyan & Huse, 2007).

## References

- Mehta, M. L. (2004). *Random Matrices*. Elsevier.
- Bohigas, O., Giannoni, M.-J., & Schmit, C. (1984). Characterization of
  chaotic quantum spectra and universality of level fluctuation laws.
  *Physical Review Letters*, 52, 1.
- Oganesyan, V., & Huse, D. A. (2007). Localization of interacting fermions at
  high temperature. *Physical Review B*, 75, 155111.
- Atas, Y. Y., Bogomolny, E., Giraud, O., & Roux, G. (2013). Distribution of
  the ratio of consecutive level spacings in random matrix ensembles.
  *Physical Review Letters*, 110, 084101.
- Finn, E. S., et al. (2015). Functional connectome fingerprinting.
  *Nature Neuroscience*, 18, 1664.
- Cole, M. W., et al. (2014). Multi-task connectivity reveals flexible hubs
  for adaptive task control. *Nature Neuroscience*, 16, 1348.
- Glasser, M. F., et al. (2016). A multi-modal parcellation of human cerebral
  cortex. *Nature*, 536, 171.
- Van Essen, D. C., et al. (2013). The WU-Minn Human Connectome Project.
  *NeuroImage*, 80, 62.
