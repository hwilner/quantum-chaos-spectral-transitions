# Current results and discussion

## Synthetic validation (complete)

All 16 unit tests pass. Headline anchors:

| System | ⟨r⟩ measured | ⟨r⟩ theory |
|---|---|---|
| Poisson (independent levels) | 0.386 ± 0.02 | 0.38629 (exact: 2 ln 2 − 1) |
| GOE (pooled over realisations) | 0.534 ± 0.02 | ≈ 0.5359 |
| Spectral form factor, Poisson | flat ≈ 1 | flat |
| Spectral form factor, GOE | ramp + plateau | ramp + plateau |

## The pitfall that shaped the design

Early testing revealed that **white-noise (Wishart) covariance matrices are
already GOE-like in the bulk** (⟨r⟩ ≈ 0.50–0.53 for p = 60, T = 2 000): level
repulsion is a property of *any* empirical correlation matrix, not proof of
interesting brain coupling. Two consequences:

1. The science question must be a **contrast** (task vs rest, patient vs
   control), never an absolute claim of "the brain is GOE".
2. Null #3 (degrees-of-freedom-matched Wishart) is not optional — it is the
   baseline. Additionally, ratio statistics must be **pooled across
   realisations** rather than averaged per-realisation (a subtle point our
   test suite now guards: concatenating independently sorted spectra mimics
   Poisson statistics and must never be done).

## Interpretation sketch for H1

If the brain's task state strengthens structured coupling, random-matrix
theory predicts motion *toward* the GOE anchor. The size of the shift we can
detect: with 360 levels, bootstrap SE(⟨r⟩) ≈ 0.01 per scan, so paired
task–rest designs detect shifts of ≈0.005–0.01 — well below the Poisson–GOE
distance of 0.15, giving ample dynamic range.

## Risks

- **Haemodynamic smoothing** couples neighbouring parcels trivially → inflated
  repulsion. Mitigation: motion/confound regression, Schaefer-200 replication,
  and checking that the shift survives within-subject normalisation.
- **Global signal regression** (controversial) can introduce negative
  correlations; we run the pipeline with and without it and report both.
- **Multiple comparisons** across 7 tasks × several statistics: Bonferroni
  within each hypothesis, permutation tests for H3.

## Open questions

- Does ⟨r⟩ vary more across *states* within a subject or across *subjects*
  within a state? (This variance ratio decides whether the fingerprint or the
  state signal dominates.)
- Is the transition coordinate related to the criticality/synergy measures of
  our sister repository `criticality-synergy-brainwide`? A joint analysis is
  queued in the rotation.
