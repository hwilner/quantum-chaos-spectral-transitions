# Current results and discussion

## Synthetic validation (complete)

All 43 unit tests pass. Headline anchors:

| System | ⟨r⟩ measured | ⟨r⟩ theory |
|---|---|---|
| Poisson (independent levels) | 0.386 ± 0.02 | 0.38629 (exact: 2 ln 2 − 1) |
| GOE (pooled over realisations) | 0.534 ± 0.02 | ≈ 0.5359 |
| Spectral form factor, Poisson | flat ≈ 1 | flat |
| Spectral form factor, GOE | ramp + plateau | ramp + plateau |

### Software verification record (2026-10-30)

What was executed, and what it covers. This is a record about the
*repository*, not about the brain.

| Item | Command | Outcome |
|---|---|---|
| Test suite, NumPy 2.4.6 | `python -m pytest -q` | 43 passed |
| Test suite, NumPy 1.26.4 | `python -m pytest -q` | 43 passed |
| Documented test count | `python tools/check_doc_test_count.py` | matches the suite (43) |
| CI, Python 3.11/3.12 × NumPy 1.x/2.x | GitHub Actions, workflow `Tests` | all jobs succeeded |
| Clean-checkout install | `pip install .` | imported `spectral_chaos` 0.1.0 |

The synthetic anchors of the table above were re-derived with the code as
committed and reproduce within the stated tolerances: Poisson ⟨r⟩ =
0.38778 (table: 0.386 ± 0.02), GOE ⟨r⟩ pooled over 40 realisations =
0.535899 (table: 0.534 ± 0.02), and the Poisson form factor has median
0.811816 over τ ∈ [0.2, 3.0] (table: flat ≈ 1). The form-factor row for
GOE is now anchored by a test that averages the statistic over 20
realisations: the GOE spectrum dips at small τ and then plateaus, whereas
the Poisson spectrum does not dip.

**What this record does not cover:** it contains no empirical brain data.
The HCP cohort has not been downloaded, no null model has been run, no
task-versus-rest contrast has been computed, and no quantity in the table
above was measured on a brain scan. The table remains a statement about
synthetic reference ensembles only.

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
