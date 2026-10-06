# Status and plan

## Done

- [x] Repository scaffold, license, requirements.
- [x] `spectral_chaos` package: correlation matrices, eigen-spectra,
      unfolding, spacing ratios, mean ratio, spectral form factor, Wigner
      surmise, GOE/Poisson synthetic anchors.
- [x] 16-test suite, all anchors green (Poisson ⟨r⟩ ≈ 0.38629; GOE
      ⟨r⟩ ≈ 0.536; SFF flat vs ramp; Wishart-repulsion pitfall documented).
- [x] Publication-quality INTRODUCTION and high-school EXTENDED_INTRODUCTION.

## In progress

- [ ] HCP data access and downloader (Issue #1).
- [ ] Dynamic/task FC pipeline (Issue #2).
- [ ] Spectral statistic computation over the cohort (Issue #3).

## Planned (see Issues #1–#7)

| # | Milestone | Depends on |
|---|---|---|
| 1 | HCP downloader + parcellation | — |
| 2 | Rest/task FC pipeline | 1 |
| 3 | Spectral statistics over cohort | 2 |
| 4 | Task-vs-rest transition test with nulls | 3 |
| 5 | Fingerprint reliability (ICC + ID rate) | 3 |
| 6 | Behaviour prediction (ΔR²) | 3 |
| 7 | Null models consolidation + figures | 4–6 |

## Next quarter target

H1 decided (task-minus-rest ⟨r⟩ shift, null-controlled) on the full S1200
sample, with the fingerprint analysis running in parallel.
