# Quantum-Chaos Spectral Transitions in Brain Networks

**Do the eigenvalue statistics of functional brain networks shift between
"regular" (Poisson-like) and "chaotic" (level-repelling, random-matrix)
regimes across rest and task — and does that position predict cognition?**

A deep result of quantum chaos — the Bohigas–Giannoni–Schmit conjecture —
says the *fine structure* of a spectrum betrays the underlying dynamics:
classically chaotic systems show random-matrix level repulsion (levels
"avoid" each other), while integrable systems show uncorrelated Poisson
levels. The same diagnostics apply to *any* symmetric matrix — including
the functional-connectivity (FC) matrices of the human brain. This
project asks whether cognitively engaged brain states sit closer to the
level-repelling end of the spectrum than rest, whether individuals occupy
reproducible positions on this regular↔chaotic axis, and whether the
position predicts cognitive ability out of sample. No physical quantum
claim is made: random-matrix statistics are used as a principled,
parameter-light fingerprint of dynamical mixing.

- **Core hypothesis.** Task states show stronger level repulsion than
  rest (H1); individuals have reliable spectral fingerprints (H2);
  intermediate spectral chaos — not the maximum — predicts cognitive
  flexibility (H3).
- **Falsifier.** After matching null models (preserving autocorrelation,
  variance profile, matrix size and sample length), task and rest
  spectra are statistically indistinguishable, and spectral statistics
  add no out-of-sample predictive value for cognition.

## Repository layout

```
src/spectral_chaos/
    spectral.py  Correlation matrices, eigenvalue unfolding, adjacent-gap
                 ratios (Poisson 0.386 vs GOE 0.536), spectral form
                 factor, Wigner surmise, Poisson/GOE reference ensembles.
tests/
    test_spectral.py  15 unit tests against random-matrix ground truth.
docs/
    INTRODUCTION.md              Publication-quality motivation and theory.
    EXTENDED_INTRODUCTION.md     High-school-level walkthrough with links.
    METHODS.md                   Materials and methods (data, pipeline).
    STATUS_AND_PLAN.md           Task board snapshot.
    CURRENT_RESULTS_AND_DISCUSSION.md
```

## Quickstart

```bash
pip install -r requirements.txt
export PYTHONPATH=src        # or: pip install -e .
python -m unittest discover -s tests -v
```

```python
from spectral_chaos import (goe_matrix, poisson_spectrum,
                            sorted_eigenvalues, mean_spacing_ratio)

print(mean_spacing_ratio(poisson_spectrum(2000, seed=0)))        # ~0.386
print(mean_spacing_ratio(sorted_eigenvalues(goe_matrix(300))))   # ~0.536
```

## Empirical target

[HCP Young Adult](https://www.humanconnectome.org/study/hcp-young-adult):
1,200 participants with resting-state and 7 task fMRI sessions plus
extensive cognitive batteries. Per subject × condition we build dynamic
FC matrices, unfold their spectra, and compute spacing ratios, rigidity
and form factors — then test task-vs-rest shifts, fingerprint
reliability across scan days, and incremental prediction of cognition
over conventional FC features (see `docs/METHODS.md`).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). All code uses Google-style
docstrings and ships with unit tests.

## License

MIT (see LICENSE).
