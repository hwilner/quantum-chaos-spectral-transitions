"""Random-matrix / quantum-chaos spectral statistics.

The Bohigas-Giannoni-Schmit conjecture links the eigenvalue fluctuations
of classically chaotic quantum systems to random-matrix ensembles (GOE
for time-reversal-symmetric systems), while integrable systems show
Poisson statistics.  This module computes the standard diagnostics for
any real symmetric matrix (functional-connectivity, covariance,
Laplacian):

* nearest-neighbour spacing ratios  r_i = min(s_i, s_{i+1}) / max(...),
* the mean ratio (Poisson 0.3863 vs GOE 0.5359, Oganesyan-Huse 2007),
* spectral unfolding via the smooth part of the staircase,
* the spectral form factor  K(tau),
* reference ensembles and the Wigner surmise.
"""

from __future__ import annotations

import numpy as np

#: Mean adjacent-gap ratio of the Poisson (integrable) ensemble.
POISSON_R = 2.0 * np.log(2.0) - 1.0          # 0.38629...
#: Mean adjacent-gap ratio of the Gaussian Orthogonal Ensemble.
GOE_R = 0.53590


def correlation_matrix(timeseries: np.ndarray) -> np.ndarray:
    """Pearson correlation matrix of a multivariate time series.

    Args:
        timeseries: Array of shape ``(n_timepoints, n_channels)``.

    Returns:
        Symmetric correlation matrix of shape ``(n_channels, n_channels)``.
    """
    x = np.asarray(timeseries, dtype=float)
    if x.ndim != 2:
        raise ValueError("timeseries must be 2-D")
    return np.corrcoef(x, rowvar=False)


def sorted_eigenvalues(matrix: np.ndarray) -> np.ndarray:
    """Sorted eigenvalues of a real symmetric matrix.

    Args:
        matrix: Square symmetric matrix.

    Returns:
        Eigenvalues in ascending order.
    """
    m = np.asarray(matrix, dtype=float)
    if m.ndim != 2 or m.shape[0] != m.shape[1]:
        raise ValueError("matrix must be square")
    return np.linalg.eigvalsh(m)


def unfold_spectrum(eigs: np.ndarray, poly_degree: int = 3) -> np.ndarray:
    """Unfold a spectrum to unit mean level density.

    The staircase ``N(lambda)`` (number of eigenvalues below lambda) is
    fitted with a polynomial of degree ``poly_degree``; unfolded levels
    are ``e_i = N_fit(lambda_i)``, which have asymptotically constant
    unit density by construction.

    Args:
        eigs: Sorted eigenvalues (1-D).
        poly_degree: Degree of the smooth staircase fit.

    Returns:
        Unfolded eigenvalues with the same length.
    """
    eigs = np.sort(np.asarray(eigs, dtype=float))
    if eigs.size < poly_degree + 2:
        raise ValueError("too few eigenvalues to unfold")
    staircase = np.arange(1, eigs.size + 1)
    coeffs = np.polyfit(eigs, staircase, poly_degree)
    return np.polyval(coeffs, eigs)


def spacing_ratios(eigs: np.ndarray, unfolded: bool = False) -> np.ndarray:
    """Adjacent-gap ratios of a spectrum.

    ``r_i = min(s_i, s_{i+1}) / max(s_i, s_{i+1})`` with
    ``s_i = e_{i+1} - e_i``.  Ratios are unfolding-independent by
    construction and lie in ``[0, 1]``.

    Args:
        eigs: Sorted eigenvalues (1-D).
        unfolded: If True, ``eigs`` are already unfolded (no effect on
            the ratio; kept for API clarity).

    Returns:
        Array of ``len(eigs) - 2`` ratios.
    """
    eigs = np.sort(np.asarray(eigs, dtype=float))
    s = np.diff(eigs)
    s = s[s > 0]
    if s.size < 2:
        raise ValueError("need at least two positive spacings")
    lo = np.minimum(s[:-1], s[1:])
    hi = np.maximum(s[:-1], s[1:])
    return lo / hi


def mean_spacing_ratio(eigs: np.ndarray) -> float:
    """Mean adjacent-gap ratio (Poisson 0.386, GOE 0.536).

    Args:
        eigs: Sorted eigenvalues (1-D).

    Returns:
        Mean of :func:`spacing_ratios`.
    """
    return float(np.mean(spacing_ratios(eigs)))


def spectral_form_factor(eigs: np.ndarray, taus: np.ndarray) -> np.ndarray:
    """Spectral form factor ``K(tau)`` of an unfolded spectrum.

    ``K(tau) = (1/N) |sum_j exp(-2 pi i e_j tau)|^2``.  For a Poisson
    spectrum ``K(tau) -> 1`` (up to the tau = 0 peak); random-matrix
    ensembles show a dip-ramp-plateau.

    Args:
        eigs: Unfolded eigenvalues (1-D).
        taus: Non-negative times at which to evaluate K.

    Returns:
        ``K(tau)`` evaluated at each tau.
    """
    e = np.asarray(eigs, dtype=float)
    taus = np.asarray(taus, dtype=float)
    phases = np.exp(-2j * np.pi * np.outer(taus, e))
    return np.abs(phases.sum(axis=1)) ** 2 / e.size


def wigner_surmise_pdf(s: np.ndarray, beta: int = 1) -> np.ndarray:
    """Wigner-surmise spacing pdf for the GOE (beta=1) or GUE (beta=2).

    ``p(s) = a s^beta exp(-b s^2)`` with ``a, b`` chosen so the
    distribution has unit mean and unit normalisation.

    Args:
        s: Spacings (>= 0).
        beta: Dyson index (1 = GOE, 2 = GUE).

    Returns:
        Probability density at each ``s``.
    """
    from scipy.special import gamma as gamma_fn

    s = np.asarray(s, dtype=float)
    b = (gamma_fn((beta + 2) / 2) / gamma_fn((beta + 1) / 2)) ** 2
    a = 2 * b ** ((beta + 1) / 2) / gamma_fn((beta + 1) / 2)
    return a * s ** beta * np.exp(-b * s**2)


def goe_matrix(n: int, seed: int = 0) -> np.ndarray:
    """Sample an n x n Gaussian Orthogonal Ensemble matrix.

    Args:
        n: Matrix dimension.
        seed: Random seed.

    Returns:
        Real symmetric matrix ``(A + A^T) / sqrt(2)`` with standard
        normal ``A``.
    """
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n))
    return (a + a.T) / np.sqrt(2.0)


def poisson_spectrum(n: int, seed: int = 0) -> np.ndarray:
    """Sample n Poisson-distributed levels (i.i.d. exponential spacings).

    Args:
        n: Number of levels.
        seed: Random seed.

    Returns:
        Sorted levels with independent exponential spacings (unit mean).
    """
    rng = np.random.default_rng(seed)
    return np.cumsum(rng.exponential(1.0, size=n))
