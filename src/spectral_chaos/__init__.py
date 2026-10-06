"""Quantum-chaos spectral statistics of functional brain networks.

Public API of :mod:`spectral_chaos`: correlation matrices, eigenvalue
unfolding, spacing ratios, spectral form factor, and reference ensembles
(Poisson / GOE) for null-model comparisons.
"""

from spectral_chaos.spectral import (
    GOE_R,
    POISSON_R,
    correlation_matrix,
    goe_matrix,
    mean_spacing_ratio,
    poisson_spectrum,
    sorted_eigenvalues,
    spacing_ratios,
    spectral_form_factor,
    unfold_spectrum,
    wigner_surmise_pdf,
)

__all__ = [
    "GOE_R",
    "POISSON_R",
    "correlation_matrix",
    "goe_matrix",
    "mean_spacing_ratio",
    "poisson_spectrum",
    "sorted_eigenvalues",
    "spacing_ratios",
    "spectral_form_factor",
    "unfold_spectrum",
    "wigner_surmise_pdf",
]

__version__ = "0.1.0"
