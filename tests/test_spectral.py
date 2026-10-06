"""Unit tests for the spectral_chaos package."""

import unittest

import numpy as np

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


class TestReferenceEnsembles(unittest.TestCase):
    """Spacing statistics on analytic reference ensembles."""

    def test_poisson_ratio_constant(self):
        self.assertAlmostEqual(POISSON_R, 2 * np.log(2) - 1, places=12)
        self.assertAlmostEqual(POISSON_R, 0.38629, places=4)

    def test_poisson_spectrum_recovers_ratio(self):
        eigs = poisson_spectrum(6000, seed=0)
        self.assertAlmostEqual(mean_spacing_ratio(eigs), POISSON_R,
                               delta=0.02)

    def test_goe_recovers_ratio(self):
        # ratios are pooled per realisation: concatenating sorted spectra
        # from independent matrices would interleave levels and mimic
        # Poisson statistics (a classic pitfall).
        rs = np.concatenate([spacing_ratios(sorted_eigenvalues(
            goe_matrix(60, seed=s))[10:-10]) for s in range(40)])
        self.assertAlmostEqual(float(rs.mean()), GOE_R, delta=0.02)

    def test_goe_and_poisson_are_discriminated(self):
        goe_r = mean_spacing_ratio(sorted_eigenvalues(goe_matrix(200, 1)))
        pois_r = mean_spacing_ratio(poisson_spectrum(200, 1))
        self.assertGreater(goe_r - pois_r, 0.08)


class TestSpacingMechanics(unittest.TestCase):
    """Basic properties of the estimators."""

    def test_ratios_in_unit_interval(self):
        r = spacing_ratios(poisson_spectrum(500, seed=2))
        self.assertTrue(np.all(r >= 0.0))
        self.assertTrue(np.all(r <= 1.0))

    def test_equally_spaced_gives_one(self):
        eigs = np.linspace(0, 10, 50)
        np.testing.assert_allclose(spacing_ratios(eigs), 1.0)

    def test_unfold_unit_density(self):
        eigs = np.sort(np.random.default_rng(3).normal(size=200))
        u = unfold_spectrum(eigs, poly_degree=5)
        self.assertEqual(u.shape, eigs.shape)
        # mean density of unfolded levels ~ 1: spacings average to ~1
        self.assertAlmostEqual(np.diff(u).mean(), 1.0, delta=0.15)

    def test_too_few_levels_raise(self):
        with self.assertRaises(ValueError):
            unfold_spectrum(np.array([0.1, 0.2]), poly_degree=3)
        with self.assertRaises(ValueError):
            spacing_ratios(np.array([0.0]))


class TestFormFactor(unittest.TestCase):
    """Spectral form factor behaviour."""

    def test_poisson_form_factor_is_flat(self):
        eigs = poisson_spectrum(4000, seed=4)
        u = unfold_spectrum(eigs, poly_degree=3)
        taus = np.linspace(0.2, 3.0, 20)
        k = spectral_form_factor(u, taus)
        # Poisson: K(tau) / 1 -> 1 away from tau = 0 (plus O(1) noise)
        self.assertAlmostEqual(float(np.median(k)), 1.0, delta=0.35)

    def test_form_factor_at_zero(self):
        eigs = poisson_spectrum(500, seed=5)
        k0 = spectral_form_factor(eigs, np.array([0.0]))
        self.assertAlmostEqual(k0[0], 500.0, places=6)  # N^2 / N


class TestWignerSurmise(unittest.TestCase):
    """Reference spacing distribution."""

    def test_normalisation_and_mean(self):
        s = np.linspace(0, 4, 4001)
        p = wigner_surmise_pdf(s, beta=1)
        ds = s[1] - s[0]
        self.assertAlmostEqual(float(np.trapezoid(p, s)), 1.0, delta=1e-3)
        self.assertAlmostEqual(float(np.trapezoid(s * p, s)), 1.0, delta=1e-3)

    def test_level_repulsion_at_zero(self):
        self.assertEqual(wigner_surmise_pdf(np.array([0.0]))[0], 0.0)


class TestCorrelationPipeline(unittest.TestCase):
    """fMRI-style entry point."""

    def test_correlation_matrix_shape_and_symmetry(self):
        rng = np.random.default_rng(7)
        x = rng.normal(size=(300, 25))
        c = correlation_matrix(x)
        self.assertEqual(c.shape, (25, 25))
        np.testing.assert_allclose(c, c.T, atol=1e-12)
        np.testing.assert_allclose(np.diag(c), 1.0, atol=1e-10)

    def test_independent_levels_are_poisson_like(self):
        # a diagonal matrix with i.i.d. entries = independent levels
        rng = np.random.default_rng(8)
        d = np.diag(rng.normal(size=400))
        r = mean_spacing_ratio(sorted_eigenvalues(d))
        self.assertAlmostEqual(r, POISSON_R, delta=0.04)

    def test_white_noise_covariance_repels_levels(self):
        # Wishart sample covariance of white noise is in the same
        # universality class as GOE in the bulk: it shows level repulsion
        # even though nothing is "chaotic".  This is precisely why the
        # project needs matrix-ensemble nulls (see docs/METHODS.md).
        # pooled over realisations to average out single-matrix noise
        rs = []
        for seed in range(10):
            rng = np.random.default_rng(seed)
            x = rng.normal(size=(400, 80))
            w = x.T @ x / 400.0
            rs.append(spacing_ratios(sorted_eigenvalues(w)[10:-10]))
        r = float(np.concatenate(rs).mean())
        self.assertGreater(r, POISSON_R + 0.1)


if __name__ == "__main__":
    unittest.main()
