"""Unit tests for the input validation of the spectral_chaos package.

Each estimator in this package must reject an input it cannot handle rather
than return a well-formed number that describes a different calculation than
the one the caller asked for.  These tests pin that contract; the statistical
anchors live in ``test_spectral.py``.
"""

import unittest

import numpy as np

from spectral_chaos.spectral import (
    correlation_matrix,
    sorted_eigenvalues,
)


class TestCorrelationMatrixValidation(unittest.TestCase):
    """The correlation-matrix entry point must not return NaN."""

    def test_constant_channel_is_rejected(self):
        # A dead parcel has zero variance; np.corrcoef would put NaN in the
        # whole row and column, and the failure would surface in the
        # eigensolver rather than here.
        x = np.random.default_rng(0).normal(size=(200, 4))
        x[:, 2] = 3.0
        with self.assertRaises(ValueError) as ctx:
            correlation_matrix(x)
        self.assertIn("2", str(ctx.exception))

    def test_several_constant_channels_are_all_named(self):
        x = np.random.default_rng(1).normal(size=(200, 4))
        x[:, 0] = 1.0
        x[:, 3] = -2.0
        with self.assertRaises(ValueError) as ctx:
            correlation_matrix(x)
        message = str(ctx.exception)
        self.assertIn("0", message)
        self.assertIn("3", message)

    def test_non_finite_input_is_rejected(self):
        x = np.random.default_rng(2).normal(size=(200, 3))
        x[5, 1] = np.nan
        with self.assertRaises(ValueError):
            correlation_matrix(x)

    def test_infinite_input_is_rejected(self):
        x = np.random.default_rng(3).normal(size=(200, 3))
        x[5, 1] = np.inf
        with self.assertRaises(ValueError):
            correlation_matrix(x)

    def test_single_channel_returns_a_matrix(self):
        # np.corrcoef collapses one channel to the scalar 1.0, but the
        # documented contract is a (n_channels, n_channels) matrix.
        c = correlation_matrix(np.random.default_rng(4).normal(size=(50, 1)))
        self.assertEqual(c.shape, (1, 1))
        np.testing.assert_allclose(c, [[1.0]], atol=1e-12)

    def test_single_timepoint_is_rejected(self):
        # A correlation needs at least two observations per channel.
        with self.assertRaises(ValueError):
            correlation_matrix(np.random.default_rng(5).normal(size=(1, 4)))

    def test_no_channel_is_rejected(self):
        with self.assertRaises(ValueError):
            correlation_matrix(np.zeros((10, 0)))

    def test_one_dimensional_input_is_rejected(self):
        with self.assertRaises(ValueError):
            correlation_matrix(np.random.default_rng(6).normal(size=10))

    def test_valid_input_is_unaffected(self):
        x = np.random.default_rng(7).normal(size=(300, 25))
        c = correlation_matrix(x)
        self.assertEqual(c.shape, (25, 25))
        np.testing.assert_allclose(c, c.T, atol=1e-12)


class TestSortedEigenvaluesValidation(unittest.TestCase):
    """The eigensolver entry point must not read half of a matrix."""

    def test_asymmetric_matrix_is_rejected(self):
        # numpy.linalg.eigvalsh reads only the lower triangle and never
        # reports asymmetry, so this would silently return the eigenvalues of
        # a different matrix.
        rng = np.random.default_rng(8)
        a = rng.normal(size=(4, 4))
        a = a + a.T
        a[0, 1] += 50.0
        with self.assertRaises(ValueError):
            sorted_eigenvalues(a)

    def test_non_finite_matrix_is_rejected(self):
        with self.assertRaises(ValueError):
            sorted_eigenvalues(np.array([[1.0, np.nan], [np.nan, 1.0]]))

    def test_non_square_matrix_is_rejected(self):
        with self.assertRaises(ValueError):
            sorted_eigenvalues(np.zeros((3, 4)))

    def test_symmetric_matrix_is_accepted(self):
        rng = np.random.default_rng(9)
        a = rng.normal(size=(20, 20))
        a = a + a.T
        eigs = sorted_eigenvalues(a)
        np.testing.assert_allclose(eigs, np.linalg.eigvalsh(a), atol=1e-10)

    def test_matrix_symmetric_to_roundoff_is_accepted(self):
        # The correlation matrix that the package itself produces is
        # symmetric only to roundoff, and must still be accepted.
        x = np.random.default_rng(10).normal(size=(400, 30))
        c = correlation_matrix(x)
        self.assertFalse(np.array_equal(c, c.T))
        eigs = sorted_eigenvalues(c)
        self.assertEqual(eigs.size, 30)


if __name__ == "__main__":
    unittest.main()
