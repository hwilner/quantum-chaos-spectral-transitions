# Contributing

1. Pick an open Issue (they are ordered; check "Depends on").
2. Branch `issue-N-short-name`, commit with messages referencing `#N`.
3. All new statistics need unit tests anchored to a known theoretical value —
   no exceptions; see `tests/test_spectral.py` for the pattern.
4. Docstrings: Google style, with a worked micro-example where possible.
5. Run `pytest` before opening a PR; CI must be green. CI runs the suite on
   NumPy 1.x and 2.x, and checks that the test count stated in the
   documentation still matches the suite.
6. Discuss analysis choices (trim fractions, null variants) in the Issue
   thread *before* implementing — spectral statistics are subtle.
