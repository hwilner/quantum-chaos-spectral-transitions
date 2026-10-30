"""Check that the documentation's test count matches the test suite.

The project states its number of unit tests in several files.  That number
drifts whenever a test is added, and a stale count in a results document is
a published claim that the code does not support, so it is checked here.

Run it as ``python tools/check_doc_test_count.py``; it exits non-zero when a
documented count disagrees with the number of tests actually collected.
"""

from __future__ import annotations

import pathlib
import re
import sys
import unittest

#: Repository root, derived from the location of this file.
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]

#: Directory holding the package, added to the import path so that the
#: checker runs from a fresh checkout without an editable install.
SOURCE_DIR = REPO_ROOT / "src"

#: Documents that state a test count, with the pattern that extracts it.
#: The pattern must have exactly one capturing group holding the count.
DOCUMENTS = {
    "README.md": r"(\d+) unit tests",
    "docs/CURRENT_RESULTS_AND_DISCUSSION.md": r"All (\d+) unit tests pass",
    "docs/STATUS_AND_PLAN.md": r"(\d+)-test suite",
    "docs/METHODS.md": r"\((\d+) anchors\)",
}


def collect_test_count() -> int:
    """Return the number of tests collected from the ``tests`` directory.

    Returns:
        The number of test cases that the test suite contains.
    """
    if str(SOURCE_DIR) not in sys.path:
        sys.path.insert(0, str(SOURCE_DIR))
    suite = unittest.TestLoader().discover(
        str(REPO_ROOT / "tests"), top_level_dir=str(REPO_ROOT / "tests"))
    return suite.countTestCases()


def documented_counts() -> dict[str, list[int]]:
    """Return every test count stated in the documentation.

    Returns:
        A mapping from file name to the list of counts found in that file.
    """
    found: dict[str, list[int]] = {}
    for relative_path, pattern in DOCUMENTS.items():
        text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
        found[relative_path] = [int(value) for value in re.findall(pattern, text)]
    return found


def main() -> int:
    """Compare the documented test counts against the collected test count.

    Returns:
        ``0`` when every documented count matches the suite, ``1`` otherwise.
    """
    actual = collect_test_count()
    mismatches = []
    for relative_path, counts in documented_counts().items():
        if not counts:
            mismatches.append(
                f"{relative_path}: no test count found (expected a match for "
                f"{DOCUMENTS[relative_path]!r})")
            continue
        for count in counts:
            if count != actual:
                mismatches.append(
                    f"{relative_path}: states {count} tests, suite has {actual}")
    if mismatches:
        print("documented test count does not match the test suite:")
        for line in mismatches:
            print(f"  {line}")
        return 1
    print(f"documented test count matches the suite ({actual} tests)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
