"""
Tests for create_samples_file.py.

The behaviour under test is that a filter matching nothing is an error. Writing
a root-only samples file instead makes matUtils extract and taxonium produce a
single-node "subtree" that Snakemake records as a successful build, so the
mistake surfaces only when someone opens the visualization.
"""

import lzma
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))

import create_samples_file


ROOT_ID = "EPI_ISL_ROOT"
REFERENCE_ID = "REFERENCE"


def write_inputs(tmpdir, rows):
    """
    Build the three inputs the script reads.

    rows: list of (isolate_id, geographic_group) for the metadata. Every
    isolate_id also goes into the curated MSA, behind the reference.
    """
    msa = os.path.join(tmpdir, "curated_msa.fasta.xz")
    with lzma.open(msa, "wt") as f:
        f.write(f">{REFERENCE_ID}\nACGT\n")
        for isolate_id, _ in rows:
            f.write(f">{isolate_id}\nACGT\n")

    root = os.path.join(tmpdir, "curated_root.fasta")
    with open(root, "w") as f:
        f.write(f">{ROOT_ID}\nACGT\n")

    metadata = os.path.join(tmpdir, "metadata.csv")
    with open(metadata, "w") as f:
        f.write("isolate_id,geographic_group\n")
        f.write(f"{ROOT_ID},north_america\n")
        for isolate_id, group in rows:
            f.write(f"{isolate_id},{group}\n")

    return msa, root, metadata


def run_script(tmpdir, rows, value, column="geographic_group"):
    """Invoke main() as the rule does. Returns the output path."""
    msa, root, metadata = write_inputs(tmpdir, rows)
    output = os.path.join(tmpdir, "samples.txt")
    argv = [
        "create_samples_file.py",
        "--curated-msa", msa,
        "--metadata", metadata,
        "--column", column,
        "--value", value,
        "--root", root,
        "--output", output,
    ]
    with mock.patch.object(sys, "argv", argv):
        create_samples_file.main()
    return output


class TestZeroMatchesIsAnError(unittest.TestCase):
    """A filter matching no samples must fail, not emit a root-only file."""

    def test_no_matching_samples_exits_nonzero(self):
        rows = [("EPI_ISL_1", "europe"), ("EPI_ISL_2", "europe")]
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(SystemExit) as cm:
                run_script(tmpdir, rows, value="asia")
            self.assertEqual(cm.exception.code, 1)

    def test_no_matching_samples_does_not_write_root_only_file(self):
        """The regression: a 1-line file that renders as a single-node tree.

        Asserted unconditionally. Guarding this on `if os.path.exists(output)`
        made it vacuous -- the script exits before opening the output, so the
        branch never ran and this duplicated the test above. The point is to
        catch a future reordering that writes first and validates after.
        """
        rows = [("EPI_ISL_1", "europe")]
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "samples.txt")
            with self.assertRaises(SystemExit):
                run_script(tmpdir, rows, value="asia")
            self.assertFalse(
                os.path.exists(output),
                "wrote a samples file instead of failing outright",
            )

    def test_misspelled_value_exits_nonzero(self):
        """The realistic trigger: a typo in a configured geographic group."""
        rows = [("EPI_ISL_1", "north_america")]
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(SystemExit) as cm:
                run_script(tmpdir, rows, value="north_americaa")
            self.assertEqual(cm.exception.code, 1)


class TestMatchingSamplesStillSucceed(unittest.TestCase):
    """The failure path must not have broken the normal path."""

    def test_matching_samples_written_with_root_first(self):
        rows = [
            ("EPI_ISL_3", "asia"),
            ("EPI_ISL_1", "asia"),
            ("EPI_ISL_2", "europe"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            output = run_script(tmpdir, rows, value="asia")
            with open(output) as f:
                lines = f.read().splitlines()
        # Root first, then the matching samples sorted for reproducibility.
        self.assertEqual(lines, [ROOT_ID, "EPI_ISL_1", "EPI_ISL_3"])

    def test_single_match_succeeds(self):
        """One match is enough; only zero is an error."""
        rows = [("EPI_ISL_1", "asia"), ("EPI_ISL_2", "europe")]
        with tempfile.TemporaryDirectory() as tmpdir:
            output = run_script(tmpdir, rows, value="asia")
            with open(output) as f:
                lines = f.read().splitlines()
        self.assertEqual(lines, [ROOT_ID, "EPI_ISL_1"])

    def test_missing_column_still_exits_nonzero(self):
        """Pre-existing validation must be unaffected."""
        rows = [("EPI_ISL_1", "asia")]
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(SystemExit) as cm:
                run_script(tmpdir, rows, value="asia", column="host_group")
            self.assertEqual(cm.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
