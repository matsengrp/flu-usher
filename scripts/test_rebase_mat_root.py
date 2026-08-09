"""Tests for rebase_mat_root.py.

Only the pure helpers are covered here; the protobuf round-trip needs the MAT
schema, which lives in envs/taxonium, and is exercised by the pipeline itself.
"""
import os
import tempfile
import unittest
from types import SimpleNamespace

from rebase_mat_root import (NUCLEOTIDES, read_reference, root_sequence,
                             zero_root_branch_length)


def mutation(position, par_nuc, mut_nuc):
    """A stand-in for the protobuf message, which only needs three fields."""
    return SimpleNamespace(position=position, par_nuc=par_nuc, mut_nuc=[mut_nuc])


A, C, G, T = (NUCLEOTIDES.index(b) for b in "ACGT")


class ReadReferenceTestCase(unittest.TestCase):
    def write(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
        handle.write(text)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_single_record(self):
        self.assertEqual(read_reference(self.write(">r\nACGT\n")), "ACGT")

    def test_wrapped_lines_are_joined(self):
        self.assertEqual(read_reference(self.write(">r\nACGT\nAAGG\n")), "ACGTAAGG")

    def test_stops_at_second_record(self):
        self.assertEqual(read_reference(self.write(">a\nACGT\n>b\nTTTT\n")), "ACGT")

    def test_uppercases(self):
        self.assertEqual(read_reference(self.write(">r\nacgt\n")), "ACGT")


class RootSequenceTestCase(unittest.TestCase):
    """Applying the root's mutation list to the reference."""

    REFERENCE = "ACGTACGTAC"   # A(1) C(2) G(3) T(4) A(5) C(6) G(7) T(8) A(9) C(10)

    def test_no_mutations_returns_reference(self):
        self.assertEqual(root_sequence(self.REFERENCE, []), self.REFERENCE)

    def test_single_mutation(self):
        got = root_sequence(self.REFERENCE, [mutation(1, A, G)])
        self.assertEqual(got, "GCGTACGTAC")

    def test_several_mutations(self):
        got = root_sequence(
            self.REFERENCE,
            [mutation(1, A, G), mutation(4, T, C), mutation(10, C, T)],
        )
        self.assertEqual(got, "GCGCACGTAT")
        self.assertEqual(len(got), len(self.REFERENCE))

    def test_parent_base_mismatch_exits_nonzero(self):
        """par_nuc disagreeing with the reference means they are not a pair."""
        with self.assertRaises(SystemExit) as cm:
            # Position 1 is A, but this mutation claims its parent was T.
            root_sequence(self.REFERENCE, [mutation(1, T, G)])
        self.assertEqual(cm.exception.code, 1)

    def test_position_past_end_exits_nonzero(self):
        with self.assertRaises(SystemExit) as cm:
            root_sequence(self.REFERENCE, [mutation(99, A, G)])
        self.assertEqual(cm.exception.code, 1)

    def test_position_zero_exits_nonzero(self):
        """Positions are 1-based; 0 would silently write to index -1."""
        with self.assertRaises(SystemExit) as cm:
            root_sequence(self.REFERENCE, [mutation(0, A, G)])
        self.assertEqual(cm.exception.code, 1)


class NucleotideEncodingTestCase(unittest.TestCase):
    def test_encoding_matches_usher(self):
        """Verified against matUtils paths: pos 33 par_nuc=0 mut_nuc=3 is A33T."""
        self.assertEqual(NUCLEOTIDES, "ACGT")


class ZeroRootBranchLengthTestCase(unittest.TestCase):
    """Branch lengths are mutation counts, so the emptied root must read 0.

    Caught by comparing against bte, which does this and my first version
    did not: the tree claimed a 113-mutation root branch it no longer held.
    """

    def test_replaces_root_branch_length(self):
        self.assertEqual(
            zero_root_branch_length("(a:1,(b:2,c:3):4):113;"),
            "(a:1,(b:2,c:3):4):0;",
        )

    def test_leaves_inner_branch_lengths_alone(self):
        got = zero_root_branch_length("(a:1,(b:2,c:3):4):113;")
        self.assertIn("a:1", got)
        self.assertIn("b:2", got)
        self.assertIn(":4)", got)

    def test_named_root(self):
        self.assertEqual(
            zero_root_branch_length("(a:1,b:2)node_1:9;"),
            "(a:1,b:2)node_1:0;",
        )

    def test_no_root_branch_length_is_left_alone(self):
        self.assertEqual(zero_root_branch_length("(a:1,b:2);"), "(a:1,b:2);")

    def test_missing_terminator(self):
        self.assertEqual(zero_root_branch_length("(a:1,b:2):7"), "(a:1,b:2):0")

    def test_non_newick_is_returned_unchanged(self):
        self.assertEqual(zero_root_branch_length("nonsense"), "nonsense")


if __name__ == "__main__":
    unittest.main()
