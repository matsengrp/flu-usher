"""Tests for rebase_mat_root.py.

Covers the helpers and the rebase itself, via SimpleNamespace stand-ins for the
mutation messages -- rebase_onto_root() takes anything with .mutation lists, so
it needs neither the compiled schema (which lives in envs/taxonium) nor a real
protobuf. What is still uncovered here is the file-level round trip: load_mat,
save_mat, and that the fields this script never touches (condensed_nodes,
metadata) survive serialization. Those are exercised only by the pipeline.
"""
import os
import tempfile
import unittest
from types import SimpleNamespace

from rebase_mat_root import (NUCLEOTIDES, read_reference, rebase_onto_root,
                             root_branch_length, root_sequence,
                             zero_root_branch_length)


def mutation(position, par_nuc, *mut_nuc):
    """A stand-in for the protobuf message, which only needs three fields.

    mut_nuc is variadic because the schema declares it repeated, so it can
    carry an ambiguous state rather than a single base.
    """
    return SimpleNamespace(position=position, par_nuc=par_nuc, mut_nuc=list(mut_nuc))


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

    def test_ambiguous_mut_nuc_exits_nonzero(self):
        """mut_nuc is repeated; taking [0] would silently pick one base."""
        with self.assertRaises(SystemExit) as cm:
            root_sequence(self.REFERENCE, [mutation(1, A, G, T)])
        self.assertEqual(cm.exception.code, 1)

    def test_empty_mut_nuc_exits_nonzero(self):
        """An empty repeated field would raise IndexError rather than explain."""
        with self.assertRaises(SystemExit) as cm:
            root_sequence(self.REFERENCE, [mutation(1, A)])
        self.assertEqual(cm.exception.code, 1)


class RebaseOntoRootTestCase(unittest.TestCase):
    """The bookkeeping itself: empty the root, repoint ref_nuc, touch nothing else.

    The operation this module exists to perform, and the one a production run
    was previously the first thing to exercise. Uses the same SimpleNamespace
    stand-ins as the rest of the file, so it runs without the compiled schema.
    """

    def node(self, *mutations):
        return SimpleNamespace(mutation=list(mutations))

    def with_ref(self, position, ref_nuc, par_nuc, mut_nuc):
        m = mutation(position, par_nuc, mut_nuc)
        m.ref_nuc = ref_nuc
        return m

    def setUp(self):
        # Root moves the origin at positions 1 (A->G) and 4 (T->C).
        self.root = self.node(self.with_ref(1, A, A, G), self.with_ref(4, T, T, C))
        # A descendant mutating at position 1, which the origin moved...
        self.child = self.node(self.with_ref(1, A, G, T))
        # ...and one at position 7, which it did not.
        self.other = self.node(self.with_ref(7, G, G, A))
        self.nodes = [self.root, self.child, self.other]

    def test_root_mutation_list_is_emptied(self):
        rebase_onto_root(self.nodes)
        self.assertEqual(list(self.root.mutation), [])

    def test_reports_the_positions_that_moved(self):
        moved, _ = rebase_onto_root(self.nodes)
        self.assertEqual(moved, {1: G, 4: C})

    def test_repoints_ref_nuc_at_moved_positions(self):
        rebase_onto_root(self.nodes)
        self.assertEqual(self.child.mutation[0].ref_nuc, G)

    def test_leaves_ref_nuc_alone_elsewhere(self):
        rebase_onto_root(self.nodes)
        self.assertEqual(self.other.mutation[0].ref_nuc, G)  # unchanged

    def test_does_not_touch_par_nuc_or_mut_nuc(self):
        """The claim that this is bookkeeping: no parent-to-child change moves."""
        before = [(m.par_nuc, list(m.mut_nuc))
                  for n in (self.child, self.other) for m in n.mutation]
        rebase_onto_root(self.nodes)
        after = [(m.par_nuc, list(m.mut_nuc))
                 for n in (self.child, self.other) for m in n.mutation]
        self.assertEqual(before, after)

    def test_counts_every_repointed_annotation(self):
        """Counts records, not positions: two nodes mutating at position 1.

        Includes the root's own two, which are repointed and then deleted with
        the rest of its list. So the logged count is of annotations touched,
        not of annotations surviving in the output.
        """
        second = self.node(self.with_ref(1, A, G, C))
        _, repointed = rebase_onto_root([self.root, self.child, second, self.other])
        self.assertEqual(repointed, 4)  # root's pos 1 and 4, plus two children

    def test_root_with_no_mutations_is_a_no_op(self):
        nodes = [self.node(), self.child]
        moved, repointed = rebase_onto_root(nodes)
        self.assertEqual((moved, repointed), ({}, 0))
        self.assertEqual(self.child.mutation[0].ref_nuc, A)

    def test_is_idempotent(self):
        """Rebasing an already-rebased tree changes nothing further."""
        rebase_onto_root(self.nodes)
        snapshot = [(m.position, m.ref_nuc, m.par_nuc) for n in self.nodes
                    for m in n.mutation]
        moved, repointed = rebase_onto_root(self.nodes)
        self.assertEqual((moved, repointed), ({}, 0))
        self.assertEqual([(m.position, m.ref_nuc, m.par_nuc) for n in self.nodes
                          for m in n.mutation], snapshot)


class RootBranchLengthTestCase(unittest.TestCase):
    """The newick's independent statement of the root's mutation count."""

    def test_reads_the_root_branch_length(self):
        self.assertEqual(root_branch_length("(a:1,(b:2,c:3):4):113;"), 113)

    def test_zero(self):
        self.assertEqual(root_branch_length("(a:1,b:2):0;"), 0)

    def test_named_root(self):
        self.assertEqual(root_branch_length("(a:1,b:2)node_1:9;"), 9)

    def test_absent_branch_length_is_none(self):
        self.assertIsNone(root_branch_length("(a:1,b:2);"))

    def test_non_integer_is_none(self):
        """Only usher's mutation counts are meaningful; floats are not ours."""
        self.assertIsNone(root_branch_length("(a:1,b:2):0.5;"))

    def test_non_newick_is_none(self):
        self.assertIsNone(root_branch_length("nonsense"))


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
