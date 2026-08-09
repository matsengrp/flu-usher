"""Tests for extract_root_sequence.py.

The module had no tests at all while carrying the fail-fast branches that
decide whether curated_root.fasta is trustworthy: a mutation path that does not
match the reference it is applied to produces a wrong root sequence, and the
root sequence is what every downstream tree is compared against.

These run main() in-process, as the rule does, rather than reconstructing its
argument handling.
"""
import lzma
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))

import extract_root_sequence
from extract_root_sequence import apply_mutations, parse_mutation_path

ROOT = "EPI_ISL_ROOT"
# A(1) C(2) G(3) T(4) A(5) C(6) G(7) T(8) A(9) C(10)
REFERENCE = "ACGTACGTAC"


class ParseMutationPathTestCase(unittest.TestCase):
    def write(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
        handle.write(text)
        handle.close()
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_returns_parsed_triples_in_root_to_tip_order(self):
        path = self.write(f"{ROOT}\tnode_1:A1G node_2:C2T,G3A\n")
        self.assertEqual(
            parse_mutation_path(path, ROOT),
            [("A", 1, "G"), ("C", 2, "T"), ("G", 3, "A")],
        )

    def test_empty_path_is_no_mutations(self):
        self.assertEqual(parse_mutation_path(self.write(f"{ROOT}\t\n"), ROOT), [])

    def test_missing_sample_raises(self):
        with self.assertRaises(ValueError):
            parse_mutation_path(self.write("OTHER\tnode_1:A1G\n"), ROOT)

    def test_malformed_token_raises_rather_than_mis_slicing(self):
        """The shared parser's job; this asserts it is actually reached."""
        with self.assertRaises(ValueError):
            parse_mutation_path(self.write(f"{ROOT}\tnode_1:nonsense\n"), ROOT)


class ApplyMutationsTestCase(unittest.TestCase):
    def test_applies_in_order(self):
        got = apply_mutations(REFERENCE, [("A", 1, "G"), ("T", 4, "C")])
        self.assertEqual(got, "GCGCACGTAC")

    def test_same_position_twice_uses_the_running_state(self):
        """The second mutation's parent is the first one's result, not the ref."""
        got = apply_mutations(REFERENCE, [("A", 1, "G"), ("G", 1, "T")])
        self.assertEqual(got, "TCGTACGTAC")

    def test_parent_base_mismatch_exits_nonzero(self):
        with self.assertRaises(SystemExit) as cm:
            apply_mutations(REFERENCE, [("T", 1, "G")])   # position 1 is A
        self.assertEqual(cm.exception.code, 1)

    def test_position_past_end_exits_nonzero(self):
        with self.assertRaises(SystemExit) as cm:
            apply_mutations(REFERENCE, [("A", 99, "G")])
        self.assertEqual(cm.exception.code, 1)

    def test_position_zero_exits_nonzero(self):
        """1-based; 0 would write to index -1, silently hitting the last base."""
        with self.assertRaises(SystemExit) as cm:
            apply_mutations(REFERENCE, [("A", 0, "G")])
        self.assertEqual(cm.exception.code, 1)


class MainTestCase(unittest.TestCase):
    """End to end, as the create_root_fasta rule invokes it."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def build(self, path_field, msa_seq=None):
        d = self.tmp.name
        ref = os.path.join(d, "ref.fasta")
        with open(ref, "w") as f:
            f.write(f">reference\n{REFERENCE}\n")
        msa = os.path.join(d, "msa.fasta.xz")
        with lzma.open(msa, "wt") as f:
            f.write(f">{ROOT}\n{msa_seq if msa_seq else REFERENCE}\n")
        paths = os.path.join(d, "paths.txt")
        with open(paths, "w") as f:
            f.write(f"{ROOT}\t{path_field}\n")
        out = os.path.join(d, "root.fasta")
        return ["extract_root_sequence.py", "--reference", ref, "--msa", msa,
                "--paths", paths, "--new-root-name", ROOT, "--output", out], out

    def test_no_mutations_reproduces_the_reference(self):
        argv, out = self.build("")
        with mock.patch.object(sys, "argv", argv):
            extract_root_sequence.main()
        with open(out) as f:
            self.assertEqual(f.read().split("\n", 1)[1].replace("\n", ""), REFERENCE)

    def test_malformed_token_exits_nonzero_without_a_traceback(self):
        """A ValueError from the parser must become a clean exit, like the rest."""
        argv, _ = self.build("node_1:nonsense")
        with mock.patch.object(sys, "argv", argv):
            with self.assertRaises(SystemExit) as cm:
                extract_root_sequence.main()
        self.assertEqual(cm.exception.code, 1)

    def test_missing_sample_exits_nonzero(self):
        argv, _ = self.build("")
        argv[argv.index("--new-root-name") + 1] = "NOT_PRESENT"
        with mock.patch.object(sys, "argv", argv):
            with self.assertRaises(SystemExit) as cm:
                extract_root_sequence.main()
        self.assertEqual(cm.exception.code, 1)

    def test_mismatch_against_the_msa_is_fatal(self):
        """validate_sequences guards the result even when the path applied."""
        argv, _ = self.build("node_1:A1G", msa_seq="ACGTACGTAC")
        with mock.patch.object(sys, "argv", argv):
            with self.assertRaises(ValueError):
                extract_root_sequence.main()


if __name__ == "__main__":
    unittest.main()
