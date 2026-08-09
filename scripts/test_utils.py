"""Tests for the shared helpers in utils.py.

parse_mutation and iter_path_mutations replaced three hand-rolled parsers of
the same `<parent base><position><new base>` format that disagreed about how
much to validate. These cover the validation the strictest of them had, so the
consolidation cannot quietly lose it.
"""
import os
import tempfile
import unittest

from utils import iter_path_mutations, parse_mutation, read_reference


class ParseMutationTestCase(unittest.TestCase):
    def test_single_digit_position(self):
        self.assertEqual(parse_mutation("A1G"), ("A", 1, "G"))

    def test_multi_digit_position(self):
        """The slice a hand-rolled mut[1:2] would get wrong."""
        self.assertEqual(parse_mutation("C1042T"), ("C", 1042, "T"))

    def test_gap_bases_are_allowed(self):
        self.assertEqual(parse_mutation("-33T"), ("-", 33, "T"))
        self.assertEqual(parse_mutation("A33-"), ("A", 33, "-"))

    def test_iupac_ambiguity_is_allowed(self):
        self.assertEqual(parse_mutation("A33R"), ("A", 33, "R"))

    def test_too_short_raises(self):
        for token in ("", "A", "A1"):
            with self.subTest(token=token):
                with self.assertRaises(ValueError):
                    parse_mutation(token)

    def test_non_numeric_position_raises(self):
        with self.assertRaises(ValueError):
            parse_mutation("AxxG")

    def test_digit_base_raises(self):
        """A digit where a base belongs means the position slice was wrong."""
        with self.assertRaises(ValueError):
            parse_mutation("133G")


class IterPathMutationsTestCase(unittest.TestCase):
    def test_empty_path(self):
        self.assertEqual(list(iter_path_mutations("")), [])

    def test_single_chunk(self):
        self.assertEqual(list(iter_path_mutations("node_1:A1G")),
                         [("A", 1, "G")])

    def test_several_mutations_in_one_chunk(self):
        self.assertEqual(list(iter_path_mutations("node_1:A1G,C3T")),
                         [("A", 1, "G"), ("C", 3, "T")])

    def test_root_to_tip_order_is_preserved(self):
        """Later chunks must come last, so they can supersede earlier ones."""
        got = list(iter_path_mutations("node_1:A1G node_2:G1T"))
        self.assertEqual(got, [("A", 1, "G"), ("G", 1, "T")])

    def test_chunk_with_no_mutations_is_skipped(self):
        self.assertEqual(list(iter_path_mutations("node_1: node_2:A1G")),
                         [("A", 1, "G")])

    def test_node_name_containing_a_colon_raises(self):
        """partition(':') splits once, so the rest is not a mutation token.

        Node names here are EPI_ISL_<n> or internal_<id> and carry no colons,
        so this is unreachable in practice; the point is that it fails loudly
        rather than silently mis-parsing if the naming ever changes.
        """
        with self.assertRaises(ValueError):
            list(iter_path_mutations("node:1:A1G"))

    def test_malformed_mutation_raises(self):
        with self.assertRaises(ValueError):
            list(iter_path_mutations("node_1:nonsense"))


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


if __name__ == "__main__":
    unittest.main()
