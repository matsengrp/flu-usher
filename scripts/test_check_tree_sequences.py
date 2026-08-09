"""Tests for check_tree_sequences.py.

The end-to-end cases run the script as a subprocess so that exit codes are
tested as the pipeline sees them: a gate that logs an error but exits 0 is
worse than no gate, because it reports success.
"""
import os
import subprocess
import sys
import tempfile
import unittest

from check_tree_sequences import compose, root_children

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "check_tree_sequences.py")

# A(1) C(2) G(3) T(4) A(5) C(6) G(7) T(8) A(9) C(10)
REFERENCE = ">ref\nACGTACGTAC\n"

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEPI_A\tEPI_B\n"
)


class ComposeTestCase(unittest.TestCase):
    """The root-to-tip path collapse."""

    def setUp(self):
        self.ref = {i: b for i, b in enumerate("ACGTACGTAC", start=1)}

    def test_single_mutation(self):
        self.assertEqual(compose("node_1:A1G", self.ref), {1: "G"})

    def test_back_mutation_cancels(self):
        """A1G then G1A returns to reference, so the site must drop out."""
        self.assertEqual(compose("node_1:A1G node_2:G1A", self.ref), {})

    def test_later_chunk_overrides_earlier(self):
        self.assertEqual(compose("node_1:A1G node_2:G1C", self.ref), {1: "C"})

    def test_mutation_matching_reference_is_dropped(self):
        """A path may assert the reference base; that is not a difference."""
        self.assertEqual(compose("node_1:G1A", self.ref), {})

    def test_multiple_mutations_in_one_chunk(self):
        self.assertEqual(compose("node_1:A1G,G3T", self.ref), {1: "G", 3: "T"})

    def test_empty_path(self):
        self.assertEqual(compose("", self.ref), {})

    def test_chunk_with_no_mutations(self):
        self.assertEqual(compose("node_1: node_2:A1G", self.ref), {1: "G"})


class RootChildrenTestCase(unittest.TestCase):
    """Reading the root's direct children out of a newick string."""

    def test_leaf_and_internal_children(self):
        self.assertEqual(
            root_children("(EPI_A:55,(EPI_B:1,EPI_C:1)node_2:3);"),
            ["EPI_A", "node_2"],
        )

    def test_polytomy(self):
        self.assertEqual(
            root_children("(a:1,b:1,c:1,d:1);"), ["a", "b", "c", "d"]
        )

    def test_nested_commas_do_not_split(self):
        self.assertEqual(
            root_children("(((p,q),(r,s))node_i:2,t:1);"), ["node_i", "t"]
        )

    def test_no_branch_lengths(self):
        self.assertEqual(root_children("(a,(b,c)n);"), ["a", "n"])

    def test_rejects_non_newick(self):
        with self.assertRaises(ValueError):
            root_children("not a tree")


class EndToEndTestCase(unittest.TestCase):
    """Exit codes, as Snakemake sees them."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name, text):
        path = os.path.join(self.tmp.name, name)
        with open(path, "w") as handle:
            handle.write(text)
        return path

    def run_check(self, sampled, final, vcf_rows="", newick=None,
                  expect_root=None):
        args = [
            sys.executable, SCRIPT,
            "--sampled-paths", self.write("sampled.tsv", sampled),
            "--final-paths", self.write("final.tsv", final),
            "--reference", self.write("ref.fasta", REFERENCE),
            "--input-vcf", self.write("in.vcf", VCF_HEADER + vcf_rows),
            "--output", os.path.join(self.tmp.name, "report.txt"),
        ]
        if newick is not None:
            args += ["--final-newick", self.write("final.nwk", newick)]
        if expect_root is not None:
            args += ["--expect-root", expect_root]
        return subprocess.run(args, capture_output=True, text=True,
                              cwd=os.path.dirname(SCRIPT))

    def test_identical_sequences_pass(self):
        sampled = "EPI_A\tnode_1:A1G\nEPI_B\tnode_1:A1G node_3:G3T\n"
        final = "EPI_A\tnode_9:A1G\nEPI_B\tnode_9:A1G node_8:G3T\n"
        result = self.run_check(sampled, final)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_differing_called_site_fails(self):
        """The negative control: a real genotype change must exit nonzero."""
        sampled = "EPI_A\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:A1C\n"
        result = self.run_check(sampled, final)
        self.assertEqual(result.returncode, 1)
        self.assertIn("EPI_A", result.stderr)

    def test_empty_paths_do_not_pass(self):
        """An empty matUtils -S file compared nothing; that is not success."""
        result = self.run_check("", "")
        self.assertEqual(result.returncode, 1)
        self.assertIn("compared nothing", result.stderr)

    def test_difference_at_uncalled_site_is_masked(self):
        """Gaps are missing data, so the imputed base may move with the root."""
        sampled = "EPI_A\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:\n"
        vcf_rows = "chr\t1\t.\tA\tG\t.\t.\t.\tGT\t.\t1\n"
        result = self.run_check(sampled, final, vcf_rows=vcf_rows)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_sample_set_mismatch_fails(self):
        sampled = "EPI_A\tnode_1:A1G\nEPI_B\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:A1G\n"
        result = self.run_check(sampled, final)
        self.assertEqual(result.returncode, 1)
        self.assertIn("sample sets differ", result.stderr)

    def test_wrong_root_fails(self):
        sampled = "EPI_A\tnode_1:A1G\n"
        result = self.run_check(
            sampled, sampled,
            newick="(EPI_A:1,(EPI_B:1,EPI_C:1)node_2:2);",
            expect_root="EPI_B",
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("not rooted at", result.stderr)

    def test_expected_root_passes(self):
        sampled = "EPI_A\tnode_1:A1G\n"
        result = self.run_check(
            sampled, sampled,
            newick="(EPI_A:1,(EPI_B:1,EPI_C:1)node_2:2);",
            expect_root="EPI_A",
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_expect_root_without_newick_fails(self):
        sampled = "EPI_A\tnode_1:A1G\n"
        result = self.run_check(sampled, sampled, expect_root="EPI_A")
        self.assertEqual(result.returncode, 1)


if __name__ == "__main__":
    unittest.main()
