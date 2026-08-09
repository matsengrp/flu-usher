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

    def test_multi_digit_position(self):
        """Pins the mut[1:-1] slice, which single-digit fixtures can't."""
        self.assertEqual(compose("node_1:C10G", self.ref), {10: "G"})

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
                  expect_root=None, final_origin=None):
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
        if final_origin is not None:
            args += ["--final-origin", self.write("origin.fasta", final_origin)]
        return subprocess.run(args, capture_output=True, text=True,
                              cwd=os.path.dirname(SCRIPT))

    def report(self):
        """The report file's counters, parsed into a dict.

        Exit code alone does not say the gate measured the right thing: a
        script that always exited 0 while writing nonsense would satisfy every
        passing case here. These assert on what it actually counted.
        """
        with open(os.path.join(self.tmp.name, "report.txt")) as handle:
            counters = {}
            for line in handle:
                key, _, value = line.partition(":")
                if value.strip().isdigit():
                    counters[key.strip()] = int(value.strip())
        return counters

    def test_identical_sequences_pass(self):
        sampled = "EPI_A\tnode_1:A1G\nEPI_B\tnode_1:A1G node_3:G3T\n"
        final = "EPI_A\tnode_9:A1G\nEPI_B\tnode_9:A1G node_8:G3T\n"
        result = self.run_check(sampled, final)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.report(), {
            "samples compared": 2,
            "identical sequences": 2,
            "differing sequences": 0,
            "samples differing only at uncalled sites": 0,
            "uncalled sites excluded": 0,
        })

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
        # Masked, not silently ignored: the report must say a site was excluded.
        self.assertEqual(self.report()["uncalled sites excluded"], 1)
        self.assertEqual(self.report()["samples differing only at uncalled sites"], 1)

    def test_called_differing_site_in_the_vcf_still_fails(self):
        """Separates 'no VCF row' from 'the VCF says this base was called'."""
        sampled = "EPI_A\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:A1C\n"
        vcf_rows = "chr\t1\t.\tA\tG\t.\t.\t.\tGT\t1\t1\n"
        result = self.run_check(sampled, final, vcf_rows=vcf_rows)
        self.assertEqual(result.returncode, 1)

    def test_sample_set_mismatch_fails(self):
        sampled = "EPI_A\tnode_1:A1G\nEPI_B\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:A1G\n"
        result = self.run_check(sampled, final)
        self.assertEqual(result.returncode, 1)
        self.assertIn("sample sets differ", result.stderr)

    def test_report_lists_both_kinds_of_differing_site(self):
        """A key-set difference must not hide a same-key value difference.

        EPI_A differs at position 1 (G vs C, called in both) and at position 3
        (present only in the final tree). Reporting only the latter would send
        someone debugging this to the wrong site.
        """
        sampled = "EPI_A\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:A1C,G3T\n"
        result = self.run_check(sampled, final)
        self.assertEqual(result.returncode, 1)
        # Assert the exact list from the report, not a substring of stderr:
        # log lines carry a timestamp, and "14:30:16,636" contains both a "1"
        # and a "3", so substring checks pass no matter what was reported.
        # The old expression yields [3] here; this must distinguish them.
        with open(os.path.join(self.tmp.name, "report.txt")) as handle:
            report = handle.read()
        self.assertIn("EPI_A: [1, 3]", report, report)

    def test_extra_sample_in_final_tree_fails(self):
        """The other half of the mismatch check, OR'd into the same branch."""
        sampled = "EPI_A\tnode_1:A1G\n"
        final = "EPI_A\tnode_9:A1G\nEPI_B\tnode_9:A1G\n"
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

    def test_rebased_origin_is_reconciled_end_to_end(self):
        """The two trees record the same sequence against different origins.

        The final tree is rebased onto a root that differs from the reference
        at position 1, so its path omits the A1G the sampled tree records.
        Same sequence, and the gate must see that through --final-origin.
        """
        sampled = "EPI_A\tnode_1:A1G node_2:G3T\n"
        final = "EPI_A\tnode_9:G3T\n"
        result = self.run_check(sampled, final,
                                final_origin=">root\nGCGTACGTAC\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.report()["identical sequences"], 1)

    def test_rebased_origin_does_not_mask_a_real_difference(self):
        """The negative control: --final-origin must not paper over a change."""
        sampled = "EPI_A\tnode_1:A1G node_2:G3T\n"
        final = "EPI_A\tnode_9:G3C\n"
        result = self.run_check(sampled, final,
                                final_origin=">root\nGCGTACGTAC\n")
        self.assertEqual(result.returncode, 1)

    def test_origin_of_wrong_length_fails(self):
        sampled = "EPI_A\tnode_1:A1G\n"
        result = self.run_check(sampled, sampled, final_origin=">root\nACGT\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn("positions", result.stderr)


class ComposeWithRebasedOriginTestCase(unittest.TestCase):
    """A tree rebased onto its own root is reported against the reference."""

    def setUp(self):
        self.ref = {i: b for i, b in enumerate("ACGTACGTAC", start=1)}
        # Root differs from the reference at position 1 (A -> G).
        self.origin_diff = {1: "G"}

    def test_origin_difference_shows_up_with_an_empty_path(self):
        """A sample with no mutations still inherits the root's base."""
        self.assertEqual(compose("", self.ref, self.origin_diff), {1: "G"})

    def test_path_overrides_the_origin(self):
        """A back-mutation to the reference base cancels the origin difference."""
        self.assertEqual(compose("node_1:G1A", self.ref, self.origin_diff), {})

    def test_origin_and_path_combine(self):
        got = compose("node_1:G3T", self.ref, self.origin_diff)
        self.assertEqual(got, {1: "G", 3: "T"})

    def test_matches_the_unrebased_tree(self):
        """The point of the flag: same sequence, two ways of recording it."""
        on_reference = compose("node_1:A1G,G3T", self.ref)
        on_root = compose("node_1:G3T", self.ref, self.origin_diff)
        self.assertEqual(on_reference, on_root)


if __name__ == "__main__":
    unittest.main()
