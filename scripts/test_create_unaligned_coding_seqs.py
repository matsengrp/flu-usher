"""
Tests for create_unaligned_coding_seqs.py
"""

import unittest
import tempfile
import os
import lzma
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Import functions to test
import sys
sys.path.insert(0, os.path.dirname(__file__))
from create_unaligned_coding_seqs import (
    parse_insertions_from_tsv,
    remove_gaps,
    insert_nucleotides,
    validate_against_raw_sequences
)
from utils import sanitize_id, get_coding_region_coords


class TestSanitizeId(unittest.TestCase):
    """Tests for sanitize_id function"""

    def test_sanitize_basic(self):
        """Test basic sanitization of problematic characters"""
        self.assertEqual(sanitize_id("EPI_ISL_123456"), "EPI_ISL_123456")

    def test_sanitize_brackets(self):
        """Test removal of brackets"""
        self.assertEqual(sanitize_id("seq[variant]"), "seqvariant")
        self.assertEqual(sanitize_id("seq(2023)"), "seq2023")

    def test_sanitize_punctuation(self):
        """Test removal of various punctuation"""
        self.assertEqual(sanitize_id("seq:123;456,789"), "seq123456789")
        self.assertEqual(sanitize_id("seq'name.ext"), "seqnameext")

    def test_sanitize_multiple_chars(self):
        """Test removal of multiple problematic characters"""
        self.assertEqual(sanitize_id("A/H1N1[2023]:variant.1"), "A/H1N12023variant1")


class TestRemoveGaps(unittest.TestCase):
    """Tests for remove_gaps function"""

    def test_no_gaps(self):
        """Test sequence with no gaps"""
        self.assertEqual(remove_gaps("ACGTACGT"), "ACGTACGT")

    def test_with_gaps(self):
        """Test sequence with gaps"""
        self.assertEqual(remove_gaps("ACG-TAC-GT"), "ACGTACGT")

    def test_multiple_gaps(self):
        """Test sequence with multiple consecutive gaps"""
        self.assertEqual(remove_gaps("ACG---TAC---GT"), "ACGTACGT")

    def test_terminal_gaps(self):
        """Test sequence with gaps at ends"""
        self.assertEqual(remove_gaps("---ACGTACGT---"), "ACGTACGT")

    def test_all_gaps(self):
        """Test sequence that is all gaps"""
        self.assertEqual(remove_gaps("------"), "")


class TestInsertNucleotides(unittest.TestCase):
    """Tests for insert_nucleotides function"""

    def test_single_insertion(self):
        """Test single insertion"""
        seq = "ACGTACGT"
        # Insert "GGG" at position 4 (between T and A, neither of which is G)
        insertions = [(4, "GGG")]
        result = insert_nucleotides(seq, insertions)
        self.assertEqual(result, "ACGTGGGACGT")

    def test_multiple_insertions_descending(self):
        """Test multiple insertions in descending order"""
        seq = "ACGTACGT"
        # Insert "TTT" at position 6 (between C and G, neither of which is T)
        # Insert "AAA" at position 3 (between G and T, neither of which is A)
        insertions = [(6, "TTT"), (3, "AAA")]  # Already in descending order
        result = insert_nucleotides(seq, insertions)
        self.assertEqual(result, "ACGAAATACTTTGT")

    def test_insertion_at_start(self):
        """Test insertion at position 0 (beginning)"""
        seq = "ACGTACGT"
        insertions = [(0, "TTT")]
        result = insert_nucleotides(seq, insertions)
        self.assertEqual(result, "TTTACGTACGT")

    def test_insertion_at_end(self):
        """Test insertion at end"""
        seq = "ACGTACGT"
        # Insert "AAA" at position 8 (at the end, after T which is not A)
        insertions = [(8, "AAA")]
        result = insert_nucleotides(seq, insertions)
        self.assertEqual(result, "ACGTACGTAAA")

    def test_no_insertions(self):
        """Test with no insertions"""
        seq = "ACGTACGT"
        insertions = []
        result = insert_nucleotides(seq, insertions)
        self.assertEqual(result, "ACGTACGT")

    def test_insertion_order_with_gaps(self):
        """Test that insertions work correctly when applied to aligned sequences with gaps

        This test validates the bug fix: insertions must be applied to the aligned sequence
        BEFORE removing gaps. If insertions are applied after removing gaps, the positions
        will be incorrect because gaps shift all downstream positions.
        """
        # Aligned sequence with gaps at positions 5-10
        aligned_seq = "ACTC------ATTGC"

        # Insertion at position 11 (1-based, meaning "insert after position 11")
        # Position 11 is the 'A' character right after the gap region
        insertions = [(11, "XXX")]

        # Correct workflow: Apply insertions to aligned sequence FIRST
        seq_with_insertions = insert_nucleotides(aligned_seq, insertions)
        # Expected: "ACTC------AXXXTTGC" (insertions go after position 11 in aligned coords)
        self.assertEqual(seq_with_insertions, "ACTC------AXXXTTGC")

        # Then remove gaps
        final_seq = remove_gaps(seq_with_insertions)
        # Expected: "ACTCAXXXTTGC" (12 chars: 4 + 1 + 3 + 4)
        self.assertEqual(final_seq, "ACTCAXXXTTGC")

        # Wrong workflow would be: remove gaps first (8 chars), then try to insert at
        # position 11, which is past the end or at the wrong location


class TestGetCodingRegionCoords(unittest.TestCase):
    """Tests for get_coding_region_coords function"""

    def test_single_cds(self):
        """Test GFF with single CDS feature"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("##gff-version 3\n")
            f.write("ref\tNCBI\tCDS\t1\t1683\t.\t+\t0\tID=cds1;Name=HA\n")
            f.flush()

            try:
                min_start, max_end = get_coding_region_coords(f.name)
                self.assertEqual(min_start, 1)
                self.assertEqual(max_end, 1683)
            finally:
                os.unlink(f.name)

    def test_multiple_cds(self):
        """Test GFF with multiple CDS features"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("##gff-version 3\n")
            f.write("ref\tNCBI\tCDS\t10\t500\t.\t+\t0\tID=cds1\n")
            f.write("ref\tNCBI\tCDS\t600\t1000\t.\t+\t0\tID=cds2\n")
            f.flush()

            try:
                min_start, max_end = get_coding_region_coords(f.name)
                self.assertEqual(min_start, 10)
                self.assertEqual(max_end, 1000)
            finally:
                os.unlink(f.name)

    def test_gene_feature(self):
        """Test GFF with gene feature"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("##gff-version 3\n")
            f.write("ref\tNCBI\tgene\t50\t1500\t.\t+\t.\tID=gene1\n")
            f.flush()

            try:
                min_start, max_end = get_coding_region_coords(f.name)
                self.assertEqual(min_start, 50)
                self.assertEqual(max_end, 1500)
            finally:
                os.unlink(f.name)


class TestParseInsertionsFromTsv(unittest.TestCase):
    """Tests for parse_insertions_from_tsv function"""

    def test_basic_insertion(self):
        """Test basic insertion parsing"""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.tsv.xz', delete=False) as f:
            content = "seqName\tinsertions\n"
            content += "EPI_ISL_123\t100:ACG\n"
            content += "EPI_ISL_456\t\n"

            with lzma.open(f.name, 'wt') as xz_file:
                xz_file.write(content)

            try:
                result = parse_insertions_from_tsv(f.name, 1, 1683)
                self.assertIn("EPI_ISL_123", result)
                self.assertEqual(result["EPI_ISL_123"], [(100, "ACG")])
                self.assertNotIn("EPI_ISL_456", result)
            finally:
                os.unlink(f.name)

    def test_multiple_insertions(self):
        """Test parsing multiple insertions for one sequence"""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.tsv.xz', delete=False) as f:
            content = "seqName\tinsertions\n"
            content += "EPI_ISL_123\t100:ACG,500:TTT,800:GGG\n"

            with lzma.open(f.name, 'wt') as xz_file:
                xz_file.write(content)

            try:
                result = parse_insertions_from_tsv(f.name, 1, 1683)
                self.assertIn("EPI_ISL_123", result)
                # Should be sorted in descending order
                self.assertEqual(result["EPI_ISL_123"], [(800, "GGG"), (500, "TTT"), (100, "ACG")])
            finally:
                os.unlink(f.name)

    def test_filter_insertions_outside_cds(self):
        """Test filtering of insertions outside coding region"""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.tsv.xz', delete=False) as f:
            content = "seqName\tinsertions\n"
            # Position 0: before CDS (min_start=1)
            # Position 100: inside CDS
            # Position 2000: after CDS (max_end=1683)
            content += "EPI_ISL_123\t0:ACG,100:TTT,2000:GGG\n"

            with lzma.open(f.name, 'wt') as xz_file:
                xz_file.write(content)

            try:
                result = parse_insertions_from_tsv(f.name, 1, 1683)
                self.assertIn("EPI_ISL_123", result)
                # Only position 100 should be kept
                self.assertEqual(result["EPI_ISL_123"], [(100, "TTT")])
            finally:
                os.unlink(f.name)

    def test_adjust_positions_for_cds_slice(self):
        """Test that positions are adjusted for CDS slice"""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.tsv.xz', delete=False) as f:
            content = "seqName\tinsertions\n"
            # If min_start=50, position 100 should become (100-50+1)=51
            content += "EPI_ISL_123\t100:ACG\n"

            with lzma.open(f.name, 'wt') as xz_file:
                xz_file.write(content)

            try:
                result = parse_insertions_from_tsv(f.name, 50, 1683)
                self.assertIn("EPI_ISL_123", result)
                self.assertEqual(result["EPI_ISL_123"], [(51, "ACG")])
            finally:
                os.unlink(f.name)

    def test_sanitize_seq_ids(self):
        """Test that sequence IDs are sanitized"""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.tsv.xz', delete=False) as f:
            content = "seqName\tinsertions\n"
            content += "EPI[ISL]_123:456\t100:ACG\n"

            with lzma.open(f.name, 'wt') as xz_file:
                xz_file.write(content)

            try:
                result = parse_insertions_from_tsv(f.name, 1, 1683)
                # Sanitized ID should remove brackets and colon
                self.assertIn("EPIISL_123456", result)
                self.assertEqual(result["EPIISL_123456"], [(100, "ACG")])
            finally:
                os.unlink(f.name)


class TestValidateAgainstRawSequences(unittest.TestCase):
    """Tests for validate_against_raw_sequences function"""

    def test_valid_substring(self):
        """Test that unaligned sequence is substring of raw sequence"""
        # Create test sequences
        unaligned_records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1", description="seq1")
        ]

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.fasta.xz', delete=False) as f:
            with lzma.open(f.name, 'wt') as handle:
                # Raw sequence contains the unaligned sequence
                SeqIO.write([SeqRecord(Seq("TTTTACGTACGTGGGG"), id="seq1", description="seq1")], handle, 'fasta')

            try:
                num_validated, num_failed = validate_against_raw_sequences(unaligned_records, f.name)
                self.assertEqual(num_validated, 1)
                self.assertEqual(num_failed, 0)
            finally:
                os.unlink(f.name)

    def test_invalid_not_substring(self):
        """Test that validation fails when unaligned sequence is not substring"""
        # Create test sequences
        unaligned_records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1", description="seq1")
        ]

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.fasta.xz', delete=False) as f:
            with lzma.open(f.name, 'wt') as handle:
                # Raw sequence does NOT contain the unaligned sequence
                SeqIO.write([SeqRecord(Seq("TTTTGGGGCCCC"), id="seq1", description="seq1")], handle, 'fasta')

            try:
                num_validated, num_failed = validate_against_raw_sequences(unaligned_records, f.name)
                self.assertEqual(num_validated, 0)
                self.assertEqual(num_failed, 1)
            finally:
                os.unlink(f.name)

    def test_case_insensitive(self):
        """Test that validation is case-insensitive"""
        # Create test sequences with mixed case
        unaligned_records = [
            SeqRecord(Seq("acgtACGT"), id="seq1", description="seq1")
        ]

        with tempfile.NamedTemporaryFile(mode='wb', suffix='.fasta.xz', delete=False) as f:
            with lzma.open(f.name, 'wt') as handle:
                # Raw sequence in different case
                SeqIO.write([SeqRecord(Seq("TTTTacgtACGTGGGG"), id="seq1", description="seq1")], handle, 'fasta')

            try:
                num_validated, num_failed = validate_against_raw_sequences(unaligned_records, f.name)
                self.assertEqual(num_validated, 1)
                self.assertEqual(num_failed, 0)
            finally:
                os.unlink(f.name)


if __name__ == "__main__":
    unittest.main()
