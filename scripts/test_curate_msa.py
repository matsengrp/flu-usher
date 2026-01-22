"""
Tests for curate_msa.py
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
from curate_and_extract_coding_seqs import (
    create_matching_gff_and_gtf,
    slice_record,
    get_ambiguous_chars,
    analyze_record,
    filter_sequences
)
from utils import sanitize_id, extract_all_genes_and_cds


class TestSliceRecord(unittest.TestCase):
    """Tests for slice_record function"""

    def test_basic_slicing(self):
        """Test basic sequence slicing"""
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1", description="test")
        sliced = slice_record(record, 3, 8)
        # Positions 3-8 (1-based) = indices 2-8 (0-based exclusive) = "GTACGT"
        self.assertEqual(str(sliced.seq), "GTACGT")
        self.assertEqual(sliced.id, "seq1")

    def test_slicing_with_sanitization(self):
        """Test that sequence IDs are sanitized"""
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq[1]:test", description="test")
        sliced = slice_record(record, 1, 12)
        # ID should be sanitized
        self.assertEqual(sliced.id, "seq1test")

    def test_full_sequence(self):
        """Test slicing entire sequence"""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1", description="test")
        sliced = slice_record(record, 1, 8)
        self.assertEqual(str(sliced.seq), "ACGTACGT")

    def test_single_position(self):
        """Test slicing single position"""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1", description="test")
        sliced = slice_record(record, 3, 3)
        self.assertEqual(str(sliced.seq), "G")


class TestGetAmbiguousChars(unittest.TestCase):
    """Tests for get_ambiguous_chars function"""

    def test_no_ambiguous(self):
        """Test with only standard nucleotides"""
        records = [
            SeqRecord(Seq("ACGT-ACGT"), id="seq1"),
            SeqRecord(Seq("GCTA-GCTA"), id="seq2")
        ]
        ambig = get_ambiguous_chars(records)
        self.assertEqual(ambig, set())

    def test_with_n(self):
        """Test with N characters"""
        records = [
            SeqRecord(Seq("ACGTNACGT"), id="seq1"),
            SeqRecord(Seq("GCTA-GCTA"), id="seq2")
        ]
        ambig = get_ambiguous_chars(records)
        self.assertEqual(ambig, {'N'})

    def test_multiple_ambiguous(self):
        """Test with multiple ambiguous characters"""
        records = [
            SeqRecord(Seq("ACGTNRYWSKM"), id="seq1"),
        ]
        ambig = get_ambiguous_chars(records)
        # N, R, Y, W, S, K, M are all ambiguous
        self.assertEqual(ambig, {'N', 'R', 'Y', 'W', 'S', 'K', 'M'})


class TestAnalyzeRecord(unittest.TestCase):
    """Tests for analyze_record function"""

    def test_no_gaps_no_ambiguous(self):
        """Test sequence with no gaps or ambiguous characters"""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        ambig_chars = set()
        result = analyze_record(record, ambig_chars)

        self.assertEqual(result['id'], "seq1")
        self.assertEqual(result['length'], 8)
        self.assertEqual(result['frac_gaps'], 0.0)
        self.assertEqual(result['frac_ambiguous'], 0.0)

    def test_with_gaps(self):
        """Test sequence with gaps"""
        record = SeqRecord(Seq("ACGT--ACGT"), id="seq1")
        ambig_chars = set()
        result = analyze_record(record, ambig_chars)

        self.assertEqual(result['length'], 10)
        self.assertEqual(result['frac_gaps'], 0.2)  # 2 gaps out of 10
        self.assertEqual(result['frac_ambiguous'], 0.0)

    def test_with_ambiguous(self):
        """Test sequence with ambiguous characters"""
        record = SeqRecord(Seq("ACGTNNNACGT"), id="seq1")
        ambig_chars = {'N'}
        result = analyze_record(record, ambig_chars)

        self.assertEqual(result['length'], 11)
        self.assertEqual(result['frac_gaps'], 0.0)
        self.assertAlmostEqual(result['frac_ambiguous'], 3/11, places=5)

    def test_with_gaps_and_ambiguous(self):
        """Test sequence with both gaps and ambiguous characters"""
        record = SeqRecord(Seq("ACGT--NNNACGT"), id="seq1")
        ambig_chars = {'N'}
        result = analyze_record(record, ambig_chars)

        self.assertEqual(result['length'], 13)
        self.assertAlmostEqual(result['frac_gaps'], 2/13, places=5)
        self.assertAlmostEqual(result['frac_ambiguous'], 3/13, places=5)


class TestFilterSequences(unittest.TestCase):
    """Tests for filter_sequences function"""

    def setUp(self):
        """Set up test logger"""
        import logging
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.WARNING)  # Reduce noise in tests

    def test_filter_by_gaps(self):
        """Test filtering sequences by gap content"""
        records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1"),        # 0% gaps - KEEP
            SeqRecord(Seq("ACGT--ACGT"), id="seq2"),      # 20% gaps - FILTER
            SeqRecord(Seq("ACGT-ACGT"), id="seq3"),       # 11% gaps - FILTER
        ]
        ambig_chars = set()

        filtered = filter_sequences(records, ambig_chars, max_frac_gaps=0.05,
                                   max_frac_ambig=0.01, logger=self.logger)

        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered[0].id, "seq1")

    def test_filter_by_ambiguous(self):
        """Test filtering sequences by ambiguous nucleotide content"""
        records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1"),        # 0% ambig - KEEP
            SeqRecord(Seq("ACGTNNNACGT"), id="seq2"),     # 27% ambig - FILTER
            SeqRecord(Seq("ACGTNACGT"), id="seq3"),       # 11% ambig - FILTER
        ]
        ambig_chars = {'N'}

        filtered = filter_sequences(records, ambig_chars, max_frac_gaps=0.05,
                                   max_frac_ambig=0.01, logger=self.logger)

        self.assertEqual(len(filtered), 1)
        self.assertEqual(filtered[0].id, "seq1")

    def test_filter_terminal_gaps(self):
        """Test filtering sequences with terminal gaps"""
        records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1"),        # No terminal gaps - KEEP
            SeqRecord(Seq("---ACGTACGT"), id="seq2"),     # Terminal gaps at start - FILTER
            SeqRecord(Seq("ACGTACGT---"), id="seq3"),     # Terminal gaps at end - FILTER
            SeqRecord(Seq("ACGT---ACGT"), id="seq4"),     # Internal gaps only - KEEP
        ]
        ambig_chars = set()

        filtered = filter_sequences(records, ambig_chars, max_frac_gaps=1.0,
                                   max_frac_ambig=1.0, logger=self.logger)

        self.assertEqual(len(filtered), 2)
        self.assertEqual(filtered[0].id, "seq1")
        self.assertEqual(filtered[1].id, "seq4")

    def test_replace_ambiguous_with_n(self):
        """Test that ambiguous characters are replaced with N"""
        records = [
            SeqRecord(Seq("ACGTRYWACGT"), id="seq1"),
        ]
        ambig_chars = {'R', 'Y', 'W'}

        filtered = filter_sequences(records, ambig_chars, max_frac_gaps=1.0,
                                   max_frac_ambig=1.0, logger=self.logger)

        self.assertEqual(len(filtered), 1)
        # R, Y, W should all be replaced with N
        self.assertEqual(str(filtered[0].seq), "ACGTNNNACGT")

    def test_filter_duplicates(self):
        """Test filtering duplicate sequences"""
        records = [
            SeqRecord(Seq("ACGTACGT"), id="seq1"),        # First occurrence - KEEP
            SeqRecord(Seq("GCTAGCTA"), id="seq2"),        # Different - KEEP
            SeqRecord(Seq("ACGTACGT"), id="seq3"),        # Duplicate - FILTER
            SeqRecord(Seq("GCTAGCTA"), id="seq4"),        # Duplicate - FILTER
        ]
        ambig_chars = set()

        filtered = filter_sequences(records, ambig_chars, max_frac_gaps=1.0,
                                   max_frac_ambig=1.0, logger=self.logger,
                                   filter_duplicates=True)

        self.assertEqual(len(filtered), 2)
        self.assertEqual(filtered[0].id, "seq1")
        self.assertEqual(filtered[1].id, "seq2")

    def test_replace_gaps_with_ref(self):
        """Test replacing gaps with reference nucleotides"""
        records = [
            SeqRecord(Seq("ACGTACGT"), id="ref"),         # Reference (no gaps)
            SeqRecord(Seq("ACGT--GT"), id="seq1"),        # Has gaps
            SeqRecord(Seq("AC--ACGT"), id="seq2"),        # Has gaps
        ]
        ambig_chars = set()

        filtered = filter_sequences(records, ambig_chars, max_frac_gaps=1.0,
                                   max_frac_ambig=1.0, logger=self.logger,
                                   replace_gaps_with_ref=True)

        self.assertEqual(len(filtered), 3)
        # Reference should be unchanged
        self.assertEqual(str(filtered[0].seq), "ACGTACGT")
        # Gaps should be replaced with reference nucleotides
        self.assertEqual(str(filtered[1].seq), "ACGTACGT")
        self.assertEqual(str(filtered[2].seq), "ACGTACGT")


class TestCreateMatchingGffAndGtf(unittest.TestCase):
    """Tests for create_matching_gff_and_gtf function"""

    def test_single_feature(self):
        """Test creating GFF/GTF with single feature"""
        features = [{
            'id': 'cds1',
            'name': 'HA',
            'type': 'CDS',
            'start': 1,
            'end': 1683,
            'strand': '+',
            'phase': '0',
            'attributes': 'ID=cds1;Name=HA'
        }]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as gff_file:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as gtf_file:
                try:
                    create_matching_gff_and_gtf(gff_file.name, gtf_file.name, features, min_start=1)

                    # Check GFF file
                    with open(gff_file.name, 'r') as f:
                        gff_content = f.read()
                        self.assertIn('##gff-version 3', gff_content)
                        self.assertIn('##sequence-region Reference 1 1683', gff_content)
                        self.assertIn('CDS\t1\t1683', gff_content)

                    # Check GTF file
                    with open(gtf_file.name, 'r') as f:
                        gtf_content = f.read()
                        self.assertIn('CDS\t1\t1683', gtf_content)
                        self.assertIn('gene_id "cds1"', gtf_content)
                        self.assertIn('gene_name "HA"', gtf_content)

                finally:
                    os.unlink(gff_file.name)
                    os.unlink(gtf_file.name)

    def test_coordinate_adjustment(self):
        """Test that coordinates are adjusted correctly"""
        features = [{
            'id': 'cds1',
            'name': 'HA',
            'type': 'CDS',
            'start': 50,
            'end': 1733,
            'strand': '+',
            'phase': '0',
            'attributes': 'ID=cds1;Name=HA'
        }]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as gff_file:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as gtf_file:
                try:
                    # min_start=50 means position 50 becomes position 1
                    create_matching_gff_and_gtf(gff_file.name, gtf_file.name, features, min_start=50)

                    # Check GFF file - coordinates should be adjusted
                    with open(gff_file.name, 'r') as f:
                        gff_content = f.read()
                        self.assertIn('##sequence-region Reference 1 1684', gff_content)
                        self.assertIn('CDS\t1\t1684', gff_content)

                finally:
                    os.unlink(gff_file.name)
                    os.unlink(gtf_file.name)

    def test_multiple_features(self):
        """Test creating GFF/GTF with multiple features"""
        features = [
            {
                'id': 'gene1',
                'name': 'HA',
                'type': 'gene',
                'start': 1,
                'end': 1683,
                'strand': '+',
                'phase': '.',
                'attributes': 'ID=gene1;Name=HA'
            },
            {
                'id': 'cds1',
                'name': 'HA_CDS',
                'type': 'CDS',
                'start': 1,
                'end': 1683,
                'strand': '+',
                'phase': '0',
                'attributes': 'ID=cds1;Name=HA_CDS;Parent=gene1'
            }
        ]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as gff_file:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as gtf_file:
                try:
                    create_matching_gff_and_gtf(gff_file.name, gtf_file.name, features, min_start=1)

                    # Check that both features are in GFF
                    with open(gff_file.name, 'r') as f:
                        gff_content = f.read()
                        self.assertIn('gene\t1\t1683', gff_content)
                        self.assertIn('CDS\t1\t1683', gff_content)

                    # Check that both features are in GTF
                    with open(gtf_file.name, 'r') as f:
                        gtf_content = f.read()
                        self.assertIn('gene\t1\t1683', gtf_content)
                        self.assertIn('CDS\t1\t1683', gtf_content)

                finally:
                    os.unlink(gff_file.name)
                    os.unlink(gtf_file.name)


if __name__ == "__main__":
    unittest.main()
