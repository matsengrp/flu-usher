"""
Tests for create_unaligned_coding_seqs.py
"""

import unittest
import tempfile
import os
import lzma
import logging
from io import StringIO
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
        # Using GGG (doesn't match adjacent 'A' to avoid ambiguity)
        insertions = [(11, "GGG")]

        # Correct workflow: Apply insertions to aligned sequence FIRST
        seq_with_insertions = insert_nucleotides(aligned_seq, insertions)
        # Expected: "ACTC------AGGGTTGC" (insertions go after position 11 in aligned coords)
        self.assertEqual(seq_with_insertions, "ACTC------AGGGTTGC")

        # Then remove gaps
        final_seq = remove_gaps(seq_with_insertions)
        # Expected: "ACTCAGGGTTGC" (12 chars: 4 + 1 + 3 + 4)
        self.assertEqual(final_seq, "ACTCAGGGTTGC")

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
                result = parse_insertions_from_tsv(f.name)
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
                result = parse_insertions_from_tsv(f.name)
                self.assertIn("EPI_ISL_123", result)
                # Should be sorted in descending order
                self.assertEqual(result["EPI_ISL_123"], [(800, "GGG"), (500, "TTT"), (100, "ACG")])
            finally:
                os.unlink(f.name)

    def test_all_insertions_kept(self):
        """Test that all insertions are kept (no filtering)"""
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.tsv.xz', delete=False) as f:
            content = "seqName\tinsertions\n"
            # All positions should be kept (no filtering)
            content += "EPI_ISL_123\t10:ACG,100:TTT,2000:GGG\n"

            with lzma.open(f.name, 'wt') as xz_file:
                xz_file.write(content)

            try:
                result = parse_insertions_from_tsv(f.name)
                self.assertIn("EPI_ISL_123", result)
                # All insertions should be kept in descending order
                self.assertEqual(result["EPI_ISL_123"], [(2000, "GGG"), (100, "TTT"), (10, "ACG")])
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
                result = parse_insertions_from_tsv(f.name)
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


# ============================================================================
# TESTS FOR PER-GENE CDS EXTRACTION
# ============================================================================
#
# These tests cover functions for extracting individual CDS per gene with
# support for spliced genes and biological validation:
# - group_cds_by_gene(): Group CDS features by gene/protein (in utils.py)
# - filter_insertions_for_cds(): Filter insertions for specific CDS regions
# - extract_cds_from_aligned(): Extract CDS from aligned sequences
# - extract_gene_cds(): Complete gene CDS extraction with inline validation
# - validate_cds(): CDS biological validation (frame, start codon, stop codon)
# ============================================================================

# Test data constants
VALID_CDS_PERGENE = "ATGCGATCGTAA"  # 12 bp, starts with ATG, ends with TAA
ALIGNED_CDS_WITH_GAPS_PERGENE = "ATG---CGATCG---TAA"  # 12 bp ungapped
LONG_CDS_PERGENE = "ATGCGATCGAAACCGTTCGGTTGA"  # 24 bp
LONG_ALIGNED_CDS_PERGENE = "ATG---CGATCG---AAACCG---TTCGGTTGA"  # 24 bp ungapped
INVALID_FRAME_CDS_PERGENE = "ATGCGATCGTAAG"  # 13 bp, not divisible by 3
NEP_FRAGMENT1_PERGENE = "ATGGATTCCAACACTGTGTCAAGCTTTCAA"  # 30 bp
NEP_FRAGMENT2_PERGENE = "AAACCGTTCTAA"  # 12 bp, total 42 bp when concatenated


class TestGroupCdsByGene(unittest.TestCase):
    """Tests for group_cds_by_gene() function"""

    def test_single_gene_single_cds(self):
        """Single CDS feature for HA gene"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'CDS',
                'name': 'HA',
                'start': 1,
                'end': 1683,
                'attributes': 'gene=HA;protein_id=P12345'
            }
        ]
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 1)
        self.assertIn('P12345', result)
        self.assertEqual(result['P12345']['gene_name'], 'HA')
        self.assertEqual(len(result['P12345']['cds_list']), 1)

    def test_spliced_gene_nep(self):
        """Two CDS features with same protein_id (NEP case)"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'CDS',
                'name': 'NEP',
                'start': 1,
                'end': 30,
                'attributes': 'gene=NEP;protein_id=P67890'
            },
            {
                'type': 'CDS',
                'name': 'NEP',
                'start': 503,
                'end': 838,
                'attributes': 'gene=NEP;protein_id=P67890'
            }
        ]
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 1)
        self.assertIn('P67890', result)
        self.assertEqual(result['P67890']['gene_name'], 'NEP')
        self.assertEqual(len(result['P67890']['cds_list']), 2)

    def test_multiple_genes_no_splicing(self):
        """NS1 CDS and NEP CDS fragments (different protein_ids)"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'CDS',
                'name': 'NS1',
                'start': 1,
                'end': 693,
                'attributes': 'gene=NS1;protein_id=P11111'
            },
            {
                'type': 'CDS',
                'name': 'NEP',
                'start': 1,
                'end': 30,
                'attributes': 'gene=NEP;protein_id=P22222'
            },
            {
                'type': 'CDS',
                'name': 'NEP',
                'start': 503,
                'end': 838,
                'attributes': 'gene=NEP;protein_id=P22222'
            }
        ]
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 2)
        self.assertIn('P11111', result)
        self.assertIn('P22222', result)
        self.assertEqual(len(result['P11111']['cds_list']), 1)
        self.assertEqual(len(result['P22222']['cds_list']), 2)

    def test_cds_sorting_by_start_position(self):
        """CDS features provided out of order should be sorted"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'CDS',
                'name': 'NEP',
                'start': 503,
                'end': 838,
                'attributes': 'gene=NEP;protein_id=P67890'
            },
            {
                'type': 'CDS',
                'name': 'NEP',
                'start': 1,
                'end': 30,
                'attributes': 'gene=NEP;protein_id=P67890'
            }
        ]
        result = group_cds_by_gene(features)

        cds_list = result['P67890']['cds_list']
        self.assertEqual(cds_list[0]['start'], 1)
        self.assertEqual(cds_list[1]['start'], 503)

    def test_gene_without_protein_id(self):
        """CDS feature with gene attribute but no protein_id"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'CDS',
                'name': 'HA',
                'start': 1,
                'end': 1683,
                'attributes': 'gene=HA'
            }
        ]
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 1)
        # Should use gene name as key when protein_id is missing
        self.assertIn('HA', result)

    def test_mixed_protein_ids(self):
        """Some CDS with protein_id, some without"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'CDS',
                'name': 'HA',
                'start': 1,
                'end': 1683,
                'attributes': 'gene=HA;protein_id=P12345'
            },
            {
                'type': 'CDS',
                'name': 'NA',
                'start': 1,
                'end': 1410,
                'attributes': 'gene=NA'
            }
        ]
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 2)
        self.assertIn('P12345', result)
        self.assertIn('NA', result)

    def test_non_cds_features_ignored(self):
        """Mix of CDS, gene, mRNA features - only CDS should be grouped"""
        from utils import group_cds_by_gene

        features = [
            {
                'type': 'gene',
                'name': 'HA',
                'start': 1,
                'end': 1683,
                'attributes': 'gene=HA'
            },
            {
                'type': 'CDS',
                'name': 'HA',
                'start': 1,
                'end': 1683,
                'attributes': 'gene=HA;protein_id=P12345'
            },
            {
                'type': 'mRNA',
                'name': 'HA',
                'start': 1,
                'end': 1683,
                'attributes': 'gene=HA'
            }
        ]
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 1)
        self.assertIn('P12345', result)

    def test_empty_feature_list(self):
        """Empty list should return empty dict"""
        from utils import group_cds_by_gene

        features = []
        result = group_cds_by_gene(features)

        self.assertEqual(len(result), 0)
        self.assertIsInstance(result, dict)


class TestFilterInsertionsForCds(unittest.TestCase):
    """Tests for filter_insertions_for_cds() function"""

    def test_no_insertions(self):
        """Empty insertion list should return empty list"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = []
        offset = 0  # No offset for this test
        result = filter_insertions_for_cds(insertions, 50, 150, offset)

        self.assertEqual(len(result), 0)
        self.assertIsInstance(result, list)

    def test_single_codon_insertion_within_cds(self):
        """Insertion at position 50 within CDS 1-100"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(50, "AAA")]
        offset = 0  # Original coords = curated coords
        result = filter_insertions_for_cds(insertions, 1, 100, offset)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (50, "AAA"))

    def test_insertion_before_cds(self):
        """Insertion before CDS start should be filtered out"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(10, "GGG")]
        offset = 0
        result = filter_insertions_for_cds(insertions, 50, 150, offset)

        self.assertEqual(len(result), 0)

    def test_insertion_after_cds(self):
        """Insertion after CDS end should be filtered out"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(200, "TTT")]
        offset = 0
        result = filter_insertions_for_cds(insertions, 50, 150, offset)

        self.assertEqual(len(result), 0)

    def test_insertion_at_cds_start_boundary(self):
        """Insertion at CDS start boundary should be included"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(50, "CCC")]
        offset = 49  # offset = min_start - 1 = 50 - 1 = 49
        result = filter_insertions_for_cds(insertions, 50, 150, offset)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (1, "CCC"))  # 50 - 49 = 1 (curated position)

    def test_insertion_at_cds_end_boundary(self):
        """Insertion at position 149 (before end 150) should be included"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(149, "GGG")]
        offset = 49  # offset = 50 - 1 = 49
        result = filter_insertions_for_cds(insertions, 50, 150, offset)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (100, "GGG"))  # 149 - 49 = 100 (curated position)

    def test_multiple_codon_insertions_mixed(self):
        """Multiple insertions, some within CDS, some outside"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(200, "AAA"), (100, "GGG"), (60, "TTT"), (10, "CCC")]
        offset = 49  # offset = 50 - 1 = 49
        result = filter_insertions_for_cds(insertions, 50, 150, offset)

        self.assertEqual(len(result), 2)
        # Should keep positions 60 and 100, transformed to curated space and in descending order
        self.assertEqual(result[0], (51, "GGG"))  # 100 - 49 = 51
        self.assertEqual(result[1], (11, "TTT"))  # 60 - 49 = 11

    def test_maintains_descending_order(self):
        """Output should maintain descending order"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(90, "AAA"), (60, "GGG"), (30, "TTT")]
        offset = 0  # No offset
        result = filter_insertions_for_cds(insertions, 1, 100, offset)

        self.assertEqual(len(result), 3)
        self.assertEqual(result[0][0], 90)  # Highest position first
        self.assertEqual(result[1][0], 60)
        self.assertEqual(result[2][0], 30)

    def test_position_adjustment_calculation(self):
        """Position adjustment should transform from original to curated coords"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        insertions = [(150, "CCC")]
        offset = 99  # offset = 100 - 1 = 99
        result = filter_insertions_for_cds(insertions, 100, 200, offset)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (51, "CCC"))  # 150 - 99 = 51

    def test_realistic_ha_coordinates_with_offset(self):
        """Test with realistic HA coordinates (5' UTR causes offset)"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        # Realistic scenario: HA CDS starts at position 58, ends at 1708 in original ref
        # After curation, MSA is positions 1-1651
        # offset = 57 (58 - 1)
        original_cds_start = 58
        original_cds_end = 1708
        offset = 57

        # Insertion at position 100 in original reference
        insertions = [(100, "AAA")]
        result = filter_insertions_for_cds(insertions, original_cds_start, original_cds_end, offset)

        # Should transform to position 43 in curated space (100 - 57)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (43, "AAA"))

    def test_insertion_at_boundaries_with_large_offset(self):
        """Test insertion filtering at boundaries with large offset"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        # Original CDS: 100-500, offset = 99
        # Curated CDS: 1-401
        original_start = 100
        original_end = 500
        offset = 99

        # Insertions at various positions
        insertions = [
            (99, "AAA"),   # Before CDS start (should be filtered)
            (100, "GGG"),  # At CDS start (should be included → position 1)
            (499, "TTT"),  # Just before CDS end (should be included → position 400)
            (500, "CCC"),  # At CDS end (should be filtered, exclusive boundary)
        ]

        result = filter_insertions_for_cds(insertions, original_start, original_end, offset)

        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], (400, "TTT"))
        self.assertEqual(result[1], (1, "GGG"))

    def test_zero_offset_no_transformation(self):
        """Test that zero offset doesn't break coordinate transformation"""
        from create_unaligned_coding_seqs import filter_insertions_for_cds

        offset = 0

        # When offset is 0, original coords = curated coords
        insertions = [(50, "AAA")]
        result = filter_insertions_for_cds(insertions, 1, 100, offset)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], (50, "AAA"))  # No transformation


class TestExtractCdsFromAligned(unittest.TestCase):
    """Tests for extract_cds_from_aligned() function"""

    def test_extract_cds_from_middle(self):
        """Extract CDS from middle of aligned sequence in curated MSA"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        # Curated MSA already has non-coding regions removed
        # This represents the curated MSA (coding region only)
        aligned_seq = "ATG---CGATCG---TAANNNNNN"
        # Original reference had CDS at positions 7-24
        # With offset = 6, curated coords are 1-18
        offset = 6
        result = extract_cds_from_aligned(aligned_seq, 7, 24, offset)

        expected = "ATG---CGATCG---TAA"
        self.assertEqual(result, expected)

    def test_extract_cds_from_beginning(self):
        """Extract CDS from beginning of sequence"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        aligned_seq = "ATG---CGATCG---TAANNNNNN"
        offset = 0  # CDS starts at position 1
        result = extract_cds_from_aligned(aligned_seq, 1, 18, offset)

        expected = "ATG---CGATCG---TAA"
        self.assertEqual(result, expected)

    def test_extract_cds_to_end(self):
        """Extract CDS to end of curated MSA sequence"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        # Curated MSA (coding region only)
        aligned_seq = "ATG---CGATCG---TAA"
        # Original coords: 7-24, offset = 6
        # Curated coords: 1-18
        offset = 6
        result = extract_cds_from_aligned(aligned_seq, 7, 24, offset)

        expected = "ATG---CGATCG---TAA"
        self.assertEqual(result, expected)

    def test_extract_entire_cds(self):
        """Extract entire sequence as CDS"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        aligned_seq = "ATG---CGATCG---TAA"
        offset = 0  # No offset
        result = extract_cds_from_aligned(aligned_seq, 1, 19, offset)

        self.assertEqual(result, aligned_seq)

    def test_extract_longer_cds_with_gaps(self):
        """Extract longer CDS with multiple gap regions"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        aligned_seq = "ATGCGA---TCG---AAACCG---TTCTAA"
        offset = 0  # No offset
        result = extract_cds_from_aligned(aligned_seq, 1, 33, offset)

        self.assertEqual(result, aligned_seq)
        # After gap removal would be 24 bp

    def test_preserve_gaps_in_extraction(self):
        """Gaps should be preserved in extracted sequence"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        aligned_seq = "ATG------CGATCGTAA"
        offset = 0  # No offset
        result = extract_cds_from_aligned(aligned_seq, 1, 18, offset)

        self.assertEqual(result, aligned_seq)
        self.assertIn("------", result)

    def test_zero_based_conversion(self):
        """1-based coordinates should be converted correctly to 0-based"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        aligned_seq = "ATGCGATCGTAA"
        offset = 0  # No offset
        result = extract_cds_from_aligned(aligned_seq, 1, 12, offset)

        # Should extract entire sequence [0:12]
        self.assertEqual(result, aligned_seq)

    def test_extract_with_zero_offset(self):
        """Test extraction with zero offset (no coordinate transformation)"""
        from create_unaligned_coding_seqs import extract_cds_from_aligned

        aligned_seq = "ATG---CGATCG---TAA"
        offset = 0
        result = extract_cds_from_aligned(aligned_seq, 1, 19, offset)

        # With zero offset, coords don't change
        self.assertEqual(result, aligned_seq)


class TestExtractGeneCds(unittest.TestCase):
    """Tests for extract_gene_cds() function with inline validation"""

    def setUp(self):
        """Set up logger for tests"""
        self.logger = logging.getLogger('test_logger_new')
        self.logger.setLevel(logging.DEBUG)
        # Add handler to capture log messages
        self.log_capture = StringIO()
        handler = logging.StreamHandler(self.log_capture)
        handler.setLevel(logging.DEBUG)
        self.logger.addHandler(handler)

    def tearDown(self):
        """Clean up logger handlers"""
        self.logger.handlers.clear()

    def test_single_cds_fragment_no_insertions_with_validation(self):
        """Single CDS fragment with successful validation"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = "ATG---CGATCG---TAA"
        raw_seq = "NNNNATGCGATCGTAANNNN"  # Contains the CDS
        cds_fragments = [{'start': 1, 'end': 19}]
        insertions_list = []
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "HA", offset, self.logger
        )

        self.assertEqual(result, VALID_CDS_PERGENE)
        self.assertTrue(valid)

    def test_single_cds_fragment_validation_failure(self):
        """Single CDS fragment with validation failure"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = "ATG---CGATCG---TAA"
        raw_seq = "NNNNTTTGGGCCCAAANNNN"  # Does NOT contain the CDS
        cds_fragments = [{'start': 1, 'end': 19}]
        insertions_list = []
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "HA", offset, self.logger
        )

        self.assertEqual(result, VALID_CDS_PERGENE)
        self.assertFalse(valid)
        # Check that error was logged
        log_output = self.log_capture.getvalue()
        self.assertIn("FAIL", log_output)

    def test_spliced_gene_both_fragments_valid(self):
        """Spliced gene with both fragments valid"""
        from create_unaligned_coding_seqs import extract_gene_cds

        # Fragment 1: positions 1-30
        # Fragment 2: positions 503-514
        aligned_seq = NEP_FRAGMENT1_PERGENE + ("N" * 472) + NEP_FRAGMENT2_PERGENE
        raw_seq = NEP_FRAGMENT1_PERGENE + "NNNNN" + NEP_FRAGMENT2_PERGENE  # Contains both fragments
        cds_fragments = [
            {'start': 1, 'end': 30},
            {'start': 503, 'end': 514}
        ]
        insertions_list = []
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "NEP", offset, self.logger
        )

        expected = NEP_FRAGMENT1_PERGENE + NEP_FRAGMENT2_PERGENE
        self.assertEqual(result, expected)
        self.assertTrue(valid)

    def test_spliced_gene_first_fragment_invalid(self):
        """Spliced gene with first fragment validation failure"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = NEP_FRAGMENT1_PERGENE + ("N" * 472) + NEP_FRAGMENT2_PERGENE
        raw_seq = "NNNNN" + NEP_FRAGMENT2_PERGENE  # Missing first fragment
        cds_fragments = [
            {'start': 1, 'end': 30},
            {'start': 503, 'end': 514}
        ]
        insertions_list = []
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "NEP", offset, self.logger
        )

        self.assertFalse(valid)
        log_output = self.log_capture.getvalue()
        self.assertIn("fragment", log_output.lower())
        self.assertIn("FAIL", log_output)

    def test_spliced_gene_second_fragment_invalid(self):
        """Spliced gene with second fragment validation failure"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = NEP_FRAGMENT1_PERGENE + ("N" * 472) + NEP_FRAGMENT2_PERGENE
        raw_seq = NEP_FRAGMENT1_PERGENE + "NNNNN"  # Missing second fragment
        cds_fragments = [
            {'start': 1, 'end': 30},
            {'start': 503, 'end': 514}
        ]
        insertions_list = []
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "NEP", offset, self.logger
        )

        self.assertFalse(valid)
        log_output = self.log_capture.getvalue()
        self.assertIn("fragment", log_output.lower())
        self.assertIn("FAIL", log_output)

    def test_case_insensitive_validation(self):
        """Fragment validation should be case-insensitive"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = "atg---cgatcg---taa"
        raw_seq = "NNNNATGCGATCGTAANNNN"  # Uppercase
        cds_fragments = [{'start': 1, 'end': 19}]
        insertions_list = []
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "HA", offset, self.logger
        )

        self.assertTrue(valid)

    def test_codon_insertion_with_validation(self):
        """CDS with codon insertion should validate correctly"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = "ATG---CGATCG---TAA"
        # After insertion at aligned position 12 (9th ungapped nucleotide 'G'):
        # Insert AAA after position 12: ATGCGATCGAAATAA (15 bp)
        # Using AAA to avoid ambiguity (left='G', right='T', neither matches 'A')
        raw_seq = "NNNNATGCGATCGAAATAANNNN"
        cds_fragments = [{'start': 1, 'end': 18}]
        insertions_list = [(12, "AAA")]  # Insert AAA after aligned position 12
        offset = 0  # No offset

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "HA", offset, self.logger
        )

        expected = "ATGCGATCGAAATAA"  # 15 bp (12 original + 3 inserted)
        self.assertEqual(result, expected)
        self.assertTrue(valid)

    def test_raw_seq_none_raises_error(self):
        """Function should raise ValueError if raw_seq is None"""
        from create_unaligned_coding_seqs import extract_gene_cds

        aligned_seq = "ATG---CGATCG---TAA"
        cds_fragments = [{'start': 1, 'end': 19}]
        insertions_list = []
        offset = 0  # No offset

        with self.assertRaises(ValueError) as context:
            extract_gene_cds(
                aligned_seq, cds_fragments, insertions_list,
                "seq1", None, "HA", offset, self.logger
            )

        self.assertIn("raw_seq is required", str(context.exception))

    def test_spliced_gene_with_offset(self):
        """Test spliced gene (NEP-like) with coordinate transformation"""
        from create_unaligned_coding_seqs import extract_gene_cds

        # Realistic NS segment: CDS starts at position 27 in original reference
        # Fragment 1: original coords 27-56 (30bp)
        # Fragment 2: original coords 529-864 (336bp)
        # offset = 26 (27 - 1)
        # After transformation:
        # Fragment 1: curated coords 1-30
        # Fragment 2: curated coords 503-838

        offset = 26
        # Curated MSA (30bp + 472bp gap + 336bp)
        aligned_seq = NEP_FRAGMENT1_PERGENE + ("N" * 472) + NEP_FRAGMENT2_PERGENE
        raw_seq = NEP_FRAGMENT1_PERGENE + "NNNNN" + NEP_FRAGMENT2_PERGENE

        cds_fragments = [
            {'start': 27, 'end': 56},    # Original coords
            {'start': 529, 'end': 864}   # Original coords
        ]
        insertions_list = []

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "seq1", raw_seq, "NEP", offset, self.logger
        )

        expected = NEP_FRAGMENT1_PERGENE + NEP_FRAGMENT2_PERGENE
        self.assertEqual(result, expected)
        self.assertTrue(valid)

    def test_full_workflow_with_coordinate_transformation(self):
        """Integration test: full CDS extraction with coordinate transformation"""
        from create_unaligned_coding_seqs import extract_gene_cds

        # Simulate realistic scenario:
        # - Original reference has 5' UTR (positions 1-57)
        # - CDS at positions 58-264 in original reference (207bp)
        # - Curated MSA trimmed to positions 1-207
        # - offset = 57

        offset = 57

        # Curated MSA: CDS only (no UTR) - 12bp with gaps
        # The alignment has gaps where the insertion was removed
        aligned_seq = "ATG---CGATCG---TAA"
        # Raw sequence contains the insertion that was removed during alignment
        raw_seq = "NNNNATGCGATCGAAATAANNNN"  # Contains "AAA" that will be re-inserted

        # CDS in original coords: 58-76 (19 positions in aligned seq)
        # In curated coords: 1-19 (after transformation)
        cds_fragments = [{'start': 58, 'end': 76}]

        # Insertion at position 70 in original reference
        # Should transform to position 13 in curated space (70 - 57)
        insertions_list = [(70, "AAA")]

        result, valid = extract_gene_cds(
            aligned_seq, cds_fragments, insertions_list,
            "test_seq", raw_seq, "HA", offset, self.logger
        )

        # Validate result
        self.assertTrue(valid)
        self.assertIn("AAA", result)  # Insertion should be present
        self.assertTrue(len(result) % 3 == 0)  # Should be in frame
        self.assertEqual(len(result), 15)  # 12 original + 3 inserted


class TestValidateCds(unittest.TestCase):
    """Tests for validate_cds() function"""

    def setUp(self):
        """Set up logger for tests"""
        self.logger = logging.getLogger('test_logger_validate')
        self.logger.setLevel(logging.DEBUG)
        self.log_capture = StringIO()
        handler = logging.StreamHandler(self.log_capture)
        handler.setLevel(logging.DEBUG)
        self.logger.addHandler(handler)

    def tearDown(self):
        """Clean up logger handlers"""
        self.logger.handlers.clear()

    def test_valid_cds_all_checks_pass(self):
        """Valid CDS with TAA stop codon"""
        from create_unaligned_coding_seqs import validate_cds

        result = validate_cds(VALID_CDS_PERGENE, "HA", "seq1", self.logger)
        self.assertTrue(result)

    def test_valid_cds_with_tag_stop(self):
        """Valid CDS with TAG stop codon"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "ATGCGATCGTAG"
        result = validate_cds(cds, "HA", "seq1", self.logger)
        self.assertTrue(result)

    def test_valid_cds_with_tga_stop(self):
        """Valid CDS with TGA stop codon"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "ATGCGATCGTGA"
        result = validate_cds(cds, "HA", "seq1", self.logger)
        self.assertTrue(result)

    def test_invalid_frame_remainder_one(self):
        """CDS with length not divisible by 3 (remainder 1)"""
        from create_unaligned_coding_seqs import validate_cds

        result = validate_cds(INVALID_FRAME_CDS_PERGENE, "HA", "seq1", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("not divisible by 3", log_output)

    def test_invalid_frame_remainder_two(self):
        """CDS with length not divisible by 3 (remainder 2)"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "ATGCGATCGTAAGG"  # 14 bp
        result = validate_cds(cds, "HA", "seq1", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("not divisible by 3", log_output)

    def test_missing_start_codon(self):
        """CDS without ATG start codon"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "TTGCGATCGTAA"
        result = validate_cds(cds, "HA", "seq1", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("does not start with ATG", log_output)

    def test_missing_stop_codon(self):
        """CDS without proper stop codon"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "ATGCGATCGTTA"  # Ends with TTA, not a stop codon
        result = validate_cds(cds, "HA", "seq1", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("does not end with stop codon", log_output)

    def test_too_short_sequence(self):
        """CDS too short (only start codon)"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "ATG"
        result = validate_cds(cds, "HA", "seq1", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("too short", log_output)

    def test_case_insensitive_validation(self):
        """Validation should be case-insensitive"""
        from create_unaligned_coding_seqs import validate_cds

        cds = "atgcgatcgtaa"  # Lowercase
        result = validate_cds(cds, "HA", "seq1", self.logger)
        self.assertTrue(result)

    def test_empty_sequence(self):
        """Empty sequence should fail validation"""
        from create_unaligned_coding_seqs import validate_cds

        cds = ""
        result = validate_cds(cds, "HA", "seq1", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("too short", log_output)

    def test_logging_contains_gene_name_and_seq_id(self):
        """Log messages should include gene name and sequence ID"""
        from create_unaligned_coding_seqs import validate_cds

        result = validate_cds(INVALID_FRAME_CDS_PERGENE, "HA", "test_seq", self.logger)

        self.assertFalse(result)
        log_output = self.log_capture.getvalue()
        self.assertIn("test_seq|HA", log_output)


if __name__ == "__main__":
    unittest.main()
