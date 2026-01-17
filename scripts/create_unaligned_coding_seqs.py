"""
Script to create unaligned coding sequences from curated aligned coding sequences.

This script takes curated aligned coding sequences and converts them back to unaligned
form by:
1. Removing deletion characters (gaps: '-')
2. Re-inserting nucleotides that were removed during the initial alignment step
"""

import argparse
import lzma
import logging
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import setup_logging, sanitize_id, get_coding_region_coords


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Create unaligned coding sequences from curated aligned coding sequences'
    )
    parser.add_argument('--curated-msa', required=True,
                        help='Input curated aligned MSA file (xz compressed)')
    parser.add_argument('--tsv', required=True,
                        help='Nextclade TSV file with insertion information (xz compressed)')
    parser.add_argument('--gff', required=True,
                        help='GFF file with gene annotation')
    parser.add_argument('--output', required=True,
                        help='Output file for unaligned coding sequences (xz compressed)')
    parser.add_argument('--raw-sequences', required=True,
                        help='Original unaligned raw sequences file for validation (xz compressed)')
    return parser.parse_args()


def parse_insertions_from_tsv(tsv_file, min_start, max_end):
    """
    Parse insertion information from Nextclade TSV file.

    Args:
        tsv_file: Path to Nextclade TSV file (xz compressed)
        min_start: Start position of coding region (1-based, inclusive)
        max_end: End position of coding region (1-based, inclusive)

    Returns:
        dict: {sanitized_seq_id: [(adjusted_position, nucleotides), ...]}
              Positions are adjusted for coding region and sorted in descending order
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Parsing insertions from {tsv_file}")

    # Read TSV file using pandas
    df = pd.read_csv(tsv_file, sep='\t', compression='xz', usecols=['seqName', 'insertions'])

    insertions_dict = {}
    total_insertions = 0
    filtered_insertions = 0

    # Iterate through rows
    for _, row in df.iterrows():
        seq_name = row['seqName']
        insertions_str = row['insertions']

        # Sanitize the sequence ID
        sanitized_id = sanitize_id(seq_name)

        # Parse insertions if present (handle NaN values from pandas)
        if pd.notna(insertions_str) and str(insertions_str).strip():
            insertions = []

            for insertion in str(insertions_str).split(','):
                insertion = insertion.strip()
                if not insertion:
                    continue

                try:
                    pos_str, nucs = insertion.split(':')
                    pos = int(pos_str)
                    total_insertions += 1

                    # Filter insertions based on coding region
                    # Insertion at position X means "insert AFTER position X"
                    # Keep if min_start <= X < max_end
                    if pos < min_start:
                        # Before coding region
                        logger.debug(f"{sanitized_id}: Ignoring insertion at position {pos} (before CDS)")
                        filtered_insertions += 1
                        continue
                    elif pos >= max_end:
                        # After coding region
                        logger.debug(f"{sanitized_id}: Ignoring insertion at position {pos} (after CDS)")
                        filtered_insertions += 1
                        continue

                    # Adjust position for coding region slice
                    # Original position X -> adjusted position (X - min_start + 1)
                    adjusted_pos = pos - min_start + 1
                    insertions.append((adjusted_pos, nucs))

                except ValueError as e:
                    logger.error(f"Could not parse insertion '{insertion}' for {seq_name}: {e}")
                    raise ValueError(f"Failed to parse insertion '{insertion}' for {seq_name}: {e}")

            # Sort insertions in descending order by position
            # This allows us to insert from highest to lowest position,
            # avoiding the need to track cumulative offsets
            if insertions:
                insertions.sort(reverse=True)
                insertions_dict[sanitized_id] = insertions

    logger.info(f"Parsed {total_insertions} total insertions from TSV")
    logger.info(f"Filtered out {filtered_insertions} insertions outside coding region")
    logger.info(f"Kept {total_insertions - filtered_insertions} insertions within coding region")
    logger.info(f"Found insertions for {len(insertions_dict)} sequences")

    return insertions_dict


def remove_gaps(seq_str):
    """
    Remove all gap characters from a sequence.

    Args:
        seq_str: Sequence string

    Returns:
        Sequence string with gaps removed
    """
    return seq_str.replace('-', '')


def insert_nucleotides(ungapped_seq, insertions):
    """
    Insert nucleotides at specified positions in an ungapped sequence.

    Args:
        ungapped_seq: Ungapped sequence string
        insertions: List of (position, nucleotides) tuples, sorted in descending order
                   Positions are 1-based, meaning "insert AFTER position X"

    Returns:
        Sequence string with insertions added
    """
    seq = ungapped_seq

    # insertions should already be sorted in descending order
    for pos, nucs in insertions:
        # Position is 1-based, meaning "insert AFTER position pos"
        # In Python 0-based indexing: seq[:pos] + nucs + seq[pos:]
        seq = seq[:pos] + nucs + seq[pos:]

    return seq


def validate_against_raw_sequences(unaligned_records, raw_sequences_file):
    """
    Validate that each unaligned coding sequence is a substring of the original raw sequence.

    Args:
        unaligned_records: List of SeqRecord objects with unaligned coding sequences
        raw_sequences_file: Path to original raw sequences file (xz compressed)

    Returns:
        tuple: (num_validated, num_failed)
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Validating unaligned sequences against original raw sequences from {raw_sequences_file}")

    # Load raw sequences and create a dictionary
    raw_seqs_dict = {}
    with lzma.open(raw_sequences_file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # Sanitize the ID to match curated sequences
            sanitized_id = sanitize_id(record.id)
            # Convert to uppercase for case-insensitive comparison
            raw_seqs_dict[sanitized_id] = str(record.seq).upper()

    logger.info(f"Loaded {len(raw_seqs_dict)} raw sequences")

    # Validate each unaligned sequence
    num_validated = 0
    num_failed = 0
    num_skipped_reference = 0

    for idx, record in enumerate(unaligned_records):
        seq_id = record.id
        unaligned_seq = str(record.seq).upper()

        if seq_id not in raw_seqs_dict:
            # Skip only the first sequence (reference)
            if idx == 0:
                num_skipped_reference += 1
                logger.info(f"{seq_id}: Skipping validation (reference sequence)")
                continue
            else:
                # Raise error for any other missing sequence
                logger.error(f"{seq_id}: Not found in raw sequences file")
                raise ValueError(f"{seq_id}: Not found in raw sequences file")

        raw_seq = raw_seqs_dict[seq_id]

        # Check if unaligned coding sequence is a substring of raw sequence
        if unaligned_seq in raw_seq:
            num_validated += 1
            logger.debug(f"{seq_id}: PASS - Unaligned coding sequence found in raw sequence")
        else:
            num_failed += 1
            logger.error(f"{seq_id}: FAIL - Unaligned coding sequence NOT found in raw sequence")
            logger.error(f"  Unaligned length: {len(unaligned_seq)}, Raw length: {len(raw_seq)}")

    logger.info(f"Validation results:")
    logger.info(f"  Validated: {num_validated}")
    logger.info(f"  Failed: {num_failed}")
    logger.info(f"  Skipped (reference): {num_skipped_reference}")

    if num_failed > 0:
        logger.error(f"VALIDATION FAILED: {num_failed} sequences are not substrings of their raw sequences")

    return num_validated, num_failed


def main():
    args = parse_args()
    logger = setup_logging()

    # Get coding region coordinates
    min_start, max_end = get_coding_region_coords(args.gff)

    # Parse insertions from TSV
    insertions_dict = parse_insertions_from_tsv(args.tsv, min_start, max_end)

    # Load curated aligned sequences
    logger.info(f"Reading curated aligned sequences from {args.curated_msa}")
    with lzma.open(args.curated_msa, 'rt') as handle:
        curated_records = list(SeqIO.parse(handle, 'fasta'))

    logger.info(f"Loaded {len(curated_records)} curated aligned sequences")

    # Process each sequence
    unaligned_records = []
    sequences_with_insertions = 0
    total_insertions_added = 0

    for record in curated_records:
        seq_id = record.id
        aligned_seq = str(record.seq)

        # Check if this sequence has insertions
        if seq_id in insertions_dict:
            insertions = insertions_dict[seq_id]
            sequences_with_insertions += 1
            total_insertions_added += len(insertions)

            # Insert nucleotides into aligned sequence first
            seq_with_insertions = insert_nucleotides(aligned_seq, insertions)

            # Then remove gaps
            final_seq = remove_gaps(seq_with_insertions)

            logger.info(f"{seq_id}: Added {len(insertions)} insertion(s), "
                       f"length {len(aligned_seq)} -> {len(final_seq)}")
        else:
            # No insertions for this sequence, just remove gaps
            final_seq = remove_gaps(aligned_seq)
            logger.debug(f"{seq_id}: No insertions, length {len(aligned_seq)} -> {len(final_seq)}")

        # Create new SeqRecord with unaligned sequence
        unaligned_record = SeqRecord(
            Seq(final_seq),
            id=seq_id,
            description=seq_id
        )
        unaligned_records.append(unaligned_record)

    # Write output
    logger.info(f"Writing {len(unaligned_records)} unaligned sequences to {args.output}")
    with lzma.open(args.output, 'wt') as handle:
        SeqIO.write(unaligned_records, handle, 'fasta')

    # Summary statistics
    logger.info(f"Summary:")
    logger.info(f"  Total sequences processed: {len(unaligned_records)}")
    logger.info(f"  Sequences with insertions: {sequences_with_insertions}")
    logger.info(f"  Total insertions added: {total_insertions_added}")

    # Validate against raw sequences
    logger.info(f"Running validation against raw sequences")
    num_validated, num_failed = validate_against_raw_sequences(unaligned_records, args.raw_sequences)

    if num_failed > 0:
        logger.error(f"Validation failed for {num_failed} sequences")
        return 1
    else:
        logger.info(f"All sequences validated successfully")

    return 0


if __name__ == "__main__":
    exit(main())
