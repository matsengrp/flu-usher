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
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import setup_logging, sanitize_id, extract_all_genes_and_cds, group_cds_by_gene


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
                        help='Original GFF file with reference coordinates (NOT curated)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for per-gene FASTA files')
    parser.add_argument('--raw-sequences', required=True,
                        help='Original unaligned raw sequences file for validation (xz compressed)')
    return parser.parse_args()


def parse_insertions_from_tsv(tsv_file):
    """
    Parse insertion information from Nextclade TSV file.
    Parses all insertions without position filtering.

    Args:
        tsv_file: Path to Nextclade TSV file (xz compressed)

    Returns:
        dict: {sanitized_seq_id: [(position, nucleotides), ...]}
              Positions are 1-based and sorted in descending order
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Parsing insertions from {tsv_file}")

    # Read TSV file using pandas
    df = pd.read_csv(tsv_file, sep='\t', compression='xz', usecols=['seqName', 'insertions'])

    insertions_dict = {}
    total_insertions = 0

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
                    insertions.append((pos, nucs))

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


def filter_insertions_for_cds(insertions_list, cds_start_original, cds_end_original, offset):
    """
    Filter insertions that fall within CDS boundaries and transform to curated coordinates.

    Args:
        insertions_list: List of (pos, nucs) tuples in descending order (original reference coords)
        cds_start_original: CDS start position in original reference (1-based, inclusive)
        cds_end_original: CDS end position in original reference (1-based, inclusive)
        offset: Offset to transform from original to curated coordinates (min_start - 1)

    Returns:
        List of (adjusted_pos, nucs) tuples in curated coordinate space
    """
    filtered = []
    for pos, nucs in insertions_list:
        # Insertion at position X means "insert AFTER position X"
        # Keep if cds_start_original <= X < cds_end_original (in original coords)
        if cds_start_original <= pos < cds_end_original:
            # Transform to curated coordinate space
            # Curated position = original position - offset
            curated_pos = pos - offset
            filtered.append((curated_pos, nucs))

    # Keep descending order
    filtered.sort(reverse=True)
    return filtered


def extract_cds_from_aligned(aligned_seq, cds_start_original, cds_end_original, offset):
    """
    Extract CDS region from aligned sequence (curated MSA).

    Args:
        aligned_seq: Aligned sequence string (with gaps) from curated MSA
        cds_start_original: CDS start in original reference coordinates (1-based, inclusive)
        cds_end_original: CDS end in original reference coordinates (1-based, inclusive)
        offset: Offset to transform from original to curated coordinates (min_start - 1)

    Returns:
        CDS subsequence (still contains gaps)
    """
    # Transform original coordinates to curated coordinate space
    cds_start_curated = cds_start_original - offset
    cds_end_curated = cds_end_original - offset

    # Convert to 0-based indexing for slicing
    return aligned_seq[cds_start_curated-1:cds_end_curated]


def extract_gene_cds(aligned_seq, cds_fragments, insertions_list, seq_id, raw_seq, gene_name, offset, logger):
    """
    Extract complete CDS for a gene (handles spliced genes).
    Validates each fragment against raw sequence as it's extracted.

    Args:
        aligned_seq: Full aligned sequence from curated MSA
        cds_fragments: List of CDS feature dicts sorted by start position (original coords)
        insertions_list: List of (pos, nucs) for this sequence (original coords)
        seq_id: Sequence identifier
        raw_seq: Raw (unaligned) sequence for validation (REQUIRED)
        gene_name: Gene name for logging
        offset: Offset to transform from original to curated coordinates (min_start - 1)
        logger: Logger instance

    Returns:
        tuple: (complete_cds, all_fragments_valid)
            - complete_cds: Unaligned CDS sequence (concatenated if spliced)
            - all_fragments_valid: True if all fragments validated against raw sequence

    Raises:
        ValueError: If raw_seq is None
    """
    if raw_seq is None:
        raise ValueError(f"{seq_id}|{gene_name}: raw_seq is required for validation")

    unaligned_fragments = []
    all_fragments_valid = True

    for idx, cds in enumerate(cds_fragments):
        # Extract aligned region (transform from original to curated coords)
        aligned_cds = extract_cds_from_aligned(aligned_seq, cds['start'], cds['end'], offset)

        # Filter insertions for this CDS fragment (transform coords)
        cds_insertions = filter_insertions_for_cds(
            insertions_list, cds['start'], cds['end'], offset
        )

        # Apply insertions to aligned sequence (BEFORE gap removal)
        if cds_insertions:
            aligned_cds = insert_nucleotides(aligned_cds, cds_insertions)
            logger.debug(f"{seq_id}: Applied {len(cds_insertions)} insertions to CDS fragment {idx+1}")

        # Remove gaps
        unaligned_cds = remove_gaps(aligned_cds)
        unaligned_fragments.append(unaligned_cds)

        logger.debug(f"{seq_id}: CDS fragment {idx+1} ({cds['start']}-{cds['end']} original coords): {len(unaligned_cds)} bp")

        # Validate fragment against raw sequence
        if unaligned_cds.upper() in raw_seq.upper():
            logger.debug(f"{seq_id}|{gene_name} fragment {idx+1}: PASS")
        else:
            all_fragments_valid = False
            logger.error(
                f"{seq_id}|{gene_name} fragment {idx+1}: FAIL - "
                f"Fragment not found in raw sequence "
                f"(Fragment: {len(unaligned_cds)} bp, Raw: {len(raw_seq)} bp)"
            )

    # Concatenate all fragments for spliced genes
    complete_cds = ''.join(unaligned_fragments)

    if len(unaligned_fragments) > 1:
        logger.info(f"{seq_id}: Concatenated {len(unaligned_fragments)} CDS fragments → {len(complete_cds)} bp")

    return complete_cds, all_fragments_valid


def validate_cds(cds_seq, gene_name, seq_id, logger):
    """
    Validate that CDS is biologically correct:
    1. Length is multiple of 3
    2. Starts with ATG (start codon)
    3. Ends with stop codon (TAA, TAG, or TGA)

    Args:
        cds_seq: CDS sequence string
        gene_name: Gene name for logging
        seq_id: Sequence ID for logging
        logger: Logger instance

    Returns:
        bool: True if valid, False otherwise
    """
    length = len(cds_seq)

    # Check length is multiple of 3
    if length % 3 != 0:
        remainder = length % 3
        logger.warning(
            f"{seq_id}|{gene_name}: CDS length {length} not divisible by 3 "
            f"(remainder: {remainder}) - excluding from output"
        )
        return False

    # Check minimum length (at least start codon + stop codon = 6 bp)
    if length < 6:
        logger.warning(
            f"{seq_id}|{gene_name}: CDS length {length} too short "
            f"(minimum 6 bp for start + stop codon) - excluding from output"
        )
        return False

    # Check starts with ATG
    cds_upper = cds_seq.upper()
    if not cds_upper.startswith('ATG'):
        logger.warning(
            f"{seq_id}|{gene_name}: CDS does not start with ATG "
            f"(starts with {cds_seq[:3]}) - excluding from output"
        )
        return False

    # Check ends with stop codon (TAA, TAG, or TGA)
    stop_codons = {'TAA', 'TAG', 'TGA'}
    last_codon = cds_upper[-3:]
    if last_codon not in stop_codons:
        logger.warning(
            f"{seq_id}|{gene_name}: CDS does not end with stop codon "
            f"(ends with {last_codon}) - excluding from output"
        )
        return False

    return True


def main():
    args = parse_args()
    logger = setup_logging()

    # Load original GFF (with reference coordinates)
    logger.info(f"Reading original GFF from {args.gff}")
    features, (min_start, max_end) = extract_all_genes_and_cds(args.gff)
    genes = group_cds_by_gene(features)

    # Calculate offset for coordinate transformation
    # Curated MSA starts at position 1, but original reference starts at min_start
    offset = min_start - 1
    logger.info(f"Coordinate offset: {offset} (min_start={min_start}, max_end={max_end})")

    logger.info(f"Found {len(genes)} genes: {', '.join([g['gene_name'] for g in genes.values()])}")
    for gene_key, gene_data in genes.items():
        gene_name = gene_data['gene_name']
        num_cds = len(gene_data['cds_list'])
        logger.info(f"  {gene_name}: {num_cds} CDS feature(s)")

    # Parse insertions (in original reference coordinates)
    insertions_dict = parse_insertions_from_tsv(args.tsv)

    # Load curated MSA
    logger.info(f"Reading curated aligned sequences from {args.curated_msa}")
    with lzma.open(args.curated_msa, 'rt') as handle:
        curated_records = list(SeqIO.parse(handle, 'fasta'))
    logger.info(f"Loaded {len(curated_records)} curated aligned sequences")

    # Load raw sequences for validation
    logger.info(f"Reading raw sequences from {args.raw_sequences}")
    raw_seqs_dict = {}
    with lzma.open(args.raw_sequences, 'rt') as handle:
        for raw_record in SeqIO.parse(handle, 'fasta'):
            raw_seqs_dict[raw_record.id] = str(raw_record.seq)
    logger.info(f"Loaded {len(raw_seqs_dict)} raw sequences")

    # Prepare output data structures
    gene_records = {gene_data['gene_name']: [] for gene_data in genes.values()}
    validation_issues = []
    frame_issues = []
    num_skipped_no_raw = 0

    # Process each sequence
    for record in curated_records:
        seq_id = record.id
        aligned_seq = str(record.seq)

        # Get insertions for this sequence (in original reference coordinates)
        insertions_list = insertions_dict.get(seq_id, [])

        # Get raw sequence for validation - skip if not available
        raw_seq = raw_seqs_dict.get(seq_id)
        if raw_seq is None:
            num_skipped_no_raw += 1
            logger.debug(f"{seq_id}: Skipping (no raw sequence for validation)")
            continue

        # Extract CDS for each gene
        for gene_key, gene_data in genes.items():
            gene_name = gene_data['gene_name']
            cds_fragments = gene_data['cds_list']

            # Extract complete CDS and validate fragments
            cds_seq, fragments_valid = extract_gene_cds(
                aligned_seq,
                cds_fragments,
                insertions_list,
                seq_id,
                raw_seq,
                gene_name,
                offset,
                logger
            )

            # Check if fragments validated
            if not fragments_valid:
                validation_issues.append((seq_id, gene_name))
                logger.warning(f"{seq_id}|{gene_name}: Excluding due to fragment validation failure")
                continue  # Skip this gene for this sequence

            # Validate CDS: frame, start codon, stop codon
            if not validate_cds(cds_seq, gene_name, seq_id, logger):
                frame_issues.append((seq_id, gene_name, len(cds_seq)))
                continue  # Skip this gene for this sequence

            # Create SeqRecord
            unaligned_record = SeqRecord(
                Seq(cds_seq),
                id=seq_id,
                description=seq_id
            )
            gene_records[gene_name].append(unaligned_record)

    # Write separate output file for each gene
    os.makedirs(args.output_dir, exist_ok=True)

    for gene_name, records in gene_records.items():
        output_file = os.path.join(
            args.output_dir,
            f"curated_unaligned_{gene_name}.fasta.xz"
        )
        logger.info(f"Writing {len(records)} sequences for {gene_name} to {output_file}")
        with lzma.open(output_file, 'wt') as handle:
            SeqIO.write(records, handle, 'fasta')

    # Summary statistics
    logger.info("Summary:")
    logger.info(f"  Total sequences processed: {len(curated_records)}")
    logger.info(f"  Genes extracted: {len(genes)}")
    for gene_name, records in gene_records.items():
        logger.info(f"    {gene_name}: {len(records)} CDS sequences")

    if num_skipped_no_raw > 0:
        logger.info(f"  {num_skipped_no_raw} sequences skipped (no raw sequence for validation)")

    if validation_issues:
        logger.warning(f"  {len(validation_issues)} sequence-gene pairs excluded due to fragment validation failures")

    if frame_issues:
        logger.warning(f"  {len(frame_issues)} sequence-gene pairs excluded due to validation failures (frame, start codon, or stop codon)")

    # Return success (validation already done inline)
    logger.info("CDS extraction and validation complete")
    return 0


if __name__ == "__main__":
    exit(main())
