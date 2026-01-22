"""
Combined script to curate MSA and extract unaligned coding sequences.

This script combines the functionality of curate_msa.py and create_unaligned_coding_seqs.py
into a single unified pipeline for improved performance and consistency.

Key features:
1. Single-pass processing: Each sequence is filtered and validated in one pass
2. Stricter filtering: Sequences must pass ALL gene validations to be included
3. Counter-based logging: Tracks statistics instead of logging every failure
4. No intermediate files: Eliminates I/O overhead of curated_msa.fasta.xz
"""

import argparse
import lzma
import logging
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import (
    setup_logging,
    sanitize_id,
    extract_all_genes_and_cds,
    group_cds_by_gene
)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Curate MSA and extract unaligned coding sequences for flu-usher pipeline'
    )
    # Input files
    parser.add_argument('--input', required=True,
                        help='Input MSA file from Nextclade (xz compressed)')
    parser.add_argument('--gff', required=True,
                        help='Original GFF file with gene annotation')
    parser.add_argument('--tsv', required=True,
                        help='Nextclade TSV file with insertion information (xz compressed)')
    parser.add_argument('--raw-sequences', required=True,
                        help='Original unaligned raw sequences file for validation (xz compressed)')

    # Output directories
    parser.add_argument('--output-dir', required=True,
                        help='Directory for curated MSA and reference files')
    parser.add_argument('--output-coding-dir', required=True,
                        help='Directory for per-gene unaligned coding sequence files')

    # Quality filtering parameters
    parser.add_argument('--max_frac_gaps', type=float, default=0.05,
                        help='Maximum fraction of gaps allowed in a sequence')
    parser.add_argument('--max_frac_ambig', type=float, default=0.01,
                        help='Maximum fraction of ambiguous nucleotides allowed')
    parser.add_argument('--filter_duplicates', action='store_true',
                        help='Filter out duplicate sequences (keeping the first occurrence)')
    parser.add_argument('--replace_gaps_with_ref', action='store_true',
                        help='Replace gap characters (-) with reference nucleotides (first sequence)')

    return parser.parse_args()


def create_matching_gff_and_gtf(output_gff, output_gtf, features, min_start):
    """
    Create GFF and GTF files with adjusted coordinates for all features

    Args:
        output_gff: Path to output GFF file
        output_gtf: Path to output GTF file
        features: List of feature dictionaries
        min_start: Minimum start coordinate (will be adjusted to 1)
    """
    logger = logging.getLogger(__name__)

    # Calculate the total length of the region
    max_end_adjusted = max(f['end'] - min_start + 1 for f in features)

    # Create GFF content
    gff_lines = ["##gff-version 3\n"]
    gff_lines.append(f"##sequence-region Reference 1 {max_end_adjusted}\n")

    # Create GTF lines
    gtf_lines = []

    for i, feature in enumerate(features):
        # Adjust coordinates (min_start becomes 1)
        adj_start = feature['start'] - min_start + 1
        adj_end = feature['end'] - min_start + 1

        # GFF line - keep all original attributes, just change coordinates
        gff_line = f"Reference\tCurated\t{feature['type']}\t{adj_start}\t{adj_end}\t.\t{feature['strand']}\t{feature['phase']}\t{feature['attributes']}\n"
        gff_lines.append(gff_line)

        # GTF line - use original ID and Name for GTF attributes
        gtf_gene_id = feature['id'] if feature['id'] else feature['name']
        gtf_gene_name = feature['name'] if feature['name'] else feature['id']

        gtf_line = f"Reference\tCurated\t{feature['type']}\t{adj_start}\t{adj_end}\t.\t{feature['strand']}\t{feature['phase']}\t"
        gtf_line += f'gene_id "{gtf_gene_id}"; transcript_id "{gtf_gene_id}_transcript"; gene_name "{gtf_gene_name}";\n'
        gtf_lines.append(gtf_line)

    # Write GFF file
    with open(output_gff, 'w') as f:
        f.writelines(gff_lines)

    # Write GTF file
    with open(output_gtf, 'w') as f:
        f.writelines(gtf_lines)

    logger.info(f"Created GFF file with {len(features)} features")
    logger.info(f"Created GTF file with {len(features)} features")
    logger.info(f"Adjusted coordinates: region spans 1-{max_end_adjusted} (original: {min_start}-{max(f['end'] for f in features)})")


def slice_record(record, min_start, max_end):
    """
    Slice a sequence record to include only the region of interest
    and sanitize the sequence ID

    Args:
        record: BioPython SeqRecord object
        min_start: Start position of the region (1-based)
        max_end: End position of the region (1-based)

    Returns:
        New SeqRecord with the sliced sequence and sanitized ID
    """
    # Slice the sequence (adjust to 0-based indexing)
    sliced_seq = record.seq[min_start-1:max_end]

    # Sanitize the sequence ID
    sanitized_seq_id = sanitize_id(record.id)

    if ' ' in sanitized_seq_id:
        raise ValueError(f"Sanitized sequence ID contains spaces: {record.id}")

    # Create a new SeqRecord with the sliced sequence and sanitized ID
    new_record = SeqRecord(
        Seq(sliced_seq),
        id=sanitized_seq_id,
        description=sanitized_seq_id
    )
    return new_record


def get_ambiguous_chars(records):
    """
    Get a set of ambiguous characters in the sequences

    Args:
        records: List of SeqRecord objects

    Returns:
        Set of ambiguous characters
    """
    # Make a list of unique characters in the sequences
    unique_characters = set()
    for record in records:
        unique_characters.update(set(record.seq))

    # Make a set of ambiguous characters that are not A, C, G, T, or -
    unambiguous_characters = set("ACGT-")
    ambiguous_characters = unique_characters - unambiguous_characters

    return ambiguous_characters


def analyze_record(record, ambiguous_characters):
    """
    Calculate metrics for a sequence record

    Args:
        record: BioPython SeqRecord object
        ambiguous_characters: Set of ambiguous characters

    Returns:
        Dictionary with sequence metrics
    """
    seq = record.seq
    length = len(seq)
    num_gaps = seq.count('-')

    num_ambiguous = 0
    for char in ambiguous_characters:
        num_ambiguous += seq.count(char)

    return {
        'id': record.id,
        'length': length,
        'frac_gaps': num_gaps/length,
        'frac_ambiguous': num_ambiguous/length,
    }


def filter_and_yield_sequences(records, ambiguous_characters, max_frac_gaps, max_frac_ambig,
                                filter_duplicates=False, replace_gaps_with_ref=False):
    """
    Generator that yields filtered and modified sequences one at a time.

    Tracks filtering statistics internally instead of logging each failure.

    Args:
        records: List of SeqRecord objects
        ambiguous_characters: Set of ambiguous characters
        max_frac_gaps: Maximum fraction of gaps allowed
        max_frac_ambig: Maximum fraction of ambiguous nucleotides allowed
        filter_duplicates: Whether to filter out duplicate sequences (keeping first occurrence)
        replace_gaps_with_ref: Whether to replace gaps with reference nucleotides (first sequence)

    Yields:
        Tuple of (modified_record, quality_stats)
        quality_stats is a dict with cumulative counts of filter failures
    """
    logger = logging.getLogger(__name__)

    seen_sequences = set()

    # Initialize quality filter counters
    quality_stats = {
        'total_processed': 0,
        'failed_gap_fraction': 0,
        'failed_ambig_fraction': 0,
        'failed_terminal_gaps': 0,
        'failed_duplicate': 0,
        'passed_quality_filters': 0
    }

    # Get reference sequence if gap replacement is requested
    reference_seq = None
    if replace_gaps_with_ref:
        reference_seq = str(records[0].seq)
        # Assert that reference sequence has no gaps
        if '-' in reference_seq:
            raise ValueError("Reference sequence (first sequence) contains gaps")
        logger.info(f"Using {records[0].id} as reference for gap replacement")

    for record in records:
        quality_stats['total_processed'] += 1
        seq_str = str(record.seq)

        # Compute QC metrics on record
        metrics = analyze_record(record, ambiguous_characters)

        # Check if the first 3 or last 3 nucleotides are gaps
        has_terminal_gaps = ('---' == seq_str[:3]) or ('---' == seq_str[-3:])

        # Replace ambiguous characters with N
        if ambiguous_characters:
            for char in ambiguous_characters:
                seq_str = seq_str.replace(char, 'N')

        # Replace gaps with reference nucleotides if requested
        if replace_gaps_with_ref:
            seq_list = list(seq_str)
            for i, char in enumerate(seq_list):
                if char == '-':
                    seq_list[i] = reference_seq[i]
            seq_str = ''.join(seq_list)

        # Create a new record with the modified sequence
        if ambiguous_characters or replace_gaps_with_ref:
            modified_record = SeqRecord(
                Seq(seq_str),
                id=record.id,
                description=record.description
            )
        else:
            modified_record = record

        # Check for duplicate sequences if requested
        if filter_duplicates:
            if seq_str in seen_sequences:
                quality_stats['failed_duplicate'] += 1
                continue
            seen_sequences.add(seq_str)

        # Filter based on gap and ambiguous nucleotide content
        if has_terminal_gaps:
            quality_stats['failed_terminal_gaps'] += 1
            continue
        elif metrics['frac_gaps'] > max_frac_gaps:
            quality_stats['failed_gap_fraction'] += 1
            continue
        elif metrics['frac_ambiguous'] > max_frac_ambig:
            quality_stats['failed_ambig_fraction'] += 1
            continue

        # Passed all filters
        quality_stats['passed_quality_filters'] += 1
        yield modified_record, quality_stats


def parse_insertions_from_tsv(tsv_file):
    """
    Parse insertion information from Nextclade TSV file.

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
            if insertions:
                insertions.sort(reverse=True)
                insertions_dict[sanitized_id] = insertions

    logger.info(f"Parsed {total_insertions} total insertions")
    logger.info(f"Found insertions for {len(insertions_dict)} sequences")

    return insertions_dict


def remove_gaps(seq_str):
    """Remove all gap characters from a sequence"""
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
        seq = seq[:pos] + nucs + seq[pos:]

    return seq


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


def extract_gene_cds(aligned_seq, cds_fragments, insertions_list, seq_id, raw_seq,
                     gene_name, offset, logger, stats=None):
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
        stats: Statistics dictionary to update (optional, for tracking failures)

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
        if unaligned_cds.upper() not in raw_seq.upper():
            all_fragments_valid = False
            if stats is not None:
                stats['gene_validation_failures'][gene_name]['fragment_validation'] += 1
            logger.debug(
                f"{seq_id}|{gene_name} fragment {idx+1}: FAIL - "
                f"Fragment not found in raw sequence"
            )

    # Concatenate all fragments for spliced genes
    complete_cds = ''.join(unaligned_fragments)

    if len(unaligned_fragments) > 1:
        logger.debug(f"{seq_id}: Concatenated {len(unaligned_fragments)} CDS fragments → {len(complete_cds)} bp")

    return complete_cds, all_fragments_valid


def validate_cds(cds_seq, gene_name, seq_id, logger, stats=None):
    """
    Validate that CDS is biologically correct:
    1. Length is multiple of 3
    2. Starts with ATG (start codon)
    3. Ends with stop codon (TAA, TAG, or TGA)

    Args:
        cds_seq: CDS sequence string
        gene_name: Gene name for tracking stats
        seq_id: Sequence ID for logging
        logger: Logger instance
        stats: Statistics dictionary to update (optional, for tracking failures)

    Returns:
        bool: True if valid, False otherwise
    """
    length = len(cds_seq)

    # Check length is multiple of 3
    if length % 3 != 0:
        if stats is not None:
            stats['gene_validation_failures'][gene_name]['wrong_length'] += 1
        logger.warning(f"{seq_id}|{gene_name}: CDS length {length} not divisible by 3 (remainder: {length % 3}) - excluding from output")
        return False

    # Check minimum length
    if length < 6:
        if stats is not None:
            stats['gene_validation_failures'][gene_name]['too_short'] += 1
        logger.warning(f"{seq_id}|{gene_name}: CDS length {length} too short (minimum 6 bp for start + stop codon) - excluding from output")
        return False

    # Check starts with ATG
    cds_upper = cds_seq.upper()
    if not cds_upper.startswith('ATG'):
        if stats is not None:
            stats['gene_validation_failures'][gene_name]['missing_start_codon'] += 1
        logger.warning(f"{seq_id}|{gene_name}: CDS does not start with ATG (starts with {cds_seq[:3]}) - excluding from output")
        return False

    # Check ends with stop codon
    stop_codons = {'TAA', 'TAG', 'TGA'}
    last_codon = cds_upper[-3:]
    if last_codon not in stop_codons:
        if stats is not None:
            stats['gene_validation_failures'][gene_name]['missing_stop_codon'] += 1
        logger.warning(f"{seq_id}|{gene_name}: CDS does not end with stop codon (ends with {last_codon}) - excluding from output")
        return False

    return True


def filter_sequences(records, ambiguous_characters, max_frac_gaps, max_frac_ambig, logger,
                     filter_duplicates=False, replace_gaps_with_ref=False):
    """
    Backward-compatible wrapper around filter_and_yield_sequences for tests.
    Returns a list of filtered records.

    Args:
        records: List of SeqRecord objects
        ambiguous_characters: Set of ambiguous characters
        max_frac_gaps: Maximum fraction of gaps allowed
        max_frac_ambig: Maximum fraction of ambiguous nucleotides allowed
        logger: Logger object
        filter_duplicates: Whether to filter out duplicate sequences
        replace_gaps_with_ref: Whether to replace gaps with reference nucleotides

    Returns:
        List of filtered SeqRecord objects
    """
    filtered_records = []
    for record, _ in filter_and_yield_sequences(
        records, ambiguous_characters, max_frac_gaps, max_frac_ambig,
        filter_duplicates, replace_gaps_with_ref
    ):
        filtered_records.append(record)
    return filtered_records


def validate_against_raw_sequences(unaligned_records, raw_sequences_file):
    """
    Validate that each unaligned coding sequence is a substring of the original raw sequence.
    Backward-compatible function for tests.

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

    logger.info("=" * 80)
    logger.info("COMBINED CURATION AND CDS EXTRACTION PIPELINE")
    logger.info("=" * 80)

    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.output_coding_dir, exist_ok=True)

    # Define output file paths
    output_curated_msa = os.path.join(args.output_dir, "curated_msa.fasta.xz")
    output_reference_fasta = os.path.join(args.output_dir, "curated_reference.fasta")
    output_reference_txt = os.path.join(args.output_dir, "curated_reference.txt")
    output_reference_gff = os.path.join(args.output_dir, "curated_reference.gff")
    output_reference_gtf = os.path.join(args.output_dir, "curated_reference.gtf")

    # ========================================================================
    # SETUP PHASE: One-time operations
    # ========================================================================
    logger.info("SETUP PHASE: Loading reference data and metadata")

    # Extract all gene/CDS coordinates from the GFF file
    try:
        features, (min_start, max_end) = extract_all_genes_and_cds(args.gff)
        logger.info(f"Using coordinate range: {min_start}-{max_end}")

        # Create GFF and GTF files with all features, reindexed to start at 1
        create_matching_gff_and_gtf(
            output_reference_gff,
            output_reference_gtf,
            features,
            min_start
        )

    except Exception as e:
        logger.error(f"Failed to extract coordinates from GFF: {e}")
        return 1

    # Group CDS by gene (for handling spliced genes)
    genes = group_cds_by_gene(features)
    offset = min_start - 1
    logger.info(f"Coordinate offset: {offset}")
    logger.info(f"Found {len(genes)} genes: {', '.join([g['gene_name'] for g in genes.values()])}")

    # Parse insertion information from TSV
    insertions_dict = parse_insertions_from_tsv(args.tsv)

    # Load raw sequences for validation
    logger.info(f"Loading raw sequences from {args.raw_sequences}")
    raw_seqs_dict = {}
    with lzma.open(args.raw_sequences, 'rt') as handle:
        for raw_record in SeqIO.parse(handle, 'fasta'):
            sanitized_id = sanitize_id(raw_record.id)
            raw_seqs_dict[sanitized_id] = str(raw_record.seq).upper()
    logger.info(f"Loaded {len(raw_seqs_dict)} raw sequences")

    # Read original MSA
    logger.info(f"Reading alignment from {args.input}")
    with lzma.open(args.input, 'rt') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    logger.info(f"Read {len(records)} sequences")

    # Slice sequences to the region spanning all genes/CDS
    logger.info(f"Slicing sequences to region {min_start}-{max_end}")
    sliced_records = [slice_record(record, min_start, max_end) for record in records]

    # Identify ambiguous characters
    ambiguous_characters = get_ambiguous_chars(sliced_records)
    logger.info(f"Identified ambiguous characters: {ambiguous_characters}")

    # ========================================================================
    # PREPARE OUTPUT STRUCTURES AND STATISTICS
    # ========================================================================
    curated_records = []
    gene_records = {g['gene_name']: [] for g in genes.values()}

    # Initialize comprehensive filtering statistics
    stats = {
        'total_sequences': len(sliced_records),
        'passed_quality_filters': 0,
        'failed_cds_validation': 0,
        'added_to_curated_msa': 0,
        'gene_validation_failures': {
            g['gene_name']: {
                'fragment_validation': 0,
                'wrong_length': 0,
                'too_short': 0,
                'missing_start_codon': 0,
                'missing_stop_codon': 0
            } for g in genes.values()
        }
    }

    # ========================================================================
    # PROCESSING PHASE: Single-pass filtering and CDS extraction
    # ========================================================================
    logger.info("=" * 80)
    logger.info("PROCESSING PHASE: Filtering and extracting CDS")
    logger.info(f"Quality thresholds: max_gaps={args.max_frac_gaps}, max_ambig={args.max_frac_ambig}")
    logger.info(f"Options: filter_duplicates={args.filter_duplicates}, replace_gaps_with_ref={args.replace_gaps_with_ref}")
    logger.info("=" * 80)

    # Process sequences (single pass)
    quality_stats = None
    for idx, (curated_record, quality_stats) in enumerate(filter_and_yield_sequences(
        sliced_records,
        ambiguous_characters,
        args.max_frac_gaps,
        args.max_frac_ambig,
        args.filter_duplicates,
        args.replace_gaps_with_ref
    )):
        seq_id = curated_record.id
        stats['passed_quality_filters'] += 1

        # Always add reference sequence (first sequence) to curated MSA
        # but skip CDS extraction since it can't be validated against raw sequences
        if idx == 0:
            curated_records.append(curated_record)
            stats['added_to_curated_msa'] += 1
            logger.info(f"Reference sequence ({seq_id}) added to curated MSA (CDS extraction skipped)")
            continue

        # For all non-reference sequences, raw sequence must be available
        raw_seq = raw_seqs_dict.get(seq_id)
        if raw_seq is None:
            logger.error(f"{seq_id}: Raw sequence not found - this should not happen!")
            raise ValueError(f"Raw sequence not found for {seq_id}")

        aligned_seq = str(curated_record.seq)
        insertions_list = insertions_dict.get(seq_id, [])

        # Extract CDS for each gene and track validity
        all_genes_valid = True
        extracted_cds = {}

        for gene_key, gene_data in genes.items():
            gene_name = gene_data['gene_name']
            cds_fragments = gene_data['cds_list']

            # Extract CDS and validate fragments
            cds_seq, fragments_valid = extract_gene_cds(
                aligned_seq,
                cds_fragments,
                insertions_list,
                seq_id,
                raw_seq,
                gene_name,
                offset,
                logger,
                stats
            )

            # Validate biological correctness if fragments passed
            if fragments_valid and validate_cds(cds_seq, gene_name, seq_id, logger, stats):
                extracted_cds[gene_name] = cds_seq
            else:
                all_genes_valid = False
                break  # Stop processing this sequence

        # Only add to curated MSA if all genes are valid
        if all_genes_valid:
            curated_records.append(curated_record)
            stats['added_to_curated_msa'] += 1
            # Add all extracted CDS to per-gene collections
            for gene_name, cds_seq in extracted_cds.items():
                gene_records[gene_name].append(
                    SeqRecord(Seq(cds_seq), id=seq_id, description=seq_id)
                )
            logger.debug(f"{seq_id}: Passed all validations")
        else:
            stats['failed_cds_validation'] += 1

    # Update quality filter stats from generator
    if quality_stats:
        stats.update({
            'quality_filter_stats': quality_stats
        })

    # ========================================================================
    # OUTPUT PHASE: Write all output files
    # ========================================================================
    logger.info("=" * 80)
    logger.info("OUTPUT PHASE: Writing output files")
    logger.info("=" * 80)

    # Write curated MSA
    logger.info(f"Writing {len(curated_records)} curated sequences to {output_curated_msa}")
    with lzma.open(output_curated_msa, 'wt') as handle:
        SeqIO.write(curated_records, handle, 'fasta')

    # Write reference files
    if curated_records:
        logger.info(f"Writing reference sequence to {output_reference_fasta}")
        with open(output_reference_fasta, 'w') as handle:
            SeqIO.write([curated_records[0]], handle, 'fasta')

        logger.info(f"Writing reference sequence (no header) to {output_reference_txt}")
        with open(output_reference_txt, 'w') as handle:
            handle.write(str(curated_records[0].seq))
    else:
        logger.error("No sequences passed filtering, cannot write reference sequence")
        return 1

    # Write per-gene coding sequence files
    for gene_name, records in gene_records.items():
        output_file = os.path.join(
            args.output_coding_dir,
            f"curated_unaligned_{gene_name}.fasta.xz"
        )
        logger.info(f"Writing {len(records)} sequences for {gene_name} to {output_file}")
        with lzma.open(output_file, 'wt') as handle:
            SeqIO.write(records, handle, 'fasta')

    # ========================================================================
    # SUMMARY PHASE: Report comprehensive statistics
    # ========================================================================
    logger.info("=" * 80)
    logger.info("FILTERING SUMMARY")
    logger.info("=" * 80)

    # Overall statistics
    logger.info(f"Total sequences processed: {stats['total_sequences']:,}")

    # Quality filtering breakdown
    if 'quality_filter_stats' in stats:
        qs = stats['quality_filter_stats']
        failed_quality = qs['total_processed'] - qs['passed_quality_filters']
        logger.info(f"Passed quality filters: {qs['passed_quality_filters']:,}")
        logger.info(f"Failed quality filters: {failed_quality:,}")
        if failed_quality > 0:
            logger.info("  Quality filter breakdown:")
            if qs['failed_gap_fraction'] > 0:
                logger.info(f"    - Excessive gaps (>{args.max_frac_gaps*100}%): {qs['failed_gap_fraction']:,}")
            if qs['failed_ambig_fraction'] > 0:
                logger.info(f"    - Excessive ambiguous nucleotides (>{args.max_frac_ambig*100}%): {qs['failed_ambig_fraction']:,}")
            if qs['failed_terminal_gaps'] > 0:
                logger.info(f"    - Terminal gaps: {qs['failed_terminal_gaps']:,}")
            if qs['failed_duplicate'] > 0:
                logger.info(f"    - Duplicate sequences: {qs['failed_duplicate']:,}")

    # CDS validation statistics
    logger.info(f"Failed CDS validation: {stats['failed_cds_validation']:,}")
    logger.info(f"Added to curated MSA: {stats['added_to_curated_msa']:,}")

    # Per-gene validation failure breakdown
    logger.info("")
    logger.info("Per-gene CDS validation failures:")
    for gene_name, failure_counts in stats['gene_validation_failures'].items():
        total_gene_failures = sum(failure_counts.values())
        logger.info(f"  {gene_name}: {total_gene_failures:,} total failures")
        if total_gene_failures > 0:
            for failure_type, count in failure_counts.items():
                if count > 0:
                    logger.info(f"    - {failure_type.replace('_', ' ').title()}: {count:,}")

    # Retention rate
    retention_rate = (stats['added_to_curated_msa'] / stats['total_sequences']) * 100
    logger.info("")
    logger.info(f"Retention rate: {retention_rate:.1f}% of input sequences")
    logger.info("=" * 80)

    logger.info("Pipeline completed successfully!")
    return 0


if __name__ == "__main__":
    exit(main())
