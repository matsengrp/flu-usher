"""
Script to curate multiple sequence alignments (MSA) for flu-usher pipeline.
"""

import argparse
import lzma
import logging
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Curate MSA for flu-usher pipeline')
    parser.add_argument('--input', required=True, help='Input MSA file (xz compressed)')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    parser.add_argument('--gff', required=True, help='GFF file with gene annotation')
    parser.add_argument('--max_frac_gaps', type=float, default=0.05, 
                        help='Maximum fraction of gaps allowed in a sequence')
    parser.add_argument('--max_frac_ambig', type=float, default=0.01, 
                        help='Maximum fraction of ambiguous nucleotides allowed')
    parser.add_argument('--filter_duplicates', action='store_true',
                        help='Filter out duplicate sequences (keeping the first occurrence)')
    parser.add_argument('--replace_gaps_with_ref', action='store_true',
                        help='Replace gap characters (-) with reference nucleotides (first sequence)')
    return parser.parse_args()

def extract_all_genes_and_cds(gff_file):
    """
    Extract information for all genes and CDS features from a GFF file
    
    Args:
        gff_file: Path to GFF file
    
    Returns:
        list: List of dictionaries containing gene/CDS information
        tuple: (min_start, max_end) coordinates for all features
    """
    logger = logging.getLogger(__name__)
    
    # Regular expressions for parsing GFF
    re_gff = re.compile(r'([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)')
    
    logger.info(f"Extracting all genes and CDS features from {gff_file}")
    
    features = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
                
            match = re_gff.match(line)
            if not match:
                continue
                
            seqid, source, feature_type, start, end, score, strand, phase, attributes = match.groups()
            
            # Convert to integers
            start, end = int(start), int(end)
            
            # Strip newline character from attributes
            attributes = attributes.strip()
            
            # Look for gene or CDS features
            if feature_type.lower() in ('gene', 'cds'):
                # Extract ID and Name from attributes
                feature_id = extract_attribute_value(attributes, 'ID')
                feature_name = extract_attribute_value(attributes, 'Name')
                
                # Use ID as fallback for name if Name is not present
                if not feature_name:
                    feature_name = feature_id
                
                # Use feature_type as fallback if neither ID nor Name is present
                if not feature_name:
                    feature_name = f"{feature_type}_feature"
                
                feature_info = {
                    'id': feature_id,
                    'name': feature_name,
                    'type': feature_type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'phase': phase,
                    'attributes': attributes,
                    'original_attributes': attributes  # Keep original for reference
                }
                features.append(feature_info)
                logger.info(f"Found {feature_type} ID='{feature_id}' Name='{feature_name}' at position {start}-{end}")
    
    if not features:
        raise ValueError(f"Could not find any genes or CDS features in {gff_file}")
    
    # Calculate overall min and max coordinates
    min_start = min(f['start'] for f in features)
    max_end = max(f['end'] for f in features)
    
    logger.info(f"Found {len(features)} features spanning {min_start}-{max_end}")
    return features, (min_start, max_end)

def extract_attribute_value(attributes, attr_name):
    """
    Extract a specific attribute value from GFF attributes string
    
    Args:
        attributes: GFF attributes string
        attr_name: Name of the attribute to extract (e.g., 'ID', 'Name')
    
    Returns:
        str: Attribute value or None if not found
    """
    attr_pattern = f"{attr_name}="
    if attr_pattern in attributes:
        start_idx = attributes.find(attr_pattern) + len(attr_pattern)
        end_idx = attributes.find(';', start_idx)
        if end_idx == -1:
            end_idx = len(attributes)
        return attributes[start_idx:end_idx]
    return None

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
    
    logger.info(f"Created GFF file at {output_gff} with {len(features)} features")
    logger.info(f"Created GTF file at {output_gtf} with {len(features)} features")
    logger.info(f"Adjusted coordinates: region spans 1-{max_end_adjusted} (original: {min_start}-{max(f['end'] for f in features)})")

def slice_record(record, min_start, max_end):
    """
    Slice a sequence record to include only the region of interest
    and sanitize the sequence ID by replacing problematic characters
    
    Args:
        record: BioPython SeqRecord object
        min_start: Start position of the region (1-based)
        max_end: End position of the region (1-based)
    
    Returns:
        New SeqRecord with the sliced sequence and sanitized ID
    """
    # Slice the sequence (adjust to 0-based indexing)
    sliced_seq = record.seq[min_start-1:max_end]
    
    # Sanitize the sequence ID by removing problematic characters
    sanitized_id = record.id
    for char in ['[', ']', '(', ')', ':', ';', ',', "'", '.']:
        sanitized_id = sanitized_id.replace(char, '')
    
    if ' ' in sanitized_id:
        raise ValueError(record.id)

    # Create a new SeqRecord with the sliced sequence and sanitized ID
    new_record = SeqRecord(
        Seq(sliced_seq),
        id=sanitized_id,
        description=sanitized_id
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

def filter_sequences(records, ambiguous_characters, max_frac_gaps, max_frac_ambig, logger, filter_duplicates=False, replace_gaps_with_ref=False):
    """
    Filter sequences based on gap and ambiguous nucleotide content
    
    Args:
        records: List of SeqRecord objects
        ambiguous_characters: Set of ambiguous characters
        max_frac_gaps: Maximum fraction of gaps allowed
        max_frac_ambig: Maximum fraction of ambiguous nucleotides allowed
        logger: Logger object
        filter_duplicates: Whether to filter out duplicate sequences (keeping first occurrence)
        replace_gaps_with_ref: Whether to replace gaps with reference nucleotides (first sequence)
    
    Returns:
        List of filtered SeqRecord objects
    """
    filtered_records = []
    seen_sequences = set()
    duplicate_count = 0
    
    # Get reference sequence if gap replacement is requested
    reference_seq = None
    if replace_gaps_with_ref:
        reference_seq = str(records[0].seq)
        # Assert that reference sequence has no gaps
        if '-' in reference_seq:
            raise ValueError("Reference sequence (first sequence) contains gaps")
        logger.info(f"Using {records[0].id} as reference for gap replacement")
    
    for record in records:
        seq_str = str(record.seq)
        
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
                duplicate_count += 1
                logger.info(f"Filtered out {record.id}: duplicate sequence")
                continue
            seen_sequences.add(seq_str)
        
        metrics = analyze_record(modified_record, ambiguous_characters)

        # Check if the first 3 or last 3 nucleotides are gaps
        has_terminal_gaps = ('---' == seq_str[:3]) or ('---' == seq_str[-3:])
        
        # Filter based on gap and ambiguous nucleotide content
        if has_terminal_gaps:
            logger.info(f"Filtered out {record.id}: has gaps in first 3 or last 3 nucleotides")
        elif metrics['frac_gaps'] <= max_frac_gaps and metrics['frac_ambiguous'] <= max_frac_ambig:
            filtered_records.append(modified_record)
        else:
            logger.info(f"Filtered out {record.id}: gaps={metrics['frac_gaps']:.3f}, ambig={metrics['frac_ambiguous']:.3f}")
    
    if filter_duplicates and duplicate_count > 0:
        logger.info(f"Filtered out {duplicate_count} duplicate sequences")
    
    return filtered_records

def main():
    args = parse_args()
    logger = setup_logging()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output file paths within the output directory
    output_curated_msa = os.path.join(args.output_dir, "curated_msa.fasta.xz")
    output_reference_fasta = os.path.join(args.output_dir, "curated_reference.fasta")
    output_reference_txt = os.path.join(args.output_dir, "curated_reference.txt")
    output_reference_gff = os.path.join(args.output_dir, "curated_reference.gff")
    output_reference_gtf = os.path.join(args.output_dir, "curated_reference.gtf")
    
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
        
    # Read in the MSA
    logger.info(f"Reading alignment from {args.input}")
    with lzma.open(args.input, 'rt') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    
    logger.info(f"Read {len(records)} sequences")
    
    # Slice the sequences to the region spanning all genes/CDS and sanitize the IDs
    logger.info(f"Slicing sequences to region {min_start}-{max_end}")
    sliced_records = [slice_record(record, min_start, max_end) for record in records]
    
    # Filter sequences with too many gaps or ambiguous characters
    ambiguous_characters = get_ambiguous_chars(sliced_records)
    logger.info(f"Identified ambiguous characters: {ambiguous_characters}")
    
    logger.info(f"Filtering sequences (max gaps: {args.max_frac_gaps}, max ambig: {args.max_frac_ambig}, filter duplicates: {args.filter_duplicates}, replace gaps with ref: {args.replace_gaps_with_ref})")
    curated_records = filter_sequences(
        sliced_records, 
        ambiguous_characters, 
        args.max_frac_gaps, 
        args.max_frac_ambig, 
        logger,
        args.filter_duplicates,
        args.replace_gaps_with_ref
    )
    
    # Write the curated sequences to an output FASTA file
    logger.info(f"Writing {len(curated_records)} curated sequences to {output_curated_msa}")
    with lzma.open(output_curated_msa, 'wt') as handle:
        SeqIO.write(curated_records, handle, 'fasta')
    
    # Write the first sequence to a separate uncompressed FASTA file as the reference
    if curated_records:
        logger.info(f"Writing reference sequence to {output_reference_fasta}")
        with open(output_reference_fasta, 'w') as handle:
            SeqIO.write([curated_records[0]], handle, 'fasta')
        
        # Write just the sequence (no header) to a text file
        logger.info(f"Writing reference sequence (no header) to {output_reference_txt}")
        with open(output_reference_txt, 'w') as handle:
            handle.write(str(curated_records[0].seq))
    else:
        logger.error("No sequences passed filtering, cannot write reference sequence")
        return 1
    
    logger.info(f"Retained {len(curated_records)}/{len(records)} sequences ({len(curated_records)/len(records)*100:.1f}%)")
    return 0

if __name__ == "__main__":
    exit(main())