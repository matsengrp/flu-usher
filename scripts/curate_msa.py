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
    parser.add_argument('--gene-name', required=True, help='Name of gene to extract from GFF')
    parser.add_argument('--max_frac_gaps', type=float, default=0.05, 
                        help='Maximum fraction of gaps allowed in a sequence')
    parser.add_argument('--max_frac_ambig', type=float, default=0.01, 
                        help='Maximum fraction of ambiguous nucleotides allowed')
    return parser.parse_args()

def extract_gene_coordinates(gff_file, gene_name):
    """
    Extract gene start and end positions from a GFF file
    
    Args:
        gff_file: Path to GFF file
        gene_name: Name of the gene to extract
    
    Returns:
        tuple: (start, end) coordinates (1-based), phase, strand, attributes
    """
    logger = logging.getLogger(__name__)
    
    # Regular expressions for parsing GFF
    re_gff = re.compile(r'([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)')
    
    logger.info(f"Extracting coordinates for gene '{gene_name}' from {gff_file}")
    
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
            
            # Look for gene or CDS features
            if feature_type.lower() in ('gene', 'cds'):
                # Check if gene name matches in the attributes
                if f"ID={gene_name}" in attributes or f"Name={gene_name}" in attributes or f"gene={gene_name}" in attributes:
                    logger.info(f"Found gene '{gene_name}' at position {start}-{end}")
                    return start, end, phase, strand, attributes
    
    # If we get here, no matching gene was found
    raise ValueError(f"Could not find gene '{gene_name}' in {gff_file}")

def create_matching_gff(output_file, gene_name, original_start, original_end, phase, strand, attributes):
    """
    Create a simplified GFF file with adjusted coordinates
    
    Args:
        output_file: Path to output GFF file
        gene_name: Name of the gene
        original_start: Original start position (1-based)
        original_end: Original end position (1-based)
        phase: Original phase
        strand: Original strand
        attributes: Original attributes
    """
    logger = logging.getLogger(__name__)
    
    # Calculate new coordinates (start at 1)
    new_length = original_end - original_start + 1
    
    # Create GFF content
    gff_content = "##gff-version 3\n"
    gff_content += f"##sequence-region Reference 1 {new_length}\n"
    gff_content += f"Reference\tCurated\tCDS\t1\t{new_length}\t.\t{strand}\t{phase}\tID={gene_name};Name={gene_name}\n"
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(gff_content)
    
    logger.info(f"Created simplified GFF file at {output_file}")
    logger.info(f"Adjusted coordinates: 1-{new_length} (original: {original_start}-{original_end})")

def slice_record(record, gene_start, gene_end):
    """
    Slice a sequence record to include only the region of interest
    and sanitize the sequence ID by replacing problematic characters
    
    Args:
        record: BioPython SeqRecord object
        gene_start: Start position of the gene (1-based)
        gene_end: End position of the gene (1-based)
    
    Returns:
        New SeqRecord with the sliced sequence and sanitized ID
    """
    # Slice the sequence (adjust to 0-based indexing)
    sliced_seq = record.seq[gene_start-1:gene_end]
    
    # Sanitize the sequence ID by replacing problematic characters with underscores
    sanitized_id = record.id
    for char in ['[', ']', '(', ')', ':', ';', ',', "'"]:
        sanitized_id = sanitized_id.replace(char, '_')
    
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

def filter_sequences(records, ambiguous_characters, max_frac_gaps, max_frac_ambig, logger):
    """
    Filter sequences based on gap and ambiguous nucleotide content
    
    Args:
        records: List of SeqRecord objects
        ambiguous_characters: Set of ambiguous characters
        max_frac_gaps: Maximum fraction of gaps allowed
        max_frac_ambig: Maximum fraction of ambiguous nucleotides allowed
        logger: Logger object
    
    Returns:
        List of filtered SeqRecord objects
    """
    filtered_records = []
    
    for record in records:
        metrics = analyze_record(record, ambiguous_characters)
        
        if metrics['frac_gaps'] < max_frac_gaps and metrics['frac_ambiguous'] < max_frac_ambig:
            filtered_records.append(record)
        else:
            logger.info(f"Filtered out {record.id}: gaps={metrics['frac_gaps']:.3f}, ambig={metrics['frac_ambiguous']:.3f}")
    
    return filtered_records

def main():
    args = parse_args()
    logger = setup_logging()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output file paths within the output directory
    output_curated_msa = os.path.join(args.output_dir, "curated_msa.fasta.xz")
    output_reference_fasta = os.path.join(args.output_dir, "curated_reference.fasta")
    output_reference_gff = os.path.join(args.output_dir, "curated_reference.gff")
    
    # Extract gene coordinates from the GFF file and write a simplified GFF file
    try:
        gene_start, gene_end, phase, strand, attributes = extract_gene_coordinates(args.gff, args.gene_name)
        logger.info(f"Using gene coordinates from GFF: {gene_start}-{gene_end}")
        
        # Create a GFF file reindexed to 1 for the gene of interest
        create_matching_gff(output_reference_gff, args.gene_name, gene_start, gene_end, phase, strand, attributes)
            
    except Exception as e:
        logger.error(f"Failed to extract coordinates for gene '{args.gene_name}': {e}")
        return 1
        
    # Read in the MSA
    logger.info(f"Reading alignment from {args.input}")
    with lzma.open(args.input, 'rt') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    
    logger.info(f"Read {len(records)} sequences")
    
    # Slice the sequences to the CDS of interest and santitize the IDs
    logger.info(f"Slicing sequences to region {gene_start}-{gene_end}")
    sliced_records = [slice_record(record, gene_start, gene_end) for record in records]
    
    # Filter sequences with too many gaps or ambiguous characters
    ambiguous_characters = get_ambiguous_chars(sliced_records)
    logger.info(f"Identified ambiguous characters: {ambiguous_characters}")
    
    logger.info(f"Filtering sequences (max gaps: {args.max_frac_gaps}, max ambig: {args.max_frac_ambig})")
    curated_records = filter_sequences(
        sliced_records, 
        ambiguous_characters, 
        args.max_frac_gaps, 
        args.max_frac_ambig, 
        logger
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
    else:
        logger.error("No sequences passed filtering, cannot write reference sequence")
        return 1
    
    logger.info(f"Retained {len(curated_records)}/{len(records)} sequences ({len(curated_records)/len(records)*100:.1f}%)")
    return 0

if __name__ == "__main__":
    exit(main())