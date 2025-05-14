#!/usr/bin/env python3
"""
Script to curate multiple sequence alignments (MSA) for flu-usher pipeline.
Based on curate_alignment.ipynb notebook.
"""

import argparse
import lzma
import logging
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
    parser.add_argument('--output', required=True, help='Output curated MSA file (xz compressed)')
    parser.add_argument('--gene_start', type=int, required=True, help='Gene start position (1-based)')
    parser.add_argument('--gene_end', type=int, required=True, help='Gene end position (1-based)')
    parser.add_argument('--max_frac_gaps', type=float, default=0.05, 
                        help='Maximum fraction of gaps allowed in a sequence')
    parser.add_argument('--max_frac_ambig', type=float, default=0.01, 
                        help='Maximum fraction of ambiguous nucleotides allowed')
    return parser.parse_args()

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
    for char in ['[', ']', '(', ')', ':', ';', ',']:
        sanitized_id = sanitized_id.replace(char, '_')
    
    # Create a new SeqRecord with the sliced sequence and sanitized ID
    new_record = SeqRecord(
        Seq(sliced_seq),
        id=sanitized_id,
        description=record.description
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
    
    logger.info(f"Reading alignment from {args.input}")
    with lzma.open(args.input, 'rt') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    
    logger.info(f"Read {len(records)} sequences")
    
    logger.info(f"Slicing sequences to region {args.gene_start}-{args.gene_end}")
    sliced_records = [slice_record(record, args.gene_start, args.gene_end) for record in records]
    
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
    
    logger.info(f"Writing {len(curated_records)} curated sequences to {args.output}")
    with lzma.open(args.output, 'wt') as handle:
        SeqIO.write(curated_records, handle, 'fasta')
    
    logger.info(f"Retained {len(curated_records)}/{len(records)} sequences ({len(curated_records)/len(records)*100:.1f}%)")

if __name__ == "__main__":
    main()