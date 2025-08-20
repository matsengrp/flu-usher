"""
Randomize the order of sequences in an alignment, keeping the first sequence (reference) at the top.
"""

import argparse
import random
from Bio import SeqIO
import gzip
import lzma


def randomize_alignment(input_file, output_file, seed=None):
    """
    Randomize the order of sequences in an alignment while keeping the first sequence at the top.
    
    Args:
        input_file: Path to input FASTA file (can be compressed with .gz or .xz)
        output_file: Path to output FASTA file (can be compressed with .gz or .xz)
        seed: Random seed for reproducibility
    """
    if seed is not None:
        random.seed(seed)
    
    # Read all sequences
    sequences = []
    
    # Handle different compression formats
    if input_file.endswith('.gz'):
        handle = gzip.open(input_file, 'rt')
    elif input_file.endswith('.xz'):
        handle = lzma.open(input_file, 'rt')
    else:
        handle = open(input_file, 'r')
    
    try:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(record)
    finally:
        handle.close()
    
    if len(sequences) < 2:
        raise ValueError("Alignment must contain at least 2 sequences")
    
    # Keep the first sequence (reference) separate
    reference = sequences[0]
    other_sequences = sequences[1:]
    
    # Randomize the order of non-reference sequences
    random.shuffle(other_sequences)
    
    # Combine reference with randomized sequences
    randomized_sequences = [reference] + other_sequences
    
    # Write output
    if output_file.endswith('.gz'):
        handle = gzip.open(output_file, 'wt')
    elif output_file.endswith('.xz'):
        handle = lzma.open(output_file, 'wt')
    else:
        handle = open(output_file, 'w')
    
    try:
        SeqIO.write(randomized_sequences, handle, "fasta")
    finally:
        handle.close()
    
    print(f"Randomized {len(other_sequences)} sequences (keeping reference at top)")
    print(f"Random seed: {seed}")


def main():
    parser = argparse.ArgumentParser(description="Randomize sequence order in alignment")
    parser.add_argument("--input", "-i", required=True,
                        help="Input alignment file (FASTA format, can be .gz or .xz compressed)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output alignment file (FASTA format, can be .gz or .xz compressed)")
    parser.add_argument("--seed", "-s", type=int, default=None,
                        help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    randomize_alignment(args.input, args.output, args.seed)


if __name__ == "__main__":
    main()