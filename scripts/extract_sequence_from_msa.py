"""
Extract a specific sequence from a compressed MSA file by sequence name.
"""
import argparse
import lzma
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msa", required=True, help="Input MSA file (.fasta.xz)")
    parser.add_argument("--sequence-name", required=True, help="Name of sequence to extract")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    # Read compressed MSA
    with lzma.open(args.msa, 'rt') as f:
        sequences = SeqIO.parse(f, 'fasta')

        # Find the target sequence
        for record in sequences:
            if record.id == args.sequence_name:
                # Write to output
                with open(args.output, 'w') as out:
                    SeqIO.write(record, out, 'fasta')
                print(f"Extracted sequence: {args.sequence_name}")
                return

    raise ValueError(f"Sequence '{args.sequence_name}' not found in MSA")

if __name__ == "__main__":
    main()
