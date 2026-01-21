"""
Create a samples file for matUtils extract with reference and new root names.
"""
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--new-root", required=True, help="Name of new root sample")
    parser.add_argument("--output", required=True, help="Output samples.txt file")
    args = parser.parse_args()

    # Extract reference name from FASTA
    with open(args.reference) as f:
        ref_name = next(SeqIO.parse(f, 'fasta')).id

    # Write samples file (reference, then new root)
    with open(args.output, 'w') as f:
        f.write(f"{ref_name}\n")
        f.write(f"{args.new_root}\n")

    print(f"Created samples file with reference '{ref_name}' and new root '{args.new_root}'")

if __name__ == "__main__":
    main()
