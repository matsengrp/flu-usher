"""
Extract and infer root sequence from tree mutation paths with validation against MSA.
"""
import argparse
import lzma
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_sequence_from_msa(msa_file, sequence_name):
    """Extract a specific sequence from a compressed MSA file by sequence name."""
    with lzma.open(msa_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id == sequence_name:
                return str(record.seq)
    raise ValueError(f"Sequence '{sequence_name}' not found in MSA")


def parse_mutation_path(paths_file, new_root_name):
    """Extract all mutations from the path to the new root."""
    with open(paths_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2 and parts[0] == new_root_name:
                # Parse mutation path
                mutations = []
                if parts[1]:  # Non-empty path
                    for node in parts[1].split():
                        node_muts = node.split(':')[1]
                        if node_muts:
                            mutations.extend(node_muts.split(','))
                return mutations
    raise ValueError(f"Sample {new_root_name} not found in paths file")


def apply_mutations(reference_seq, mutations):
    """Apply mutations to reference sequence."""
    seq = list(str(reference_seq))
    for mut in mutations:
        orig_base = mut[0]
        new_base = mut[-1]
        pos = int(mut[1:-1]) - 1  # Convert to 0-based

        # Validate original base matches
        if seq[pos] != orig_base:
            print(f"Warning: Expected {orig_base} at position {pos+1}, found {seq[pos]}")

        seq[pos] = new_base
    return ''.join(seq)


def validate_sequences(msa_seq, inferred_seq):
    """Validate inferred sequence against MSA sequence."""
    if len(msa_seq) != len(inferred_seq):
        raise ValueError(f"Length mismatch: MSA={len(msa_seq)}, inferred={len(inferred_seq)}")

    mismatches = []
    gap_fills = 0

    for i, (msa_base, inferred_base) in enumerate(zip(msa_seq, inferred_seq)):
        if msa_base == '-':
            # Gap in MSA should be a nucleotide in inferred
            if inferred_base not in 'ACGT':
                raise ValueError(f"Position {i+1}: MSA has gap but inferred is not a nucleotide: {inferred_base}")
            gap_fills += 1
        else:
            # Non-gap in MSA should match inferred exactly
            if msa_base != inferred_base:
                mismatches.append(f"Position {i+1}: MSA has {msa_base} but inferred has {inferred_base}")

    if mismatches:
        error_msg = "\n".join(mismatches)
        raise ValueError(f"Validation failed with {len(mismatches)} mismatches:\n{error_msg}")

    print(f"✓ Validation passed: Inferred sequence matches MSA except at {gap_fills} gap positions")


def main():
    parser = argparse.ArgumentParser(
        description="Extract and infer root sequence from tree mutation paths with validation against MSA"
    )
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--msa", required=True, help="MSA FASTA file (.fasta.xz)")
    parser.add_argument("--paths", required=True, help="Mutation paths file from matUtils extract")
    parser.add_argument("--new-root-name", required=True, help="Name of the new root sample")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    # Read reference sequence
    with open(args.reference) as f:
        ref_record = next(SeqIO.parse(f, 'fasta'))
    print(f"Loaded reference sequence: {ref_record.id} ({len(ref_record.seq)} bp)")

    # Parse mutation path for new root
    mutations = parse_mutation_path(args.paths, args.new_root_name)
    print(f"Found {len(mutations)} mutations along path to {args.new_root_name}")

    # Apply mutations to reference sequence
    inferred_seq = apply_mutations(ref_record.seq, mutations)
    print(f"Applied mutations to generate inferred root sequence ({len(inferred_seq)} bp)")

    # Extract MSA sequence for validation
    msa_seq = extract_sequence_from_msa(args.msa, args.new_root_name)
    print(f"Extracted MSA sequence: {args.new_root_name} ({len(msa_seq)} bp)")

    # Validate inferred sequence against MSA
    validate_sequences(msa_seq, inferred_seq)

    # Write output
    output_record = SeqRecord(
        Seq(inferred_seq),
        id=args.new_root_name,
        description=""
    )
    with open(args.output, 'w') as out:
        SeqIO.write(output_record, out, 'fasta')
    print(f"Wrote inferred root sequence to {args.output}")


if __name__ == "__main__":
    main()