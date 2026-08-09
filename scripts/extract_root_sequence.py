"""
Extract and infer root sequence from tree mutation paths with validation against MSA.
"""
import argparse
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import iter_path_mutations, open_sequence_file, setup_logging

logger = setup_logging(__name__)


def extract_sequence_from_msa(msa_file, sequence_name):
    """Extract a specific sequence from a compressed MSA file by sequence name."""
    with open_sequence_file(msa_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id == sequence_name:
                return str(record.seq)
    raise ValueError(f"Sequence '{sequence_name}' not found in MSA")


def parse_mutation_path(paths_file, new_root_name):
    """Return the root-to-tip mutations for one sample, as parsed triples.

    Tokenizing is delegated to utils.iter_path_mutations rather than repeated
    here. The hand-rolled version this replaces split on ':' by index, which
    raised a bare IndexError on a colon-less chunk instead of saying what was
    wrong, and was a fourth copy of the format the shared parser exists to own.
    """
    with open(paths_file) as f:
        for line in f:
            # rstrip('\n'), not strip(): a sample whose path is empty writes
            # "name\t\n", and strip() eats the tab, so the row split to one
            # field and the sample was reported missing rather than
            # mutation-free. Every root_paths.txt on disk has such a row -- the
            # reference against itself -- so this only escaped notice because
            # nothing ever looks that row up.
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0] == new_root_name:
                return list(iter_path_mutations(parts[1]))
    raise ValueError(f"Sample {new_root_name} not found in paths file")


def apply_mutations(reference_seq, mutations):
    """Apply a root-to-tip mutation path to the reference, 1-based positions.

    A parent-base mismatch means the path and the reference are not the pair we
    think they are, so the sequence this builds would be wrong. It used to warn
    and carry on; validate_sequences() would then fail further downstream on
    the consequence rather than the cause. Fail here instead, matching
    rebase_mat_root.root_sequence(), which does the same job from the protobuf.
    """
    seq = list(str(reference_seq))
    for parent, position, mutant in mutations:
        token = f"{parent}{position}{mutant}"
        index = position - 1
        if not 0 <= index < len(seq):
            logger.error(
                f"mutation {token} is outside the {len(seq)}-base reference"
            )
            sys.exit(1)
        if seq[index] != parent:
            logger.error(
                f"mutation {token} expects {parent} at position {position}, "
                f"but the reference has {seq[index]}"
            )
            sys.exit(1)
        seq[index] = mutant
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

    logger.info(f"Validation passed: Inferred sequence matches MSA except at {gap_fills} gap positions")


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
    logger.info(f"Loaded reference sequence: {ref_record.id} ({len(ref_record.seq)} bp)")

    # Parse mutation path for new root. Both failure modes -- the sample not
    # being in the file, and a token the shared parser rejects -- are the
    # pipeline's problem to report, not a traceback's: the other two failures
    # in this script already exit(1) with an explanation, so match them.
    try:
        mutations = parse_mutation_path(args.paths, args.new_root_name)
    except ValueError as error:
        logger.error(f"could not read the path to {args.new_root_name}: {error}")
        sys.exit(1)
    logger.info(f"Found {len(mutations)} mutations along path to {args.new_root_name}")

    # Apply mutations to reference sequence
    inferred_seq = apply_mutations(ref_record.seq, mutations)
    logger.info(f"Applied mutations to generate inferred root sequence ({len(inferred_seq)} bp)")

    # Extract MSA sequence for validation
    msa_seq = extract_sequence_from_msa(args.msa, args.new_root_name)
    logger.info(f"Extracted MSA sequence: {args.new_root_name} ({len(msa_seq)} bp)")

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
    logger.info(f"Wrote inferred root sequence to {args.output}")


if __name__ == "__main__":
    main()