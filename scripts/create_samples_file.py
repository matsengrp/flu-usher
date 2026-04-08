"""
Create a samples file for matUtils extract filtered by any metadata column.

Generalized version of create_host_samples_file.py that accepts an arbitrary
column name and value for filtering.
"""
import argparse
import lzma
import sys
from Bio import SeqIO
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Create samples file filtered by metadata column for matUtils extract"
    )
    parser.add_argument("--curated-msa", required=True,
                       help="Path to curated_msa.fasta.xz (source of truth for samples in tree)")
    parser.add_argument("--metadata", required=True,
                       help="Path to augmented metadata CSV")
    parser.add_argument("--column", required=True,
                       help="Metadata column to filter on (e.g., 'host_group', 'geographic_group')")
    parser.add_argument("--value", required=True,
                       help="Value to match in the column (e.g., 'human', 'north_america')")
    parser.add_argument("--root", required=True,
                       help="Path to curated_root.fasta (final root sequence, always included)")
    parser.add_argument("--output", required=True,
                       help="Output samples.txt file path")
    args = parser.parse_args()

    # Extract root name from curated_root.fasta
    print(f"Reading root sequence from {args.root}")
    with open(args.root) as f:
        root_name = next(SeqIO.parse(f, 'fasta')).id
    print(f"Root sequence: {root_name}")

    # Extract all sample IDs from curated MSA (these are samples in the tree)
    # The first sequence is the reference and won't be in metadata
    print(f"Reading curated MSA from {args.curated_msa}")
    msa_sample_ids = set()
    reference_id = None
    with lzma.open(args.curated_msa, 'rt') as f:
        for i, record in enumerate(SeqIO.parse(f, 'fasta')):
            if i == 0:
                reference_id = record.id
                print(f"Reference sequence (first in MSA): {reference_id}")
            msa_sample_ids.add(record.id)
    print(f"Found {len(msa_sample_ids)} samples in curated MSA")

    # Read metadata
    print(f"Reading metadata from {args.metadata}")
    df = pd.read_csv(args.metadata)
    print(f"Metadata contains {len(df)} rows")

    # Validate that the column exists
    if args.column not in df.columns:
        print(f"ERROR: Column '{args.column}' not found in metadata. "
              f"Available columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Verify all MSA samples (except reference) are present in metadata
    metadata_ids = set(df['isolate_id'])
    missing_in_metadata = msa_sample_ids - metadata_ids - {reference_id}
    if missing_in_metadata:
        print(f"ERROR: {len(missing_in_metadata)} samples in MSA are missing from metadata",
              file=sys.stderr)
        print(f"First few missing samples: {list(missing_in_metadata)[:5]}", file=sys.stderr)
        sys.exit(1)

    # Filter metadata by column value
    filtered_df = df[df[args.column] == args.value]
    print(f"Found {len(filtered_df)} samples with {args.column}='{args.value}' in metadata")

    # Match filtered metadata against MSA sample IDs
    filtered_sample_ids = set(filtered_df['isolate_id']) & msa_sample_ids
    print(f"Of those, {len(filtered_sample_ids)} are present in curated MSA")

    # Warn if no samples found beyond root
    if len(filtered_sample_ids) == 0:
        print(f"WARNING: No samples found for {args.column}='{args.value}' in curated MSA",
              file=sys.stderr)
        print("Output will contain only the root sequence", file=sys.stderr)

    # Sort for reproducibility (excluding root which goes first)
    # Remove root from filtered_sample_ids to avoid writing it twice
    sorted_samples = sorted(filtered_sample_ids - {root_name})

    # Write samples file: root first, then sorted filtered samples
    with open(args.output, 'w') as f:
        f.write(f"{root_name}\n")
        for sample_id in sorted_samples:
            f.write(f"{sample_id}\n")

    total_samples = 1 + len(sorted_samples)
    print(f"Wrote {total_samples} samples to {args.output} (1 root + {len(sorted_samples)} filtered)")


if __name__ == "__main__":
    main()
