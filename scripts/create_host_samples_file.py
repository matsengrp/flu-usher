"""
Create a samples file for matUtils extract filtered by host group.
"""
import argparse
import sys
from Bio import SeqIO
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Create samples file filtered by host group for matUtils extract"
    )
    parser.add_argument("--curated-msa", required=True,
                       help="Path to curated_msa.fasta.xz (source of truth for samples in tree)")
    parser.add_argument("--metadata", required=True,
                       help="Path to combined_metadata_with_host_groups.csv")
    parser.add_argument("--host-group", required=True,
                       help="Host group to filter (e.g., 'human', 'avian')")
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
    print(f"Reading curated MSA from {args.curated_msa}")
    msa_sample_ids = set()
    for record in SeqIO.parse(args.curated_msa, 'fasta'):
        msa_sample_ids.add(record.id)
    print(f"Found {len(msa_sample_ids)} samples in curated MSA")

    # Read metadata
    print(f"Reading metadata from {args.metadata}")
    df = pd.read_csv(args.metadata)
    print(f"Metadata contains {len(df)} rows")

    # Verify all MSA samples are present in metadata
    metadata_ids = set(df['isolate_id'])
    missing_in_metadata = msa_sample_ids - metadata_ids
    if missing_in_metadata:
        print(f"ERROR: {len(missing_in_metadata)} samples in MSA are missing from metadata",
              file=sys.stderr)
        print(f"First few missing samples: {list(missing_in_metadata)[:5]}", file=sys.stderr)
        sys.exit(1)

    # Filter metadata by host group
    host_df = df[df['host_group'] == args.host_group]
    print(f"Found {len(host_df)} samples with host_group='{args.host_group}' in metadata")

    # Match filtered metadata against MSA sample IDs
    host_sample_ids = set(host_df['isolate_id']) & msa_sample_ids
    print(f"Of those, {len(host_sample_ids)} are present in curated MSA")

    # Warn if no samples found beyond root
    if len(host_sample_ids) == 0:
        print(f"WARNING: No samples found for host_group='{args.host_group}' in curated MSA",
              file=sys.stderr)
        print("Output will contain only the root sequence", file=sys.stderr)

    # Sort for reproducibility (excluding root which goes first)
    sorted_samples = sorted(host_sample_ids)

    # Write samples file: root first, then sorted host-group samples
    with open(args.output, 'w') as f:
        f.write(f"{root_name}\n")
        for sample_id in sorted_samples:
            f.write(f"{sample_id}\n")

    total_samples = 1 + len(sorted_samples)
    print(f"Wrote {total_samples} samples to {args.output} (1 root + {len(sorted_samples)} host-specific)")

if __name__ == "__main__":
    main()
