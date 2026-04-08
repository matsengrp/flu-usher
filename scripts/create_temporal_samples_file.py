"""
Create a samples file for matUtils extract filtered by temporal group (early/late).

Unlike the generic create_samples_file.py, this script computes the median date
from the per-segment/subtype sample set to ensure balanced splits within each tree.
"""
import argparse
import lzma
import sys
from datetime import datetime
from Bio import SeqIO
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Create samples file filtered by temporal group for matUtils extract"
    )
    parser.add_argument("--curated-msa", required=True,
                       help="Path to curated_msa.fasta.xz (source of truth for samples in tree)")
    parser.add_argument("--metadata", required=True,
                       help="Path to augmented metadata CSV")
    parser.add_argument("--temporal-group", required=True, choices=["early", "late"],
                       help="Temporal group to extract")
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

    # Extract all sample IDs from curated MSA
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

    # Verify all MSA samples (except reference) are present in metadata
    metadata_ids = set(df['isolate_id'])
    missing_in_metadata = msa_sample_ids - metadata_ids - {reference_id}
    if missing_in_metadata:
        print(f"ERROR: {len(missing_in_metadata)} samples in MSA are missing from metadata",
              file=sys.stderr)
        print(f"First few missing samples: {list(missing_in_metadata)[:5]}", file=sys.stderr)
        sys.exit(1)

    # Filter metadata to only samples in this MSA
    msa_df = df[df['isolate_id'].isin(msa_sample_ids - {reference_id})].copy()
    print(f"Found {len(msa_df)} MSA samples in metadata")

    # Parse collection dates, dropping invalid/missing
    def parse_date(d):
        if pd.isna(d) or not isinstance(d, str):
            return None
        try:
            return datetime.strptime(d.strip(), "%Y-%m-%d")
        except ValueError:
            return None

    msa_df["parsed_date"] = msa_df["collection_date"].apply(parse_date)
    dated_df = msa_df.dropna(subset=["parsed_date"])
    n_excluded = len(msa_df) - len(dated_df)
    print(f"Samples with valid dates: {len(dated_df)} (excluded {n_excluded} with missing/invalid dates)")

    if len(dated_df) == 0:
        print("WARNING: No samples with valid dates found", file=sys.stderr)
        # Write root-only file
        with open(args.output, 'w') as f:
            f.write(f"{root_name}\n")
        print(f"Wrote 1 sample to {args.output} (root only)")
        return

    # Compute per-tree median date
    sorted_dates = sorted(dated_df["parsed_date"])
    median_date = sorted_dates[len(sorted_dates) // 2]
    print(f"Per-tree median date: {median_date.date()}")

    # Split into early/late
    if args.temporal_group == "early":
        group_df = dated_df[dated_df["parsed_date"] <= median_date]
    else:
        group_df = dated_df[dated_df["parsed_date"] > median_date]
    print(f"Samples in '{args.temporal_group}' group: {len(group_df)}")

    # Match against MSA sample IDs
    group_sample_ids = set(group_df['isolate_id']) & msa_sample_ids

    if len(group_sample_ids) == 0:
        print(f"WARNING: No samples in temporal group '{args.temporal_group}'",
              file=sys.stderr)
        print("Output will contain only the root sequence", file=sys.stderr)

    # Sort for reproducibility, remove root to avoid writing it twice
    sorted_samples = sorted(group_sample_ids - {root_name})

    # Write samples file: root first, then sorted temporal samples
    with open(args.output, 'w') as f:
        f.write(f"{root_name}\n")
        for sample_id in sorted_samples:
            f.write(f"{sample_id}\n")

    total_samples = 1 + len(sorted_samples)
    print(f"Wrote {total_samples} samples to {args.output} (1 root + {len(sorted_samples)} temporal)")


if __name__ == "__main__":
    main()
