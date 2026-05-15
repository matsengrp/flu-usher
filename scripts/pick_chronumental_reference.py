"""
Pick a reference sample for chronumental's --reference_node flag.

Strategy: pick a sample whose date sits at a configurable fraction of the
way through the per-segment tree's date range (default 1/3), after
restricting to dates that survive a cluster-density check.

Empirically (PR #37 H5 reference-choice experiment) anchoring on a date
in the first half of the valid range gives notably better leaf residuals
than the earliest cluster-passing date — H5 median |residual| 4 d at
1994 vs 13 d at 1968. Strict midpoint (1/2) works for most subtypes but
overflows pandas' 292-year Timedelta limit on H1 because that tree's
date range starts late (1976) and an anchor at 2001 pushes some deep
internal node 346 years from origin. The 1/3-fraction default biases
the anchor slightly earlier — early enough to keep all subtypes within
pandas' window, late enough to retain the residual improvements seen
with midpoint anchors.

Algorithm:

1. From the curated MSA, drop the curated root and any tips pruned by the
   long-branch filter (passed via --dropped-tips).
2. Keep only samples whose collection_date is a parseable YYYY-MM-DD.
3. Find all distinct dates with at least --min-cluster-size dated samples
   within +/- --cluster-window-years. Call these the "valid" dates.
4. Target date = earliest_valid + (latest_valid - earliest_valid) *
   --target-fraction.
5. Choose the valid date closest to the target date.
6. Among samples at the chosen date, tie-break by smallest isolate_id for
   determinism across runs.
"""
import argparse
import lzma
import sys
from datetime import datetime, timedelta
from Bio import SeqIO
import pandas as pd


def parse_iso_date(value):
    if pd.isna(value) or not isinstance(value, str):
        return None
    try:
        return datetime.strptime(value.strip(), "%Y-%m-%d")
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Pick a chronumental reference sample near the date-range midpoint"
    )
    parser.add_argument("--curated-msa", required=True,
                        help="Path to curated_msa.fasta.xz (source of truth for samples in tree)")
    parser.add_argument("--metadata", required=True,
                        help="Path to augmented metadata CSV")
    parser.add_argument("--root", required=True,
                        help="Path to curated_root.fasta (excluded from candidate set)")
    parser.add_argument("--dropped-tips", required=False, default=None,
                        help="Optional TSV of tips pruned by filter_long_branches "
                             "(first column = sample name, header in row 1). Any "
                             "isolate listed here is removed from the candidate set "
                             "so the chosen reference is guaranteed to exist in the "
                             "newick chronumental will actually receive.")
    parser.add_argument("--output", required=True,
                        help="Output file holding the chosen reference isolate_id")
    parser.add_argument("--cluster-window-years", type=float, default=10.0,
                        help="Half-width (in years) of the window used to count nearby "
                             "samples when judging whether a candidate date is part of "
                             "a real cluster. A date is 'valid' if it has at least "
                             "--min-cluster-size dated samples within +/- this many years.")
    parser.add_argument("--min-cluster-size", type=int, default=3,
                        help="Minimum number of dated samples (including those at the "
                             "candidate date itself) required within the cluster window. "
                             "Candidate dates with fewer nearby samples are treated as "
                             "isolated metadata-error outliers and excluded from the "
                             "valid set.")
    parser.add_argument("--target-fraction", type=float, default=1.0 / 3.0,
                        help="Fraction along the valid date range at which to target "
                             "the chosen reference: target = earliest + fraction * "
                             "(latest - earliest). 0.0 = earliest, 0.5 = midpoint, "
                             "1.0 = latest. Default 1/3 biases the anchor toward the "
                             "earlier portion of the range to avoid pandas Timedelta "
                             "overflow on deep internal-node date predictions while "
                             "retaining the residual improvements of mid-range anchors.")
    args = parser.parse_args()

    print(f"Reading root sequence from {args.root}")
    with open(args.root) as f:
        root_name = next(SeqIO.parse(f, "fasta")).id
    print(f"Root sequence: {root_name}")

    print(f"Reading curated MSA from {args.curated_msa}")
    msa_sample_ids = set()
    with lzma.open(args.curated_msa, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            msa_sample_ids.add(record.id)
    print(f"Found {len(msa_sample_ids)} samples in curated MSA")

    print(f"Reading metadata from {args.metadata}")
    df = pd.read_csv(args.metadata)
    print(f"Metadata contains {len(df)} rows")

    for col in ("isolate_id", "collection_date"):
        if col not in df.columns:
            print(f"ERROR: required column '{col}' missing from metadata", file=sys.stderr)
            sys.exit(1)

    dropped_ids = set()
    if args.dropped_tips:
        with open(args.dropped_tips) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue  # header
                name = line.split("\t", 1)[0].strip()
                if name:
                    dropped_ids.add(name)
        print(f"Excluding {len(dropped_ids)} samples pruned by long-branch filter")

    eligible_ids = msa_sample_ids - {root_name} - dropped_ids
    candidates = df[df["isolate_id"].isin(eligible_ids)].copy()
    candidates["parsed_date"] = candidates["collection_date"].apply(parse_iso_date)
    candidates = candidates.dropna(subset=["parsed_date"])
    print(f"{len(candidates)} non-root, non-pruned MSA samples have a parseable "
          f"YYYY-MM-DD date")

    if len(candidates) == 0:
        print("ERROR: no non-root samples have a parseable date; cannot pick a reference",
              file=sys.stderr)
        sys.exit(1)

    # Find all dates with sufficient cluster density.
    sample_dates = candidates["parsed_date"].values
    unique_dates = sorted(candidates["parsed_date"].unique())
    window = timedelta(days=args.cluster_window_years * 365.25)
    valid_dates = []
    for d in unique_dates:
        nearby = ((sample_dates >= d - window) & (sample_dates <= d + window)).sum()
        if nearby >= args.min_cluster_size:
            valid_dates.append(d)
    if not valid_dates:
        print(f"ERROR: no date has at least {args.min_cluster_size} samples within "
              f"+/-{args.cluster_window_years} years", file=sys.stderr)
        sys.exit(1)

    earliest_valid = valid_dates[0]
    latest_valid = valid_dates[-1]
    span = latest_valid - earliest_valid
    target = earliest_valid + span * args.target_fraction
    print(f"Valid date range: {pd.Timestamp(earliest_valid).date()} to "
          f"{pd.Timestamp(latest_valid).date()}; target "
          f"({args.target_fraction:.3f} of range) = "
          f"{pd.Timestamp(target).date()}")

    # Pick the valid date closest to the target.
    chosen_date = min(valid_dates, key=lambda d: abs((d - target).total_seconds()))
    gap_to_target = abs((chosen_date - target).days)
    print(f"Closest valid date to target: {pd.Timestamp(chosen_date).date()} "
          f"({gap_to_target} days from target)")

    at_chosen = candidates[candidates["parsed_date"] == chosen_date]
    at_chosen = at_chosen.sort_values("isolate_id")
    chosen = at_chosen.iloc[0]
    ref_id = chosen["isolate_id"]
    n_tied = len(at_chosen)
    print(f"Chosen reference sample: {ref_id} ({pd.Timestamp(chosen_date).date()})"
          + (f"; broke {n_tied}-way tie by smallest isolate_id" if n_tied > 1 else ""))

    with open(args.output, "w") as f:
        f.write(f"{ref_id}\n")
    print(f"Wrote reference sample to {args.output}")


if __name__ == "__main__":
    main()
