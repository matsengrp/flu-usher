"""
Pick the earliest non-root sample with a parseable YYYY-MM-DD collection_date
that is present in the curated MSA, and write its isolate_id to a single-line
file. Used as `chronumental --reference_node` for the full per-subtype tree so
chronumental anchors on a real, well-dated observation rather than the curated
root sequence (which may be much older or undated).

Skips isolated very-early dates: a candidate date is only accepted if at least
one other dated sample falls within `--max-gap-years` of it (default 20). This
filters out GISAID metadata errors like "1905" labels on 21st-century isolates
that would otherwise distort the time scale and break chronumental's date
arithmetic (pandas timedeltas overflow ~292 years).

Mirrors the reference-picking logic in scripts/create_cohort_samples_file.py
(lines 137-151), without the cohort host / min-date filters.
"""
import argparse
import lzma
import sys
from datetime import datetime
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
        description="Pick the earliest non-root dated sample present in the curated MSA"
    )
    parser.add_argument("--curated-msa", required=True,
                        help="Path to curated_msa.fasta.xz (source of truth for samples in tree)")
    parser.add_argument("--metadata", required=True,
                        help="Path to augmented metadata CSV")
    parser.add_argument("--root", required=True,
                        help="Path to curated_root.fasta (excluded from candidate set)")
    parser.add_argument("--output", required=True,
                        help="Output file holding the isolate_id of the earliest dated sample")
    parser.add_argument("--max-gap-years", type=float, default=20.0,
                        help="Maximum allowed gap (in years) between the chosen date and "
                             "the next-later distinct date. Earlier dates that are isolated "
                             "by a larger gap are treated as metadata errors and skipped.")
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

    candidates = df[df["isolate_id"].isin(msa_sample_ids - {root_name})].copy()
    candidates["parsed_date"] = candidates["collection_date"].apply(parse_iso_date)
    candidates = candidates.dropna(subset=["parsed_date"])
    print(f"{len(candidates)} non-root MSA samples have a parseable YYYY-MM-DD date")

    if len(candidates) == 0:
        print("ERROR: no non-root samples have a parseable date; cannot pick a reference",
              file=sys.stderr)
        sys.exit(1)

    # Walk distinct dates from oldest to newest; skip isolated outliers whose
    # gap to the next-later distinct date exceeds --max-gap-years.
    unique_dates = sorted(candidates["parsed_date"].unique())
    max_gap_days = args.max_gap_years * 365.25
    chosen_date = None
    for i, d in enumerate(unique_dates):
        if i + 1 < len(unique_dates):
            gap_days = (unique_dates[i + 1] - d).days
        else:
            gap_days = 0  # last date: no successor to gap to; never isolated
        if gap_days <= max_gap_days:
            chosen_date = d
            break
        else:
            n_at_d = (candidates["parsed_date"] == d).sum()
            print(f"  skipping isolated date {pd.Timestamp(d).date()} "
                  f"({n_at_d} sample(s)); next date is "
                  f"{pd.Timestamp(unique_dates[i + 1]).date()} ({gap_days / 365.25:.1f}-year gap)")
    if chosen_date is None:
        print("ERROR: no candidate date has a non-isolated neighbor within "
              f"{args.max_gap_years} years", file=sys.stderr)
        sys.exit(1)

    earliest = candidates[candidates["parsed_date"] == chosen_date]
    earliest = earliest.sort_values("isolate_id")
    chosen = earliest.iloc[0]
    ref_id = chosen["isolate_id"]
    ref_date = chosen["parsed_date"].date()
    n_tied = len(earliest)
    print(f"Chosen reference sample: {ref_id} ({ref_date})"
          + (f"; broke {n_tied}-way tie by smallest isolate_id" if n_tied > 1 else ""))

    with open(args.output, "w") as f:
        f.write(f"{ref_id}\n")
    print(f"Wrote reference sample to {args.output}")


if __name__ == "__main__":
    main()
