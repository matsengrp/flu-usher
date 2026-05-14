"""
Pick the earliest non-root sample with a parseable YYYY-MM-DD collection_date
that is present in the curated MSA, and write its isolate_id to a single-line
file. Used as `chronumental --reference_node` for the per-segment tree so
chronumental anchors on a real, well-dated observation rather than the curated
root sequence (which may be much older or undated).

Cluster-density filter: a candidate date is accepted only if at least
`--min-cluster-size` dated samples (including those at the candidate date
itself) fall within +/- `--cluster-window-years` of it. The gap-to-next-date
filter that this replaces was too lenient for multi-subtype internal-segment
trees, where GISAID "1905" metadata errors slipped through because legitimate
1918 samples sat within 13 years. Requiring a minimum cluster size around the
candidate forces the chosen anchor to be a "real" early sample, not a
date-typo outlier.
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
        description="Pick the earliest non-root dated sample present in the curated MSA"
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
                        help="Output file holding the isolate_id of the earliest dated sample")
    parser.add_argument("--cluster-window-years", type=float, default=10.0,
                        help="Half-width (in years) of the window used to count nearby "
                             "samples when judging whether a candidate date is part of "
                             "a real cluster. The chosen date must have at least "
                             "--min-cluster-size dated samples within +/- this many years.")
    parser.add_argument("--min-cluster-size", type=int, default=3,
                        help="Minimum number of dated samples (including those at the "
                             "candidate date itself) required within the cluster window. "
                             "Candidate dates with fewer nearby samples are treated as "
                             "isolated metadata-error outliers and skipped.")
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

    # Walk distinct dates from oldest to newest; skip any whose cluster-window
    # neighbor count is below --min-cluster-size.
    sample_dates = candidates["parsed_date"].values
    unique_dates = sorted(candidates["parsed_date"].unique())
    window = timedelta(days=args.cluster_window_years * 365.25)
    chosen_date = None
    for d in unique_dates:
        nearby = ((sample_dates >= d - window) & (sample_dates <= d + window)).sum()
        if nearby >= args.min_cluster_size:
            chosen_date = d
            break
        else:
            n_at_d = int((candidates["parsed_date"] == d).sum())
            print(f"  skipping isolated date {pd.Timestamp(d).date()} "
                  f"({n_at_d} sample(s) at this date; {nearby} sample(s) within "
                  f"+/-{args.cluster_window_years:g}y, need >= {args.min_cluster_size})")
    if chosen_date is None:
        print(f"ERROR: no candidate date has at least {args.min_cluster_size} samples "
              f"within +/-{args.cluster_window_years} years", file=sys.stderr)
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
