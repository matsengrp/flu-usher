"""
Create a samples file for matUtils extract filtered by a cohort definition
(subtype + host_group + min collection_date) applied simultaneously.

Output format matches the other samples-file scripts: the curated root
sequence on line 1, then matched sample IDs (sorted) on subsequent lines.
"""
import argparse
import lzma
import re
import sys
from datetime import datetime
from Bio import SeqIO
import pandas as pd


def normalize_subtype(value):
    if not isinstance(value, str):
        return ""
    return re.sub(r"^A\s*/\s*", "", value).strip().upper()


def parse_iso_date(value):
    if pd.isna(value) or not isinstance(value, str):
        return None
    try:
        return datetime.strptime(value.strip(), "%Y-%m-%d")
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Create samples file for cohort subtree extraction"
    )
    parser.add_argument("--curated-msa", required=True,
                        help="Path to curated_msa.fasta.xz (source of truth for samples in tree)")
    parser.add_argument("--metadata", required=True,
                        help="Path to augmented metadata CSV")
    parser.add_argument("--root", required=True,
                        help="Path to curated_root.fasta (always included as line 1)")
    parser.add_argument("--subtype", required=True,
                        help="GISAID subtype to keep (e.g. 'H1N1'); compared after stripping 'A /'")
    parser.add_argument("--host", required=True,
                        help="host_group value to keep (e.g. 'human')")
    parser.add_argument("--min-date", required=True,
                        help="ISO date 'YYYY-MM-DD'; samples with collection_date strictly after this are kept")
    parser.add_argument("--output", required=True, help="Output samples.txt file path")
    parser.add_argument("--reference-output", required=True,
                        help="Output file holding the isolate_id of the earliest non-root "
                             "cohort sample (used as chronumental --reference_node). "
                             "Empty if no cohort samples were found.")
    args = parser.parse_args()

    min_date = parse_iso_date(args.min_date)
    if min_date is None:
        print(f"ERROR: --min-date '{args.min_date}' is not a valid YYYY-MM-DD date",
              file=sys.stderr)
        sys.exit(1)
    target_subtype = normalize_subtype(args.subtype)
    if not target_subtype:
        print(f"ERROR: --subtype '{args.subtype}' normalized to empty string", file=sys.stderr)
        sys.exit(1)

    print(f"Reading root sequence from {args.root}")
    with open(args.root) as f:
        root_name = next(SeqIO.parse(f, "fasta")).id
    print(f"Root sequence: {root_name}")

    print(f"Reading curated MSA from {args.curated_msa}")
    msa_sample_ids = set()
    reference_id = None
    with lzma.open(args.curated_msa, "rt") as f:
        for i, record in enumerate(SeqIO.parse(f, "fasta")):
            if i == 0:
                reference_id = record.id
                print(f"Reference sequence (first in MSA): {reference_id}")
            msa_sample_ids.add(record.id)
    print(f"Found {len(msa_sample_ids)} samples in curated MSA")

    print(f"Reading metadata from {args.metadata}")
    df = pd.read_csv(args.metadata)
    print(f"Metadata contains {len(df)} rows")

    for col in ("isolate_id", "subtype", "host_group", "collection_date"):
        if col not in df.columns:
            print(f"ERROR: required column '{col}' missing from metadata", file=sys.stderr)
            sys.exit(1)

    metadata_ids = set(df["isolate_id"])
    missing_in_metadata = msa_sample_ids - metadata_ids - {reference_id}
    if missing_in_metadata:
        print(f"ERROR: {len(missing_in_metadata)} samples in MSA are missing from metadata",
              file=sys.stderr)
        print(f"First few missing: {list(missing_in_metadata)[:5]}", file=sys.stderr)
        sys.exit(1)

    msa_df = df[df["isolate_id"].isin(msa_sample_ids - {reference_id})].copy()
    print(f"Found {len(msa_df)} MSA samples in metadata")

    msa_df["subtype_norm"] = msa_df["subtype"].apply(normalize_subtype)
    msa_df["parsed_date"] = msa_df["collection_date"].apply(parse_iso_date)

    after_subtype = msa_df[msa_df["subtype_norm"] == target_subtype]
    print(f"Samples matching subtype='{target_subtype}': {len(after_subtype)}")
    after_host = after_subtype[after_subtype["host_group"] == args.host]
    print(f"Samples also matching host='{args.host}': {len(after_host)}")
    after_date = after_host.dropna(subset=["parsed_date"])
    after_date = after_date[after_date["parsed_date"] > min_date]
    print(f"Samples also with collection_date > {min_date.date()}: {len(after_date)}")

    cohort_sample_ids = set(after_date["isolate_id"]) & msa_sample_ids
    print(f"Of those, {len(cohort_sample_ids)} are present in curated MSA")

    if len(cohort_sample_ids) == 0:
        print("WARNING: no cohort samples found; writing root-only file (cohort will be skipped)",
              file=sys.stderr)

    sorted_samples = sorted(cohort_sample_ids - {root_name})

    with open(args.output, "w") as f:
        f.write(f"{root_name}\n")
        for sample_id in sorted_samples:
            f.write(f"{sample_id}\n")

    total = 1 + len(sorted_samples)
    print(f"Wrote {total} samples to {args.output} (1 root + {len(sorted_samples)} cohort)")

    # Pick the earliest cohort sample (excluding root) as a chronumental reference.
    # Cohort filtering already required a valid parsed_date, so earliest is well-defined
    # whenever the cohort is non-empty.
    cohort_df = after_date[after_date["isolate_id"].isin(cohort_sample_ids - {root_name})]
    if len(cohort_df) == 0:
        with open(args.reference_output, "w") as f:
            pass
        print(f"Wrote empty reference file to {args.reference_output} (no cohort samples)")
    else:
        earliest = cohort_df.loc[cohort_df["parsed_date"].idxmin()]
        ref_id = earliest["isolate_id"]
        ref_date = earliest["parsed_date"].date()
        with open(args.reference_output, "w") as f:
            f.write(f"{ref_id}\n")
        print(f"Earliest non-root cohort sample: {ref_id} ({ref_date}); wrote to {args.reference_output}")


if __name__ == "__main__":
    main()
