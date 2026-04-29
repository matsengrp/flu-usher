"""Build the global strain<TAB>date TSV consumed by chronumental.

Chronumental looks up a date for each leaf of its input tree from the
``--dates`` file and silently ignores rows whose strain is not in the tree.
That lets us emit a single global file shared by every (segment, subtype)
chronumental job rather than filtering per tree.

Only rows with a fully-resolved YYYY-MM-DD ``collection_date`` are kept;
year-only, year-month-only, and missing dates are dropped because
chronumental cannot use them.
"""

import argparse
from datetime import datetime

import pandas as pd


def parse_full_date(value):
    if pd.isna(value) or not isinstance(value, str):
        return None
    try:
        return datetime.strptime(value.strip(), "%Y-%m-%d").date()
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True,
                        help="Augmented metadata CSV (combined_metadata_augmented.csv)")
    parser.add_argument("--output", required=True,
                        help="Output TSV with header 'strain<TAB>date'")
    args = parser.parse_args()

    df = pd.read_csv(args.metadata, usecols=["isolate_id", "collection_date"])
    print(f"Read {len(df)} rows from {args.metadata}")

    df["parsed_date"] = df["collection_date"].apply(parse_full_date)
    kept = df.dropna(subset=["parsed_date"])
    dropped = len(df) - len(kept)
    print(f"Kept {len(kept)} rows with full YYYY-MM-DD dates "
          f"(dropped {dropped} with partial or missing dates)")

    out = pd.DataFrame({
        "strain": kept["isolate_id"],
        "date": kept["parsed_date"].astype(str),
    })
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out)} strain/date pairs to {args.output}")


if __name__ == "__main__":
    main()
