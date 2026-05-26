"""
Emit a 2-column (isolate_id, subtype) CSV from the augmented metadata, for
PastML's --data argument. The raw GISAID subtype column looks like "A / H5N1";
this script normalizes it to "H5N1" so PastML's state labels are compact.
"""

import argparse
import re

import pandas as pd


_SUBTYPE_RE = re.compile(r"(H\d+N\d+)")


def normalize_subtype(subtype_str):
    """Normalize a raw GISAID subtype string (e.g. 'A / H5N1') to 'H5N1'.

    Returns 'unknown' for missing/non-matching values.
    """
    if pd.isna(subtype_str) or not isinstance(subtype_str, str):
        return "unknown"
    match = _SUBTYPE_RE.search(subtype_str.replace(" ", ""))
    if match:
        return match.group(1)
    return "unknown"


def main():
    parser = argparse.ArgumentParser(
        description="Extract isolate_id and normalized H*N* subtype from augmented metadata for PastML."
    )
    parser.add_argument("input", help="Input CSV (combined_metadata_augmented.csv)")
    parser.add_argument("output", help="Output 2-column CSV (isolate_id, subtype)")
    args = parser.parse_args()

    df = pd.read_csv(args.input, usecols=["isolate_id", "subtype"])
    df["subtype"] = df["subtype"].apply(normalize_subtype)
    df.to_csv(args.output, index=False)
    print(f"Wrote {len(df)} rows to {args.output}")
    print("Subtype distribution (top 10):")
    print(df["subtype"].value_counts().head(10).to_string())


if __name__ == "__main__":
    main()
