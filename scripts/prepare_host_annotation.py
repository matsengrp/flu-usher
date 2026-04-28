"""
Emit a 2-column (isolate_id, host_group) CSV from the augmented metadata, for
PastML's --data argument.
"""

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Extract isolate_id, host_group from augmented metadata for PastML."
    )
    parser.add_argument("input", help="Input CSV (combined_metadata_augmented.csv)")
    parser.add_argument("output", help="Output 2-column CSV (isolate_id, host_group)")
    args = parser.parse_args()

    df = pd.read_csv(args.input, usecols=["isolate_id", "host_group"])
    df.to_csv(args.output, index=False)
    print(f"Wrote {len(df)} rows to {args.output}")


if __name__ == "__main__":
    main()
