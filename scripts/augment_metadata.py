"""
Augment combined metadata with host_group and geographic_group columns.

Reads the combined metadata CSV once, adds grouping columns, and writes a single output file.
"""

import argparse

import pandas as pd

from simplified_host_classifier import get_simplified_host_group


def get_geographic_group(location):
    """Extract continent from GISAID location and map to geographic group.

    GISAID location format: "Continent / Country / Region / City"
    Maps to: north_america, europe, asia, other
    """
    if pd.isna(location) or not isinstance(location, str) or not location.strip():
        return "other"

    continent = location.split("/")[0].strip()

    mapping = {
        "North America": "north_america",
        "Europe": "europe",
        "Asia": "asia",
    }
    return mapping.get(continent, "other")


def main():
    parser = argparse.ArgumentParser(
        description="Augment metadata with host_group and geographic_group columns"
    )
    parser.add_argument("input", help="Input CSV file (combined_metadata.csv)")
    parser.add_argument("output", help="Output CSV file with new columns added")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    print(f"Read {len(df)} records from {args.input}")

    # Host group
    df["host_group"] = df["host"].apply(get_simplified_host_group)
    print(f"\nHost group distribution:")
    print(df["host_group"].value_counts().to_string())

    # Geographic group
    df["geographic_group"] = df["location"].apply(get_geographic_group)
    print(f"\nGeographic group distribution:")
    print(df["geographic_group"].value_counts().to_string())

    df.to_csv(args.output, index=False)
    print(f"\nWrote augmented metadata to {args.output}")


if __name__ == "__main__":
    main()
