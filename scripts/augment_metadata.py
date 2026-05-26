"""
Augment combined metadata with host_group, geographic_group, and temporal_group columns.

Reads the combined metadata CSV once, adds all grouping columns, and writes a single output file.
"""

import argparse
from datetime import datetime

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


def get_temporal_group(collection_date, median_date):
    """Assign temporal group based on median date split.

    Returns 'early', 'late', or 'unknown' (for missing/invalid dates).
    """
    if pd.isna(collection_date) or not isinstance(collection_date, str):
        return "unknown"
    try:
        dt = datetime.strptime(collection_date.strip(), "%Y-%m-%d")
    except ValueError:
        return "unknown"

    if dt <= median_date:
        return "early"
    return "late"


def compute_global_median_date(dates_series):
    """Compute the median date from a series of date strings (YYYY-MM-DD)."""
    valid_dates = []
    for d in dates_series:
        if pd.isna(d) or not isinstance(d, str):
            continue
        try:
            valid_dates.append(datetime.strptime(d.strip(), "%Y-%m-%d"))
        except ValueError:
            continue

    valid_dates.sort()
    if not valid_dates:
        raise ValueError("No valid dates found in collection_date column")

    median_date = valid_dates[len(valid_dates) // 2]
    return median_date


def main():
    parser = argparse.ArgumentParser(
        description="Augment metadata with host_group, geographic_group, and temporal_group columns"
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

    # Temporal group
    median_date = compute_global_median_date(df["collection_date"])
    print(f"\nGlobal median collection date: {median_date.date()}")
    df["temporal_group"] = df["collection_date"].apply(
        lambda d: get_temporal_group(d, median_date)
    )
    print(f"\nTemporal group distribution:")
    print(df["temporal_group"].value_counts().to_string())

    df.to_csv(args.output, index=False)
    print(f"\nWrote augmented metadata to {args.output}")


if __name__ == "__main__":
    main()
