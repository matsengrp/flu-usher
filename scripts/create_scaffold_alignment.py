"""
Cut an alignment down to a small subset of sequences spread evenly over
collection time, to be built into a scaffold tree.

The tree search places sequences in whatever order the alignment happens to be
in, and on HA/H3 it ends up with a high-level shape that is measurably not the
most parsimonious one available: sequences collected from 2006 on sit ~70
branches deeper on the trunk than 2005's without a matching rise in divergence
(issue #53). Optimising ~1000 sequences is a far easier search than optimising
86,000, so the pipeline now builds a backbone from a time-spread subset first
and places everything else onto it.

"Evenly over time" means an equal quota per collection year, not per sequence:
the eight 1963 HA/H3 sequences say far more about the shape of the trunk than
any eight of 2024's 11,786 do. Years that cannot fill their quota hand the
remainder back to be redistributed over the years that still have sequences to
give -- repeatedly, since a redistribution can itself overfill a small year.

Sequences with no collection date are not candidates: they carry no information
about where along the time axis they belong. They are placed later, in the pass
that adds every non-scaffold sequence to the backbone, so nothing is dropped.

The first record is taken to be the curated reference and is always kept, at the
top: faToVcf reads the reference off the alignment's first record. Everything
after it is emitted in the order it appeared, so a randomised input alignment
yields a randomised scaffold and each randomisation gets a different backbone.
"""

import argparse
import csv
import random
from collections import defaultdict

from Bio import SeqIO

from utils import open_sequence_file, setup_logging

logger = setup_logging(__name__)


def read_collection_years(metadata_file, wanted_ids):
    """
    Map each wanted id to its collection year, skipping ones we cannot date.

    Ids missing from the metadata, or carrying a blank or non-numeric date, are
    left out rather than bucketed together -- they are not placeable on the time
    axis, so they are not candidates for a time-stratified sample.
    """
    years = {}
    with open(metadata_file) as handle:
        for row in csv.DictReader(handle):
            isolate_id = row["isolate_id"]
            if isolate_id not in wanted_ids:
                continue
            date = (row.get("collection_date") or "").strip()
            if len(date) < 4 or not date[:4].isdigit():
                continue
            years[isolate_id] = int(date[:4])
    return years


def allocate_quotas(available, total):
    """
    Split `total` draws over the keys of `available` as evenly as possible.

    `available` maps each year to the number of sequences it actually has.
    Returns year -> number to draw, where no value exceeds that year's
    availability and the values sum to min(total, sum of availability).

    Water-filling: give every year that still has room an equal share of what is
    left, let the years that cannot use their full share return the remainder,
    and repeat. Each round either places a sequence or exhausts a year, so this
    terminates.
    """
    quotas = {year: 0 for year in available}
    remaining = min(total, sum(available.values()))

    while remaining > 0:
        open_years = [y for y in available if quotas[y] < available[y]]
        if not open_years:
            break
        share, extra = divmod(remaining, len(open_years))
        # Years with the least room left take the +1s, so a remainder cannot be
        # handed to a year with no space for it while a roomier year sits idle.
        open_years.sort(key=lambda y: (available[y] - quotas[y], y))
        placed = 0
        for i, year in enumerate(open_years):
            want = share + (1 if i < extra else 0)
            take = min(want, available[year] - quotas[year])
            quotas[year] += take
            placed += take
        if placed == 0:
            break
        remaining -= placed

    return quotas


def select_scaffold_ids(candidate_ids, years, n_taxa, seed):
    """
    Choose `n_taxa` of `candidate_ids` spread evenly over collection year.

    `years` maps a subset of the candidates to their collection year; candidates
    absent from it are not eligible. Returns an unordered set of chosen ids --
    the caller decides the output order.
    """
    rng = random.Random(seed)

    by_year = defaultdict(list)
    for isolate_id in candidate_ids:
        if isolate_id in years:
            by_year[years[isolate_id]].append(isolate_id)
    if not by_year:
        raise ValueError("no alignment sequence has a usable collection date")
    # Sort so the draw depends on the seed alone, not on input order.
    for year in by_year:
        by_year[year].sort()

    quotas = allocate_quotas({y: len(v) for y, v in by_year.items()}, n_taxa)

    chosen = set()
    for year in sorted(by_year):
        if quotas[year]:
            chosen.update(rng.sample(by_year[year], quotas[year]))

    filled = [q for q in quotas.values() if q]
    logger.info(
        f"{len(candidate_ids)} candidates, {len(by_year)} dated years "
        f"({min(by_year)}-{max(by_year)}); drew {len(chosen)} over {len(filled)} "
        f"years, {min(filled)}-{max(filled)} per year"
    )
    if len(chosen) < n_taxa:
        logger.warning(
            f"asked for {n_taxa} scaffold taxa but only {len(chosen)} dated "
            f"sequences are available"
        )
    return chosen


def create_scaffold_alignment(alignment_file, metadata_file, out_alignment,
                              out_samples, n_taxa, seed):
    """Write the scaffold sub-alignment and the list of ids that went into it."""
    with open_sequence_file(alignment_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) < 2:
        raise ValueError("alignment must contain a reference and at least one sequence")

    reference, rest = records[0], records[1:]
    years = read_collection_years(metadata_file, {r.id for r in rest})
    chosen = select_scaffold_ids([r.id for r in rest], years, n_taxa, seed)

    # Input order, not sampling order: a randomised input alignment must yield a
    # randomised scaffold, since usher-sampled places greedily in file order.
    scaffold = [reference] + [r for r in rest if r.id in chosen]

    with open_sequence_file(out_alignment, "wt") as handle:
        SeqIO.write(scaffold, handle, "fasta")
    with open(out_samples, "w") as handle:
        for record in scaffold:
            handle.write(f"{record.id}\n")

    logger.info(
        f"wrote {len(scaffold)} sequences (reference {reference.id} + "
        f"{len(scaffold) - 1} scaffold taxa) to {out_alignment}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Subset an alignment to sequences spread evenly across collection year"
    )
    parser.add_argument("--alignment", required=True,
                        help="Input alignment, reference first (.gz/.xz ok)")
    parser.add_argument("--metadata", required=True,
                        help="CSV with isolate_id and collection_date columns")
    parser.add_argument("--output-alignment", required=True,
                        help="Where to write the scaffold alignment (.gz/.xz ok)")
    parser.add_argument("--output-samples", required=True,
                        help="Where to write the scaffold ids, one per line")
    parser.add_argument("--n-taxa", type=int, required=True,
                        help="How many sequences to draw, excluding the reference")
    parser.add_argument("--seed", type=int, required=True,
                        help="Random seed for the within-year draw")
    return parser.parse_args()


def main():
    args = parse_args()
    create_scaffold_alignment(
        args.alignment,
        args.metadata,
        args.output_alignment,
        args.output_samples,
        args.n_taxa,
        args.seed,
    )


if __name__ == "__main__":
    main()
