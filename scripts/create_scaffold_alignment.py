"""
Cut an alignment down to a small subset of sequences spread evenly over
collection time, to be built into a scaffold tree.

Placing 86,000 sequences onto an empty tree ends up with a high-level shape that
is measurably not the most parsimonious one available: on HA/H3, sequences
collected from 2006 on sit ~70 branches deeper on the trunk than 2005's without a
matching rise in divergence (issue #53). Optimising ~1000 sequences is a far
easier search than optimising 86,000, so the pipeline now builds a backbone from
a time-spread subset first and places everything else onto it.

"Evenly over time" means an equal quota per collection year, not per sequence:
the eight 1963 HA/H3 sequences say far more about the shape of the trunk than
any eight of 2024's 11,786 do. Years that cannot fill their quota hand the
remainder back to be redistributed over the years that still have sequences to
give -- repeatedly, since a redistribution can itself overfill a small year.

Sequences with no collection date are not candidates: they carry no information
about where along the time axis they belong. They are placed later, in the pass
that adds every non-scaffold sequence to the backbone, so nothing is dropped.

An equal quota per year, on its own, drains the sparse years: a year holding 21
sequences against a quota of 17 hands over almost the same 17 under every seed,
so it contributes no diversity between randomisations while still consuming its
slots. `max_year_fraction` caps what any one year may give at that fraction of
what it holds -- but never below one sequence, so no year is dropped for being
small. See "How the scaffold is sampled" in README.md for the measurements.

The first record is taken to be the curated reference and is always kept, at the
top: faToVcf reads the reference off the alignment's first record. Everything
after it is emitted in the order it appeared, so the scaffold reads back in the
same order as the alignment it came from.

What makes each randomisation's backbone different is the seed, not that order.
The draw is deliberately independent of input order -- each year's candidates are
sorted before sampling -- so drawing from a differently shuffled alignment under
the same seed returns exactly the same set. On HA/H3 any two randomisations share
~19% of their scaffold ids, which is where the diversity comes from. (An
earlier version of this docstring said usher-sampled "places sequences greedily
in alignment order", which overstates it: -A is --sort-before-placement-3, so it
sorts by ambiguous-base count and alignment order is only a tie-break.)
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

    `available` maps each year to the number of sequences it actually has, and
    `total` must be non-negative. Returns year -> number to draw, where no value
    exceeds that year's availability and the values sum to min(total, sum of
    availability).

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
        # A round always places something: every open year has room for at least
        # one more, and every open year is asked for at least one -- either
        # share >= 1, or share == 0 and the leftover goes to the first `extra`
        # years as a +1. Assert rather than `break`, because breaking here would
        # quietly return short quotas, and a round that could not advance would
        # spin forever.
        assert placed > 0, "water-filling made no progress"
        remaining -= placed

    return quotas


def select_scaffold_ids(candidate_ids, years, n_taxa, seed, max_year_fraction):
    """
    Choose `n_taxa` of `candidate_ids` spread evenly over collection year.

    `years` maps a subset of the candidates to their collection year; candidates
    absent from it are not eligible. Returns an unordered set of chosen ids --
    the caller decides the output order.

    `max_year_fraction` bounds what any single year may contribute, as a
    fraction of the candidates it holds, floored at one so that a year with a
    single sequence still contributes it. 1.0 imposes no bound, which is the
    behaviour this function had before the parameter existed.
    """
    if n_taxa < 1:
        raise ValueError(f"n_taxa must be at least 1, got {n_taxa}")
    if not 0 < max_year_fraction <= 1:
        raise ValueError(
            f"max_year_fraction must be in (0, 1], got {max_year_fraction}"
        )

    rng = random.Random(seed)

    # A repeated id counts once, or it would take more than one slot out of its
    # year's quota while the caller, matching output records by id, emits every
    # copy anyway -- an overdraw disguised as a full draw.
    duplicated = len(candidate_ids) - len(set(candidate_ids))
    if duplicated:
        logger.warning(
            f"{duplicated} of {len(candidate_ids)} alignment ids are repeats; "
            f"counting each id once and keeping its first record"
        )

    by_year = defaultdict(list)
    seen = set()
    for isolate_id in candidate_ids:
        if isolate_id in years and isolate_id not in seen:
            seen.add(isolate_id)
            by_year[years[isolate_id]].append(isolate_id)
    if not by_year:
        raise ValueError("no alignment sequence has a usable collection date")
    # Sort so the draw depends on the seed alone, not on input order.
    for year in by_year:
        by_year[year].sort()

    # Floored at one: int() alone would zero out every year holding fewer than
    # 1/max_year_fraction sequences, and those are the oldest years -- exactly
    # the ones the scaffold exists to anchor. On HA/H7 and NA/N9 the oldest year
    # holds a single sequence, so int() would have deleted it outright.
    capacity = {y: max(1, int(len(v) * max_year_fraction))
                for y, v in by_year.items()}
    quotas = allocate_quotas(capacity, n_taxa)

    chosen = set()
    for year in sorted(by_year):
        if quotas[year]:
            chosen.update(rng.sample(by_year[year], quotas[year]))

    filled = [q for q in quotas.values() if q]
    logger.info(
        f"{len(candidate_ids)} candidates, {len(by_year)} dated years "
        f"({min(by_year)}-{max(by_year)}); drew {len(chosen)} over {len(filled)} "
        f"years, {min(filled)}-{max(filled)} per year, "
        f"max_year_fraction={max_year_fraction}"
    )
    if len(chosen) < n_taxa:
        n_dated = sum(len(v) for v in by_year.values())
        if n_dated >= n_taxa:
            # The data could have filled the request and the ceiling stopped it.
            # Silently shipping a short scaffold would make n_taxa a lie, so
            # this is the caller's error to fix, in one knob or the other.
            raise ValueError(
                f"max_year_fraction={max_year_fraction} caps the draw at "
                f"{sum(capacity.values())} sequences, short of the {n_taxa} "
                f"requested, though {n_dated} dated sequences are available. "
                f"Raise max_year_fraction or lower n_taxa."
            )
        logger.warning(
            f"asked for {n_taxa} scaffold taxa but only {len(chosen)} dated "
            f"sequences are available"
        )
    return chosen


def create_scaffold_alignment(alignment_file, metadata_file, out_alignment,
                              out_samples, n_taxa, seed, max_year_fraction):
    """Write the scaffold sub-alignment and the list of ids that went into it."""
    with open_sequence_file(alignment_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) < 2:
        raise ValueError("alignment must contain a reference and at least one sequence")

    reference, rest = records[0], records[1:]
    rest_ids = [record.id for record in rest]
    years = read_collection_years(metadata_file, set(rest_ids))
    chosen = select_scaffold_ids(rest_ids, years, n_taxa, seed, max_year_fraction)

    # Input order, not sampling order: a randomised input alignment must yield a
    # randomised scaffold, since usher-sampled places greedily in file order.
    # Each drawn id contributes one record, so a repeated id cannot smuggle a
    # second copy into the scaffold and give usher two samples of the same name.
    pending = set(chosen)
    scaffold = [reference]
    for record in rest:
        if record.id in pending:
            pending.discard(record.id)
            scaffold.append(record)

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
    parser.add_argument("--max-year-fraction", type=float, required=True,
                        help="Most of one year's candidates any draw may take, "
                             "as a fraction; floored at one sequence. 1.0 for "
                             "no bound")
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
        args.max_year_fraction,
    )


if __name__ == "__main__":
    main()
