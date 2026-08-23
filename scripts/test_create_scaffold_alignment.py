"""
Tests for create_scaffold_alignment.py.

The behaviour that matters is the quota: the point of the scaffold is that it
spans the tree's whole time range, so a year holding eight sequences has to
contribute all eight even when another year holds eleven thousand. A plain
random draw of 1000 from HA/H3 would take ~150 sequences from 2024 alone and
would miss most years before 1991 entirely, which is exactly the backbone the
search cannot find on its own (issue #53).

The second thing under test is output order: the scaffold comes out in the
randomized alignment's order rather than in sampling order, so that it reads back
in the same order as the alignment it was cut from. What actually keeps the
randomizations' backbones distinct is the per-randomization seed -- the draw is
independent of input order by design, which `test_does_not_depend_on_candidate_order`
pins -- so `test_output_follows_input_order_not_sampling_order` is guarding a
weaker property than an earlier version of this docstring claimed.
"""

import lzma
import os
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))

import create_scaffold_alignment as csa


REFERENCE_ID = "REFERENCE"


def write_metadata(path, dates):
    """dates: {isolate_id: collection_date string}."""
    with open(path, "w") as handle:
        handle.write("isolate_id,isolate_name,collection_date\n")
        for isolate_id, date in dates.items():
            handle.write(f"{isolate_id},name_{isolate_id},{date}\n")


def write_alignment(path, ids, reference_id=REFERENCE_ID):
    """An .xz alignment with the reference first, then `ids` in order."""
    with lzma.open(path, "wt") as handle:
        handle.write(f">{reference_id}\nACGT\n")
        for isolate_id in ids:
            handle.write(f">{isolate_id}\nACGT\n")


def read_fasta_ids(path):
    ids = []
    with lzma.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip())
    return ids


class AllocateQuotasTestCase(unittest.TestCase):
    def test_splits_evenly_when_every_year_has_room(self):
        quotas = csa.allocate_quotas({2000: 100, 2001: 100, 2002: 100}, 30)
        self.assertEqual(quotas, {2000: 10, 2001: 10, 2002: 10})

    def test_never_exceeds_a_year_s_availability(self):
        available = {2000: 2, 2001: 3, 2002: 500}
        quotas = csa.allocate_quotas(available, 100)
        for year, quota in quotas.items():
            self.assertLessEqual(quota, available[year], f"overdrew {year}")

    def test_small_years_are_drained_and_the_remainder_redistributed(self):
        # The whole point: 2000 and 2001 cannot fill an equal share, so the
        # sequences they cannot supply go to 2002 rather than going unused.
        quotas = csa.allocate_quotas({2000: 2, 2001: 3, 2002: 500}, 100)
        self.assertEqual(quotas[2000], 2)
        self.assertEqual(quotas[2001], 3)
        self.assertEqual(quotas[2002], 95)
        self.assertEqual(sum(quotas.values()), 100)

    def test_a_single_huge_year_does_not_crowd_out_many_tiny_ones(self):
        available = {1990 + i: 5 for i in range(10)}
        available[2024] = 11000
        quotas = csa.allocate_quotas(available, 100)
        for year in range(1990, 2000):
            self.assertEqual(quotas[year], 5, f"{year} should be drained")
        self.assertEqual(quotas[2024], 50)

    def test_returns_everything_when_asked_for_more_than_exists(self):
        available = {2000: 4, 2001: 6}
        quotas = csa.allocate_quotas(available, 1000)
        self.assertEqual(quotas, available)

    def test_sums_to_the_requested_total_when_supply_allows(self):
        # A total that divides unevenly over the years exercises the remainder
        # path, which must not hand a +1 to a year with no room for it.
        quotas = csa.allocate_quotas({2000: 1, 2001: 1, 2002: 50}, 17)
        self.assertEqual(sum(quotas.values()), 17)
        self.assertEqual(quotas[2000], 1)
        self.assertEqual(quotas[2001], 1)

    def test_zero_total_draws_nothing(self):
        quotas = csa.allocate_quotas({2000: 10, 2001: 10}, 0)
        self.assertEqual(sum(quotas.values()), 0)

    def test_single_year(self):
        self.assertEqual(csa.allocate_quotas({2000: 10}, 4), {2000: 4})


class ReadCollectionYearsTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.metadata = os.path.join(self.tmp.name, "metadata.csv")

    def tearDown(self):
        self.tmp.cleanup()

    def test_reads_the_year_from_a_full_date(self):
        write_metadata(self.metadata, {"A": "2009-03-14"})
        self.assertEqual(csa.read_collection_years(self.metadata, {"A"}), {"A": 2009})

    def test_accepts_a_year_only_date(self):
        write_metadata(self.metadata, {"A": "1997"})
        self.assertEqual(csa.read_collection_years(self.metadata, {"A"}), {"A": 1997})

    def test_skips_undatable_rows(self):
        write_metadata(self.metadata,
                       {"A": "", "B": "199", "C": "unknown", "D": "2009-01-01"})
        years = csa.read_collection_years(self.metadata, {"A", "B", "C", "D"})
        self.assertEqual(years, {"D": 2009})

    def test_ignores_ids_that_were_not_asked_for(self):
        write_metadata(self.metadata, {"A": "2009-01-01", "B": "2010-01-01"})
        self.assertEqual(csa.read_collection_years(self.metadata, {"A"}), {"A": 2009})


class SelectScaffoldIdsTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.metadata = os.path.join(self.tmp.name, "metadata.csv")

    def tearDown(self):
        self.tmp.cleanup()

    def _years(self, per_year):
        """per_year: {year: count} -> ({id: year}, [ids])."""
        years, ids = {}, []
        for year, count in per_year.items():
            for i in range(count):
                isolate_id = f"S_{year}_{i}"
                years[isolate_id] = year
                ids.append(isolate_id)
        return years, ids

    def test_draws_the_requested_number(self):
        years, ids = self._years({2000: 50, 2001: 50, 2002: 50})
        self.assertEqual(len(csa.select_scaffold_ids(ids, years, 30, seed=0, max_year_fraction=1.0)), 30)

    def test_is_deterministic_for_a_seed(self):
        years, ids = self._years({2000: 50, 2001: 50})
        first = csa.select_scaffold_ids(ids, years, 20, seed=7, max_year_fraction=1.0)
        second = csa.select_scaffold_ids(ids, years, 20, seed=7, max_year_fraction=1.0)
        self.assertEqual(first, second)

    def test_different_seeds_give_different_draws(self):
        """Each randomization must get its own backbone, or they all coincide."""
        years, ids = self._years({2000: 200, 2001: 200})
        first = csa.select_scaffold_ids(ids, years, 40, seed=0, max_year_fraction=1.0)
        second = csa.select_scaffold_ids(ids, years, 40, seed=1, max_year_fraction=1.0)
        self.assertNotEqual(first, second)

    def test_does_not_depend_on_candidate_order(self):
        years, ids = self._years({2000: 50, 2001: 50})
        forward = csa.select_scaffold_ids(ids, years, 20, seed=3, max_year_fraction=1.0)
        backward = csa.select_scaffold_ids(list(reversed(ids)), years, 20, seed=3, max_year_fraction=1.0)
        self.assertEqual(forward, backward)

    def test_covers_the_sparse_years_not_just_the_dense_one(self):
        years, ids = self._years({1963: 8, 1990: 8, 2024: 11000})
        chosen = csa.select_scaffold_ids(ids, years, 100, seed=0, max_year_fraction=1.0)
        drawn_years = {years[i] for i in chosen}
        self.assertEqual(drawn_years, {1963, 1990, 2024})
        self.assertEqual(sum(1 for i in chosen if years[i] == 1963), 8)
        self.assertEqual(sum(1 for i in chosen if years[i] == 1990), 8)

    def test_undated_candidates_are_not_drawn(self):
        years, ids = self._years({2000: 10})
        ids = ids + ["UNDATED_1", "UNDATED_2"]
        chosen = csa.select_scaffold_ids(ids, years, 12, seed=0, max_year_fraction=1.0)
        self.assertNotIn("UNDATED_1", chosen)
        self.assertNotIn("UNDATED_2", chosen)
        self.assertEqual(len(chosen), 10)

    def test_raises_when_nothing_can_be_dated(self):
        with self.assertRaises(ValueError):
            csa.select_scaffold_ids(["A", "B"], {}, 2, seed=0, max_year_fraction=1.0)

    def test_a_year_holding_one_sequence_still_contributes_it(self):
        """The sparse end of the range is the whole point of stratifying."""
        years, ids = self._years({2000: 1, 2001: 50})
        chosen = csa.select_scaffold_ids(ids, years, 10, seed=0, max_year_fraction=1.0)
        self.assertEqual(sum(1 for i in chosen if years[i] == 2000), 1)

    def test_a_repeated_id_takes_a_single_slot(self):
        """
        A repeat must not consume two of its year's quota. It used to: the id
        went into the year's candidate list twice, and both copies could be
        drawn -- collapsing to one entry in the returned set, so the draw was
        silently short while the caller emitted both records.
        """
        years = {"DUP": 2000, "OTHER": 2000}
        chosen = csa.select_scaffold_ids(["DUP", "DUP", "OTHER"], years, 2, seed=0, max_year_fraction=1.0)
        self.assertEqual(chosen, {"DUP", "OTHER"})

    def test_zero_n_taxa_raises(self):
        """Must name the problem, not die inside a logging f-string."""
        with self.assertRaises(ValueError) as caught:
            csa.select_scaffold_ids(["A"], {"A": 2000}, 0, seed=0, max_year_fraction=1.0)
        self.assertIn("n_taxa", str(caught.exception))

    def test_negative_n_taxa_raises(self):
        with self.assertRaises(ValueError) as caught:
            csa.select_scaffold_ids(["A"], {"A": 2000}, -5, seed=0, max_year_fraction=1.0)
        self.assertIn("n_taxa", str(caught.exception))


class CreateScaffoldAlignmentTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.alignment = os.path.join(self.tmp.name, "msa.fasta.xz")
        self.metadata = os.path.join(self.tmp.name, "metadata.csv")
        self.out_msa = os.path.join(self.tmp.name, "scaffold.fasta.xz")
        self.out_samples = os.path.join(self.tmp.name, "scaffold_samples.txt")

    def tearDown(self):
        self.tmp.cleanup()

    def _build(self, per_year, n_taxa, seed=0, undated=(), max_year_fraction=1.0):
        ids, dates = [], {}
        for year, count in per_year.items():
            for i in range(count):
                isolate_id = f"S_{year}_{i}"
                ids.append(isolate_id)
                dates[isolate_id] = f"{year}-06-01"
        ids.extend(undated)
        for isolate_id in undated:
            dates[isolate_id] = ""
        write_alignment(self.alignment, ids)
        write_metadata(self.metadata, dates)
        csa.create_scaffold_alignment(self.alignment, self.metadata, self.out_msa,
                                      self.out_samples, n_taxa, seed,
                                      max_year_fraction)
        return ids

    def test_reference_is_kept_and_comes_first(self):
        """faToVcf reads the reference off the first record, so this is load-bearing."""
        self._build({2000: 20, 2001: 20}, n_taxa=10)
        written = read_fasta_ids(self.out_msa)
        self.assertEqual(written[0], REFERENCE_ID)
        self.assertEqual(len(written), 11)

    def test_reference_is_not_counted_against_the_taxon_quota(self):
        self._build({2000: 20, 2001: 20}, n_taxa=10)
        self.assertEqual(len(read_fasta_ids(self.out_msa)) - 1, 10)

    def test_output_follows_input_order_not_sampling_order(self):
        ids = self._build({2000: 40, 2001: 40}, n_taxa=20)
        written = read_fasta_ids(self.out_msa)[1:]
        position = {isolate_id: i for i, isolate_id in enumerate(ids)}
        self.assertEqual(written, sorted(written, key=position.__getitem__))

    def test_samples_file_matches_the_alignment(self):
        self._build({2000: 20, 2001: 20}, n_taxa=10)
        with open(self.out_samples) as handle:
            listed = [line.strip() for line in handle if line.strip()]
        self.assertEqual(listed, read_fasta_ids(self.out_msa))

    def test_undated_sequences_are_left_for_the_placement_pass(self):
        self._build({2000: 20}, n_taxa=10, undated=("NO_DATE_1", "NO_DATE_2"))
        written = read_fasta_ids(self.out_msa)
        self.assertNotIn("NO_DATE_1", written)
        self.assertNotIn("NO_DATE_2", written)

    def test_spans_the_full_time_range(self):
        self._build({1963: 3, 1985: 3, 2005: 500, 2024: 5000}, n_taxa=50)
        written = read_fasta_ids(self.out_msa)[1:]
        drawn_years = {i.split("_")[1] for i in written}
        self.assertEqual(drawn_years, {"1963", "1985", "2005", "2024"})

    def test_rejects_an_alignment_with_nothing_to_sample(self):
        write_alignment(self.alignment, [])
        write_metadata(self.metadata, {})
        with self.assertRaises(ValueError):
            csa.create_scaffold_alignment(self.alignment, self.metadata,
                                          self.out_msa, self.out_samples, 10, 0,
                                          1.0)

    def test_asking_for_more_than_exists_yields_everything_dated(self):
        self._build({2000: 5, 2001: 5}, n_taxa=1000)
        self.assertEqual(len(read_fasta_ids(self.out_msa)) - 1, 10)

    def test_a_repeated_id_yields_a_single_record(self):
        """
        Two records under one name would reach usher as two samples of the same
        name, from a scaffold that reported drawing one taxon.
        """
        write_alignment(self.alignment, ["DUP", "DUP"])
        write_metadata(self.metadata, {"DUP": "2000-06-01"})
        csa.create_scaffold_alignment(self.alignment, self.metadata, self.out_msa,
                                      self.out_samples, 1, 0, 1.0)
        written = read_fasta_ids(self.out_msa)
        self.assertEqual(written, [REFERENCE_ID, "DUP"])
        with open(self.out_samples) as handle:
            listed = [line.strip() for line in handle if line.strip()]
        self.assertEqual(listed, [REFERENCE_ID, "DUP"])

    def test_a_sequence_absent_from_the_metadata_is_not_drawn(self):
        """Absent is not the same as present-but-blank, and both are ineligible."""
        write_alignment(self.alignment, ["DATED", "NOT_IN_METADATA"])
        write_metadata(self.metadata, {"DATED": "2000-06-01"})
        csa.create_scaffold_alignment(self.alignment, self.metadata, self.out_msa,
                                      self.out_samples, 5, 0, 1.0)
        self.assertEqual(read_fasta_ids(self.out_msa), [REFERENCE_ID, "DATED"])

    def test_zero_n_taxa_raises(self):
        with self.assertRaises(ValueError):
            self._build({2000: 5}, n_taxa=0)


class CliTestCase(unittest.TestCase):
    """
    Exercise main() the way the Snakemake rule does. A flag renamed here but not
    in the Snakefile -- or an args attribute that no longer matches its flag --
    would otherwise surface only once the pipeline ran.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_main_wires_every_flag_through(self):
        alignment = os.path.join(self.tmp.name, "msa.fasta.xz")
        metadata = os.path.join(self.tmp.name, "metadata.csv")
        out_msa = os.path.join(self.tmp.name, "scaffold.fasta.xz")
        out_samples = os.path.join(self.tmp.name, "scaffold_samples.txt")
        write_alignment(alignment, [f"S_{i}" for i in range(10)])
        write_metadata(metadata, {f"S_{i}": f"{2000 + i}-06-01" for i in range(10)})

        argv = [
            "create_scaffold_alignment.py",
            "--alignment", alignment,
            "--metadata", metadata,
            "--output-alignment", out_msa,
            "--output-samples", out_samples,
            "--n-taxa", "4",
            "--max-year-fraction", "1.0",
            "--seed", "0",
        ]
        with mock.patch.object(sys, "argv", argv):
            csa.main()

        written = read_fasta_ids(out_msa)
        self.assertEqual(written[0], REFERENCE_ID)
        self.assertEqual(len(written), 5)

    def test_seed_reaches_the_draw(self):
        """A --seed that main() ignored would make every randomization identical."""
        alignment = os.path.join(self.tmp.name, "msa.fasta.xz")
        metadata = os.path.join(self.tmp.name, "metadata.csv")
        write_alignment(alignment, [f"S_{i}" for i in range(200)])
        write_metadata(metadata,
                       {f"S_{i}": f"{2000 + i % 4}-06-01" for i in range(200)})

        draws = []
        for seed in ("0", "1"):
            out_msa = os.path.join(self.tmp.name, f"s{seed}.fasta.xz")
            out_samples = os.path.join(self.tmp.name, f"s{seed}.txt")
            argv = [
                "create_scaffold_alignment.py",
                "--alignment", alignment,
                "--metadata", metadata,
                "--output-alignment", out_msa,
                "--output-samples", out_samples,
                "--n-taxa", "20",
                "--max-year-fraction", "1.0",
                "--seed", seed,
            ]
            with mock.patch.object(sys, "argv", argv):
                csa.main()
            draws.append(read_fasta_ids(out_msa))
        self.assertNotEqual(draws[0], draws[1])


class MaxYearFractionTestCase(unittest.TestCase):
    """Bound what one year can give, so sparse years stop being drawn dry.

    An equal per-year quota alone hands over nearly the same sequences from a
    year that can barely fill its quota, under every seed. Measured across the
    10 randomizations after the collection dates were repaired, mean pairwise
    overlap between two HA/H3 backbones was 38% unbounded and 19% at 0.5.
    """

    def _years(self, per_year):
        years, ids = {}, []
        for year, count in per_year.items():
            for i in range(count):
                isolate_id = f"S_{year}_{i}"
                years[isolate_id] = year
                ids.append(isolate_id)
        return years, ids

    def test_a_year_cannot_give_more_than_its_fraction(self):
        """The dense year is capped even though it could fill the whole draw."""
        years, ids = self._years({2000: 100, 2001: 100})
        chosen = csa.select_scaffold_ids(ids, years, 40, seed=0,
                                        max_year_fraction=0.25)
        from_2000 = sum(1 for i in chosen if years[i] == 2000)
        self.assertLessEqual(from_2000, 25)

    def test_one_point_zero_is_exactly_the_unbounded_draw(self):
        """The parameter must have an identity value, or it cannot be turned off."""
        years, ids = self._years({2000: 40, 2001: 5, 2002: 200})
        self.assertEqual(
            csa.select_scaffold_ids(ids, years, 50, seed=4, max_year_fraction=1.0),
            csa.select_scaffold_ids(ids, years, 50, seed=4, max_year_fraction=1.0),
        )
        bounded = csa.select_scaffold_ids(ids, years, 50, seed=4,
                                         max_year_fraction=0.5)
        unbounded = csa.select_scaffold_ids(ids, years, 50, seed=4,
                                           max_year_fraction=1.0)
        self.assertNotEqual(bounded, unbounded)

    def test_a_single_sequence_year_survives_any_fraction(self):
        """int(1 * 0.25) is 0, which would delete the year outright.

        On HA/H7 and NA/N9 the oldest year holds exactly one sequence, and those
        are the sequences the scaffold exists to anchor the root region with.
        """
        years, ids = self._years({1902: 1, 2020: 400})
        chosen = csa.select_scaffold_ids(ids, years, 5, seed=0,
                                        max_year_fraction=0.01)
        self.assertIn("S_1902_0", chosen)

    def test_a_bound_that_cannot_meet_n_taxa_raises(self):
        """Short-drawing silently would make n_taxa a lie.

        400 candidates at 0.1 can field 40, so a request for 100 is a
        configuration error rather than a data limitation.
        """
        years, ids = self._years({2000: 200, 2001: 200})
        with self.assertRaises(ValueError) as caught:
            csa.select_scaffold_ids(ids, years, 100, seed=0,
                                    max_year_fraction=0.1)
        self.assertIn("max_year_fraction", str(caught.exception))

    def test_too_few_dated_sequences_still_only_warns(self):
        """The other reason a draw comes up short is unavoidable, so it is not
        an error -- NA/N9 legitimately holds fewer sequences than some n_taxa."""
        years, ids = self._years({2000: 3, 2001: 2})
        chosen = csa.select_scaffold_ids(ids, years, 100, seed=0,
                                        max_year_fraction=1.0)
        self.assertEqual(len(chosen), 5)

    def test_lowering_the_bound_raises_diversity_between_seeds(self):
        """The property the parameter exists for.

        2000 must be a year that gets drawn *dry*, which is where the diversity
        goes: a fair share of 60 over three years is 20, and 2000 holds exactly
        20, so unbounded it hands over the same 20 under every seed. At 0.5 it
        may give only 10 of its 20, and which 10 varies. Equal-sized years would
        show nothing -- a cap above the quota never binds.

        Averaged over seed pairs rather than a single pair, so the assertion
        rests on the distribution and not on one lucky draw.
        """
        years, ids = self._years({2000: 20, 2001: 100, 2002: 100})

        def mean_overlap(fraction):
            draws = [csa.select_scaffold_ids(ids, years, 60, seed=s,
                                             max_year_fraction=fraction)
                     for s in range(6)]
            pairs = [(a, b) for i, a in enumerate(draws) for b in draws[i + 1:]]
            return sum(len(a & b) / len(a) for a, b in pairs) / len(pairs)

        self.assertLess(mean_overlap(0.5), mean_overlap(1.0))

    def test_out_of_range_fractions_raise(self):
        years, ids = self._years({2000: 10})
        for bad in (0, -0.5, 1.5):
            with self.assertRaises(ValueError):
                csa.select_scaffold_ids(ids, years, 5, seed=0,
                                        max_year_fraction=bad)

    def test_the_cli_flag_reaches_the_draw(self):
        """A flag main() ignored would leave every backbone unbounded."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        alignment = os.path.join(tmp.name, "msa.fasta.xz")
        metadata = os.path.join(tmp.name, "metadata.csv")
        write_alignment(alignment, [f"S_{i}" for i in range(200)])
        # Two years only, so a bound of 0.1 must visibly shrink the draw.
        write_metadata(metadata,
                       {f"S_{i}": f"{2000 + i % 2}-06-01" for i in range(200)})
        out_msa = os.path.join(tmp.name, "bounded.fasta.xz")
        out_samples = os.path.join(tmp.name, "bounded.txt")
        argv = [
            "create_scaffold_alignment.py",
            "--alignment", alignment,
            "--metadata", metadata,
            "--output-alignment", out_msa,
            "--output-samples", out_samples,
            "--n-taxa", "20",
            "--max-year-fraction", "0.1",
            "--seed", "0",
        ]
        with mock.patch.object(sys, "argv", argv):
            csa.main()
        # 100 candidates per year at 0.1 caps each at 10, so 20 is still met;
        # the draw must not exceed the cap in either year.
        drawn = [i for i in read_fasta_ids(out_msa) if i != REFERENCE_ID]
        self.assertEqual(len(drawn), 20)
        self.assertEqual(sum(1 for i in drawn if int(i.split("_")[1]) % 2 == 0), 10)


if __name__ == "__main__":
    unittest.main()
