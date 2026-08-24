"""
Tests for parse_gisaid_data.py.

Two behaviours are under test:

1. Only configured subtypes are written. The script previously produced 30
   segment/subtype directories against 16 configured combinations.
2. main() returns a nonzero exit code when a configured combination ends up
   empty or no metadata is written. It previously had no return statement at
   all, so `sys.exit(main())` was always `sys.exit(None)` -- exit 0 no matter
   how much was dropped.
3. collection_date survives at the precision GISAID supplied, per metadata
   file. A bare pd.to_datetime(errors='coerce') used to infer one format per
   file from its first element and empty every value that did not match it
   (issue #55).

pandas.read_excel is patched rather than writing real .xls fixtures; the Excel
parsing itself is not what these tests are about.
"""

import lzma
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(__file__))

import parse_gisaid_data


# GISAID sequence IDs are EPI|SEGMENT|NAME|EPI_ISL|SUBTYPE
def seq_id(segment, epi_isl, subtype):
    return f"EPI|{segment}|A/test/1/2020|{epi_isl}|{subtype}"


METADATA_COLUMNS = [
    "Isolate_Id", "Isolate_Name", "Subtype", "Clade",
    "Passage_History", "Location", "Host", "Collection_Date",
]


def fake_metadata(isolate_ids, dates=None):
    """One row per id. `dates` supplies Collection_Date in the same order.

    Row order matters: pandas infers a date format from the *first* element, so
    a test for issue #55 has to control which value comes first.
    """
    if dates is None:
        dates = ["2020-01-01"] * len(isolate_ids)
    if len(dates) != len(isolate_ids):
        raise ValueError("dates and isolate_ids must be the same length")
    return pd.DataFrame([
        {
            "Isolate_Id": i, "Isolate_Name": f"A/test/{i}/2020",
            "Subtype": "A / H1N1", "Clade": "x", "Passage_History": "Original",
            "Location": "North America / USA", "Host": "Human",
            "Collection_Date": d,
        }
        for i, d in zip(isolate_ids, dates)
    ], columns=METADATA_COLUMNS)


def run_parse(tmpdir, records, isolate_ids, ha_subtypes, na_subtypes,
              segments=("HA", "NA", "PB2"), with_metadata=True,
              dates=None, metadata_frames=None):
    """Write inputs, run main(), return (exit_code, output_dir).

    `metadata_frames` supplies one DataFrame per metadata file, consumed in
    sorted filename order -- the script sorts its glob, so frame i is always
    read from meta_{i:02d}.xls. The zero-padding is what keeps that true past
    ten files, where "meta_10" would otherwise sort before "meta_2". Use it to
    test that each file is parsed on its own. Otherwise a single file is built
    from `isolate_ids` and `dates`.
    """
    data_dir = os.path.join(tmpdir, "data")
    os.makedirs(data_dir, exist_ok=True)
    with open(os.path.join(data_dir, "seqs.fasta"), "w") as f:
        for rid in records:
            f.write(f">{rid}\nACGTACGTACGT\n")
    if with_metadata:
        # Content is supplied by the read_excel patch; the files must merely
        # exist for the script's glob to find them.
        n_files = 1 if metadata_frames is None else len(metadata_frames)
        for i in range(n_files):
            open(os.path.join(data_dir, f"meta_{i:02d}.xls"), "w").close()

    out_dir = os.path.join(tmpdir, "results")
    os.makedirs(out_dir, exist_ok=True)
    argv = [
        "parse_gisaid_data.py",
        "--input-dirs", data_dir,
        "--output-dir", out_dir,
        "--segments", *segments,
        "--ha-subtypes", *ha_subtypes,
        "--na-subtypes", *na_subtypes,
    ]
    if metadata_frames is None:
        excel_patch = mock.patch.object(
            pd, "read_excel", return_value=fake_metadata(isolate_ids, dates))
    else:
        excel_patch = mock.patch.object(
            pd, "read_excel", side_effect=list(metadata_frames))
    with mock.patch.object(sys, "argv", argv), excel_patch:
        code = parse_gisaid_data.main()
    return code, out_dir


def read_collection_dates(out_dir):
    """isolate_id -> collection_date as written, with blanks as ''."""
    df = pd.read_csv(os.path.join(out_dir, "combined_metadata.csv"), dtype=str)
    return {
        row.isolate_id: ("" if pd.isna(row.collection_date) else row.collection_date)
        for row in df.itertuples()
    }


def written_combinations(out_dir):
    found = set()
    for segment in os.listdir(out_dir):
        seg_path = os.path.join(out_dir, segment)
        if not os.path.isdir(seg_path):
            continue
        for group in os.listdir(seg_path):
            if os.path.exists(os.path.join(seg_path, group, "raw_sequences.fasta.xz")):
                found.add((segment, group))
    return found


def read_records(out_dir, segment, group):
    path = os.path.join(out_dir, segment, group, "raw_sequences.fasta.xz")
    with lzma.open(path, "rt") as f:
        return [r.id for r in SeqIO.parse(f, "fasta")]


class TestSubtypeAllowList(unittest.TestCase):
    """Only configured subtypes get an output directory."""

    def test_unconfigured_ha_subtype_is_not_written(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("HA", "EPI_ISL_2", "A_/_H7N9"),   # H7 not configured
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, records, ["EPI_ISL_1", "EPI_ISL_2"],
                                  ha_subtypes=["H1"], na_subtypes=["N1"])
            self.assertEqual(code, 0)
            self.assertEqual(written_combinations(out),
                             {("HA", "H1"), ("NA", "N1"), ("PB2", "all")})

    def test_unconfigured_na_subtype_is_not_written(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_2", "A_/_H7N9"),   # N9 not configured
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, records, ["EPI_ISL_1", "EPI_ISL_2"],
                                  ha_subtypes=["H1"], na_subtypes=["N1"])
            self.assertEqual(code, 0)
            self.assertNotIn(("NA", "N9"), written_combinations(out))

    def test_internal_segments_ignore_subtype(self):
        """Internal segments group under 'all' regardless of subtype."""
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_2", "A_/_H7N9"),  # unconfigured subtype, still kept
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, records, ["EPI_ISL_1", "EPI_ISL_2"],
                                  ha_subtypes=["H1"], na_subtypes=["N1"])
            self.assertEqual(code, 0)
            self.assertEqual(sorted(read_records(out, "PB2", "all")),
                             ["EPI_ISL_1", "EPI_ISL_2"])


class TestFailsLoudly(unittest.TestCase):
    """main() returns nonzero instead of silently succeeding."""

    def test_configured_combination_with_no_records_exits_nonzero(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            # nothing for PB2, which is a configured segment
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, _ = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                ha_subtypes=["H1"], na_subtypes=["N1"])
        self.assertEqual(code, 1)

    def test_configured_subtype_absent_from_data_exits_nonzero(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            # H3 is configured but absent from the data
            code, _ = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                ha_subtypes=["H1", "H3"], na_subtypes=["N1"])
        self.assertEqual(code, 1)

    def test_no_metadata_exits_nonzero(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, _ = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                ha_subtypes=["H1"], na_subtypes=["N1"],
                                with_metadata=False)
        self.assertEqual(code, 1)

    def test_all_records_dropped_exits_nonzero(self):
        """The worst case the old code reported as success."""
        records = [seq_id("HA", "EPI_ISL_1", "A_/_H7N9")]  # nothing configured matches
        with tempfile.TemporaryDirectory() as tmpdir:
            code, _ = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                ha_subtypes=["H1"], na_subtypes=["N1"])
        self.assertEqual(code, 1)

    def test_happy_path_exits_zero(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, _ = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                ha_subtypes=["H1"], na_subtypes=["N1"])
        self.assertEqual(code, 0)


class TestDuplicateHandling(unittest.TestCase):
    """Duplicates are legitimate -- overlapping download windows -- and must be
    reported rather than treated as an error."""

    def test_duplicate_epi_isl_is_dropped_but_not_fatal(self):
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),   # same id, same combination
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                  ha_subtypes=["H1"], na_subtypes=["N1"])
            self.assertEqual(code, 0)
            self.assertEqual(read_records(out, "HA", "H1"), ["EPI_ISL_1"])

    def test_same_id_in_different_segments_is_kept(self):
        """Dedup is per segment/subtype, not global."""
        records = [
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                  ha_subtypes=["H1"], na_subtypes=["N1"])
            self.assertEqual(code, 0)
            for segment, group in [("HA", "H1"), ("NA", "N1"), ("PB2", "all")]:
                self.assertEqual(read_records(out, segment, group), ["EPI_ISL_1"])

    def test_unparseable_id_is_dropped_but_not_fatal(self):
        records = [
            "not-a-gisaid-id",
            seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
            seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, records, ["EPI_ISL_1"],
                                  ha_subtypes=["H1"], na_subtypes=["N1"])
            self.assertEqual(code, 0)
            self.assertEqual(read_records(out, "HA", "H1"), ["EPI_ISL_1"])


class TestUnparseableSubtype(unittest.TestCase):
    """A subtype that is not H*N* is skipped, not crashed on and not guessed.

    extract_ha_na_subtype() returns (None, None) for anything its regex misses,
    and HA and NA records take separate branches on that. Real GISAID data
    carries entries like "B / Victoria" that land here.

    Every case supplies a full complement of parseable records, because
    parse_gisaid_data exits nonzero if any configured combination produces
    nothing -- a different guard, already covered elsewhere.
    """

    GOOD = [
        seq_id("HA", "EPI_ISL_1", "A_/_H1N1"),
        seq_id("NA", "EPI_ISL_1", "A_/_H1N1"),
        seq_id("PB2", "EPI_ISL_1", "A_/_H1N1"),
    ]
    IDS = ["EPI_ISL_1", "EPI_ISL_2"]

    def run_with(self, extra):
        tmpdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmpdir, True)
        return run_parse(tmpdir, self.GOOD + [extra], self.IDS,
                         ha_subtypes=["H1"], na_subtypes=["N1"])

    def test_ha_record_with_unparseable_subtype_is_skipped(self):
        code, out = self.run_with(seq_id("HA", "EPI_ISL_2", "B_/_Victoria"))
        self.assertEqual(code, 0)
        self.assertEqual(read_records(out, "HA", "H1"), ["EPI_ISL_1"])

    def test_na_record_with_unparseable_subtype_is_skipped(self):
        code, out = self.run_with(seq_id("NA", "EPI_ISL_2", "B_/_Yamagata"))
        self.assertEqual(code, 0)
        self.assertEqual(read_records(out, "NA", "N1"), ["EPI_ISL_1"])

    def test_partial_subtype_is_not_half_matched(self):
        """"H3" alone has no N part, so the regex must reject the whole thing.

        H3 is deliberately CONFIGURED here. With it unconfigured, a regex that
        half-matched and returned ("H3", None) would be dropped by the
        unrelated "subtype not configured" branch, producing exactly the same
        observable result as correct rejection -- so the test could not tell
        the two apart. Configuring it means a half-match would visibly write
        an HA/H3 record.
        """
        tmpdir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmpdir, True)
        records = self.GOOD + [
            seq_id("HA", "EPI_ISL_3", "A_/_H3N2"),      # keeps H3 non-empty
            seq_id("HA", "EPI_ISL_2", "A_/_H3"),        # the partial one
        ]
        code, out = run_parse(tmpdir, records,
                              ["EPI_ISL_1", "EPI_ISL_2", "EPI_ISL_3"],
                              ha_subtypes=["H1", "H3"], na_subtypes=["N1"])
        self.assertEqual(code, 0)
        self.assertEqual(read_records(out, "HA", "H3"), ["EPI_ISL_3"],
                         "EPI_ISL_2 was half-matched into H3 on a subtype "
                         "with no N part")

    def test_internal_segment_ignores_the_subtype_entirely(self):
        """Internal segments group as 'all', so the H*N* parse never applies."""
        code, out = self.run_with(seq_id("PB2", "EPI_ISL_2", "B_/_Victoria"))
        self.assertEqual(code, 0)
        self.assertEqual(sorted(read_records(out, "PB2", "all")),
                         ["EPI_ISL_1", "EPI_ISL_2"])


class TestCollectionDatePrecision(unittest.TestCase):
    """collection_date is kept verbatim at GISAID's own precision (issue #55).

    The old code called pd.to_datetime(errors='coerce') once per metadata file.
    pandas >= 2.0 infers a single format from the first non-null element, so one
    row decided the fate of that file's whole column: on the real
    data/H1N1/H1N1-window-1.xls, whose first row is the year "2009", all 17,908
    fully-resolved dates in the file became blank. Across five such files 42,619
    valid dates were destroyed, and errors='coerce' meant nothing warned.

    Every case supplies records for all three configured combinations, because
    an empty combination trips a separate nonzero-exit guard.
    """

    def segments_for(self, isolate_ids):
        return [seq_id(seg, i, "A_/_H1N1")
                for seg in ("HA", "NA", "PB2") for i in isolate_ids]

    def test_year_only_first_row_does_not_destroy_later_full_dates(self):
        """The real regression, with values taken from the real .xls.

        Rows 0 and 1 of data/H1N1/H1N1-window-1.xls are EPI_ISL_65694 dated
        "2009" and EPI_ISL_164609 (A/swine/BinhDuong/01-12/2010) dated
        "2010-03-01". Under the old code the leading year fixed the format at
        %Y and EPI_ISL_164609 came out blank in combined_metadata.csv.
        """
        ids = ["EPI_ISL_65694", "EPI_ISL_164609"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, self.segments_for(ids), ids,
                                  ha_subtypes=["H1"], na_subtypes=["N1"],
                                  dates=["2009", "2010-03-01"])
            self.assertEqual(code, 0)
            dates = read_collection_dates(out)
            self.assertEqual(dates["EPI_ISL_164609"], "2010-03-01")
            self.assertEqual(dates["EPI_ISL_65694"], "2009")

    def test_each_metadata_file_is_parsed_independently(self):
        """Two files, each led by a different precision, both fully preserved.

        The old call sat inside the per-file loop, so file A's leading "2022-11"
        emptied its own full date and file B's leading full date emptied its
        year. Nothing crossed between files, which is why the bug hit exactly
        five of 56 files.
        """
        frame_a = fake_metadata(["EPI_ISL_1", "EPI_ISL_2"],
                                dates=["2022-11", "2022-11-15"])
        frame_b = fake_metadata(["EPI_ISL_3", "EPI_ISL_4"],
                                dates=["2020-01-01", "1983"])
        ids = ["EPI_ISL_1", "EPI_ISL_2", "EPI_ISL_3", "EPI_ISL_4"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, self.segments_for(ids), ids,
                                  ha_subtypes=["H1"], na_subtypes=["N1"],
                                  metadata_frames=[frame_a, frame_b])
            self.assertEqual(code, 0)
            dates = read_collection_dates(out)
            self.assertEqual(dates, {
                "EPI_ISL_1": "2022-11", "EPI_ISL_2": "2022-11-15",
                "EPI_ISL_3": "2020-01-01", "EPI_ISL_4": "1983",
            })
            self.assertNotIn("", dates.values())

    def test_partial_dates_are_not_padded_to_a_day(self):
        """Guards against format='mixed' creeping back in.

        format='mixed' would recover all 42,619 destroyed dates, but it also
        maps "2009" to 2009-01-01 and "2022-11" to 2022-11-01, inventing a day
        for 47,090 rows. Consumers that must drop partial dates could then no
        longer tell them apart, so the padding is asserted against explicitly.
        """
        ids = ["EPI_ISL_1", "EPI_ISL_2", "EPI_ISL_3"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, self.segments_for(ids), ids,
                                  ha_subtypes=["H1"], na_subtypes=["N1"],
                                  dates=["2009", "2022-11", "2010-03-01"])
            self.assertEqual(code, 0)
            dates = read_collection_dates(out)
            self.assertEqual(dates["EPI_ISL_1"], "2009")
            self.assertEqual(dates["EPI_ISL_2"], "2022-11")
            self.assertEqual(dates["EPI_ISL_3"], "2010-03-01")

    def test_blank_collection_date_is_kept_blank_and_not_fatal(self):
        """GISAID supplies a date for every isolate today, but a blank is
        missing data rather than a format change, so it must not fail the run.

        Unlike its neighbours this one passes against the pre-fix code too: the
        old coerce-to-NaT path handled blanks correctly. It pins the invariant
        down so the new validation cannot start rejecting them; it is not
        evidence that issue #55 is fixed.
        """
        ids = ["EPI_ISL_1", "EPI_ISL_2"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, self.segments_for(ids), ids,
                                  ha_subtypes=["H1"], na_subtypes=["N1"],
                                  dates=["2010-03-01", ""])
            self.assertEqual(code, 0)
            self.assertEqual(read_collection_dates(out)["EPI_ISL_2"], "")

    def test_malformed_collection_date_exits_nonzero(self):
        """A value matching none of the three precisions means GISAID changed
        what it ships, which is exactly what errors='coerce' used to hide."""
        ids = ["EPI_ISL_1", "EPI_ISL_2"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, _ = run_parse(tmpdir, self.segments_for(ids), ids,
                                ha_subtypes=["H1"], na_subtypes=["N1"],
                                dates=["2010-03-01", "03/04/2020"])
            self.assertEqual(code, 1)

    def test_a_date_that_no_calendar_contains_exits_nonzero(self):
        """Right shape, impossible day.

        A shape-only check would pass "2019-02-30" here and leave it to a
        downstream strict %Y-%m-%d parse -- flu-dasm-antigenic-evo's
        `date_filters.parse_iso_date` returns None, dropping the sequence from
        chronumental's node dating without a word. The failure would survive,
        just quieter and one repository further away.
        """
        for bad in ("2019-02-30", "2020-13-45", "2020-00"):
            ids = ["EPI_ISL_1", "EPI_ISL_2"]
            with self.subTest(value=bad), tempfile.TemporaryDirectory() as tmpdir:
                code, _ = run_parse(tmpdir, self.segments_for(ids), ids,
                                    ha_subtypes=["H1"], na_subtypes=["N1"],
                                    dates=["2010-03-01", bad])
                self.assertEqual(code, 1)

    def test_a_leap_day_is_accepted_in_a_leap_year_only(self):
        """Guards the guard: validation strict enough to reject real dates
        would fail the pipeline on ~1/1460 of the corpus."""
        self.assertTrue(parse_gisaid_data.is_gisaid_date("2020-02-29"))
        self.assertFalse(parse_gisaid_data.is_gisaid_date("2019-02-29"))

    def test_the_error_names_the_file_the_bad_date_came_from(self):
        """With 56 input files, an error that does not say which one is nearly
        as hard to act on as no error. Two files, bad date only in the second."""
        good = fake_metadata(["EPI_ISL_1", "EPI_ISL_2"],
                             dates=["2020-01-01", "2021-01-01"])
        bad = fake_metadata(["EPI_ISL_3", "EPI_ISL_4"],
                            dates=["2020-01-01", "03/04/2020"])
        ids = ["EPI_ISL_1", "EPI_ISL_2", "EPI_ISL_3", "EPI_ISL_4"]
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertLogs("parse_gisaid_data", level="ERROR") as logged:
                code, _ = run_parse(tmpdir, self.segments_for(ids), ids,
                                    ha_subtypes=["H1"], na_subtypes=["N1"],
                                    metadata_frames=[good, bad])
            self.assertEqual(code, 1)
            reported = "\n".join(logged.output)
            self.assertIn("meta_01.xls", reported)
            self.assertNotIn("meta_00.xls", reported)
            self.assertIn("EPI_ISL_4", reported)

    def test_whitespace_only_date_counts_as_blank_not_malformed(self):
        """Stripping happens before the blank test, so "  " is missing data."""
        ids = ["EPI_ISL_1", "EPI_ISL_2"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, self.segments_for(ids), ids,
                                  ha_subtypes=["H1"], na_subtypes=["N1"],
                                  dates=["2010-03-01", "   "])
            self.assertEqual(code, 0)
            self.assertEqual(read_collection_dates(out)["EPI_ISL_2"], "")

    def test_an_entirely_blank_column_is_not_a_malformed_column(self):
        """The all-blank case must not trip the fatal check on emptiness."""
        ids = ["EPI_ISL_1", "EPI_ISL_2"]
        with tempfile.TemporaryDirectory() as tmpdir:
            code, out = run_parse(tmpdir, self.segments_for(ids), ids,
                                  ha_subtypes=["H1"], na_subtypes=["N1"],
                                  dates=["", ""])
            self.assertEqual(code, 0)
            self.assertEqual(set(read_collection_dates(out).values()), {""})


if __name__ == "__main__":
    unittest.main()
