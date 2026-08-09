"""
Tests for parse_gisaid_data.py.

Two behaviours are under test:

1. Only configured subtypes are written. The script previously produced 30
   segment/subtype directories against 16 configured combinations.
2. main() returns a nonzero exit code when a configured combination ends up
   empty or no metadata is written. It previously had no return statement at
   all, so `sys.exit(main())` was always `sys.exit(None)` -- exit 0 no matter
   how much was dropped.

pandas.read_excel is patched rather than writing real .xls fixtures; the Excel
parsing itself is not what these tests are about.
"""

import lzma
import os
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


def fake_metadata(isolate_ids):
    return pd.DataFrame([
        {
            "Isolate_Id": i, "Isolate_Name": f"A/test/{i}/2020",
            "Subtype": "A / H1N1", "Clade": "x", "Passage_History": "Original",
            "Location": "North America / USA", "Host": "Human",
            "Collection_Date": "2020-01-01",
        }
        for i in isolate_ids
    ], columns=METADATA_COLUMNS)


def run_parse(tmpdir, records, isolate_ids, ha_subtypes, na_subtypes,
              segments=("HA", "NA", "PB2"), with_metadata=True):
    """Write inputs, run main(), return (exit_code, output_dir)."""
    data_dir = os.path.join(tmpdir, "data")
    os.makedirs(data_dir, exist_ok=True)
    with open(os.path.join(data_dir, "seqs.fasta"), "w") as f:
        for rid in records:
            f.write(f">{rid}\nACGTACGTACGT\n")
    if with_metadata:
        # Content is supplied by the read_excel patch; the file must merely exist
        # for the script's glob to find it.
        open(os.path.join(data_dir, "meta.xls"), "w").close()

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
    with mock.patch.object(sys, "argv", argv), \
         mock.patch.object(pd, "read_excel", return_value=fake_metadata(isolate_ids)):
        code = parse_gisaid_data.main()
    return code, out_dir


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


if __name__ == "__main__":
    unittest.main()
