"""
Script to parse GISAID data files from multiple directories and organize sequences by segment.
For HA and NA segments, further split by subtype (e.g., H1, N1).
"""

import argparse
import os
import sys
import glob
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import re
from datetime import date
from utils import open_sequence_file, setup_logging

# GISAID ships collection dates at exactly three precisions.
GISAID_DATE_PATTERN = re.compile(r'\d{4}(-\d{2}(-\d{2})?)?')


def is_gisaid_date(value):
    """True if `value` is a real date at one of GISAID's three precisions.

    Shape alone is not enough: "2020-13-45" and "2019-02-30" match the pattern
    and name no day that exists. Letting them through would not make them
    harmless, only quiet -- flu-dasm-antigenic-evo's `date_filters.parse_iso_date`
    parses strict %Y-%m-%d and returns None on failure, so an impossible date
    would drop out of chronumental's node dating and the root-to-tip regression
    with no error, one repository downstream of here. Missing precision is
    normal; an impossible date is a transcription error.
    """
    if not GISAID_DATE_PATTERN.fullmatch(value):
        return False
    parts = value.split('-')
    try:
        # Absent fields default to the 1st, which is only ever used to ask
        # whether the fields that *are* present name a real day.
        date(int(parts[0]),
             int(parts[1]) if len(parts) > 1 else 1,
             int(parts[2]) if len(parts) > 2 else 1)
    except ValueError:
        return False
    return True

logger = setup_logging(__name__)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Parse GISAID data files from multiple directories')
    parser.add_argument('--input-dirs', nargs='+', required=True, help='Input directories containing GISAID data')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--segments', nargs='+', help='List of segments to process (e.g., HA NA)')
    parser.add_argument('--ha-subtypes', nargs='+', required=True,
                        help='HA subtypes to write (e.g. H1 H3). Records of any other HA subtype are dropped.')
    parser.add_argument('--na-subtypes', nargs='+', required=True,
                        help='NA subtypes to write (e.g. N1 N2). Records of any other NA subtype are dropped.')
    return parser.parse_args()

def extract_ha_na_subtype(subtype_str):
    """Extract HA and NA subtypes from a full subtype string like H1N1 or A_/_H1N1"""
    # First, try to extract the H*N* part from patterns like A_/_H1N1
    if '_/_' in subtype_str:
        # Split by '_/_' and take the part after it
        parts = subtype_str.split('_/_')
        if len(parts) >= 2:
            subtype_str = parts[1]
    
    # Match patterns like H1N1, H3N2, etc.
    match = re.match(r'(H\d+)(N\d+)', subtype_str)
    if match:
        return match.group(1), match.group(2)
    return None, None

def main():
    # Parse command line arguments
    args = parse_args()
    
    logger.info(f"Processing data from {len(args.input_dirs)} input directories")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # List of segments to keep (from config file or all segments if not specified)
    valid_segments = set(args.segments) if args.segments else None

    # Only configured subtypes are written. Without this the script produced 30
    # segment/subtype directories against 16 configured combinations, and the 14
    # extras were consumed by nothing.
    valid_ha_subtypes = set(args.ha_subtypes)
    valid_na_subtypes = set(args.na_subtypes)
    expected_combinations = (
        {('HA', st) for st in valid_ha_subtypes}
        | {('NA', st) for st in valid_na_subtypes}
        | {(seg, 'all') for seg in (valid_segments or set()) if seg not in ('HA', 'NA')}
    )

    # Every record that does not reach an output file is counted under exactly one
    # reason, so the totals below reconcile against the number read.
    skipped = defaultdict(int)
    total_records_read = 0
    
    # Iterate over sequences, group by segment (and subtype for HA/NA)
    # For HA/NA: segment_records[segment][subtype] = [records]
    # For others: segment_records[segment]['all'] = [records]
    segment_records = defaultdict(lambda: defaultdict(list))
    
    # Track which EPI_ISL IDs have been seen for each segment-subtype combination
    segment_subtype_epi_isl_ids = defaultdict(lambda: defaultdict(set))
    
    # Collect all metadata dataframes
    all_metadata_dfs = []

    # Collection dates that are neither blank nor one of GISAID's three
    # precisions. Reported per file with example isolates, and fatal, because the
    # only way to get one is a change in what GISAID ships.
    malformed_dates = []
    
    # Process each input directory
    for data_dir in args.input_dirs:
        if not os.path.exists(data_dir):
            logger.warning(f"Directory {data_dir} does not exist, skipping")
            continue
            
        logger.info(f"Processing directory: {data_dir}")
        
        # sorted(): glob returns readdir order, which is arbitrary. Record order
        # in raw_sequences.fasta.xz feeds randomize_alignment --seed {n} and
        # decides which duplicate survives the keep-first dedup below, so an
        # unsorted glob makes the trees depend on filesystem layout.
        # Snakefile's INPUT_DATA_FILES already applies this discipline.
        fasta_files = sorted(glob.glob(os.path.join(data_dir, "*.fasta")))
        metadata_files = sorted(glob.glob(os.path.join(data_dir, "*.xls")))
        
        if len(fasta_files) == 0 or len(metadata_files) == 0:
            logger.warning(f"Expected at least one FASTA file and one XLS file in {data_dir}")
            logger.warning(f"Found {len(fasta_files)} FASTA files and {len(metadata_files)} XLS files")
            continue
        
        # Process FASTA files
        for fasta_file in fasta_files:
            records = list(SeqIO.parse(fasta_file, "fasta"))
            total_records_read += len(records)
            logger.info(f"Read {len(records)} records from {fasta_file}")
            
            for record in records:
                try:
                    # Parse the sequence ID
                    (epi, segment, name, epi_isl, seq_subtype) = record.id.split('|')
                    
                    # Skip segments not in the config file
                    if valid_segments and segment not in valid_segments:
                        skipped['segment_not_configured'] += 1
                        continue
                    
                    # Determine the grouping key based on segment type
                    if segment == 'HA':
                        # Extract H subtype (e.g., H1 from H1N1)
                        ha_subtype, _ = extract_ha_na_subtype(seq_subtype)
                        if not ha_subtype:
                            skipped['ha_subtype_unparseable'] += 1
                            continue
                        if ha_subtype not in valid_ha_subtypes:
                            skipped['ha_subtype_not_configured'] += 1
                            continue
                        group_key = ha_subtype
                    elif segment == 'NA':
                        # Extract N subtype (e.g., N1 from H1N1)
                        _, na_subtype = extract_ha_na_subtype(seq_subtype)
                        if not na_subtype:
                            skipped['na_subtype_unparseable'] += 1
                            continue
                        if na_subtype not in valid_na_subtypes:
                            skipped['na_subtype_not_configured'] += 1
                            continue
                        group_key = na_subtype
                    else:
                        # For non-HA/NA segments, group all together
                        group_key = 'all'
                    
                    # Skip record if we've already seen this EPI_ISL ID for this segment-group combination
                    if epi_isl in segment_subtype_epi_isl_ids[segment][group_key]:
                        # Legitimate and expected: the download windows in
                        # data/H3N2/ and data/H1N1/ overlap by design. Reported in
                        # aggregate rather than one warning per record.
                        skipped['duplicate_epi_isl'] += 1
                        continue
                    
                    segment_subtype_epi_isl_ids[segment][group_key].add(epi_isl)
                    
                    # Update record ID and description to just use the EPI_ISL
                    record.id = epi_isl
                    record.description = epi_isl
                    
                    # Store record by segment and group
                    segment_records[segment][group_key].append(record)
                    
                except ValueError:
                    skipped['unparseable_sequence_id'] += 1
                    continue
        
        # Process metadata files
        for metadata_file in metadata_files:
            df = pd.read_excel(metadata_file, sheet_name=0)
            logger.info(f"Read {len(df)} metadata records from {metadata_file}")
            
            # Subset to columns of interest
            cols = [
                'Isolate_Id', 'Isolate_Name', 'Subtype', 'Clade', 'Passage_History',
                'Location', 'Host', 'Collection_Date'
            ]
            df = (
                df[cols]
                .rename(columns={col : col.lower() for col in cols})
            )
            # Keep collection_date exactly as GISAID supplies it, at whichever
            # of its three precisions that is (YYYY-MM-DD, YYYY-MM, YYYY). It
            # used to go through a bare pd.to_datetime(errors='coerce'), which
            # silently emptied 87,854 of 629,106 dates -- see "Collection dates"
            # in the Design Notes of README.md for why, and issue #55.
            raw = df['collection_date'].astype('string').str.strip()
            blank = raw.isna() | (raw == '')
            # fillna('') so the check is plain bool rather than nullable: a
            # blank is judged by `blank`, not here.
            malformed = ~blank & ~raw.fillna('').map(is_gisaid_date)
            if malformed.any():
                examples = ', '.join(df.loc[malformed, 'isolate_id'].astype(str).head(5))
                malformed_dates.append(
                    f"{metadata_file}: {int(malformed.sum())} collection_date "
                    f"values are not a real date at one of YYYY-MM-DD, "
                    f"YYYY-MM, YYYY "
                    f"(e.g. {examples})"
                )
            if blank.any():
                logger.info(
                    f"{int(blank.sum())} of {len(df)} records in {metadata_file} "
                    f"have no collection_date"
                )
            df['collection_date'] = raw
            all_metadata_dfs.append(df)
    
    # Concatenate all metadata dataframes
    if all_metadata_dfs:
        metadata_df = pd.concat(all_metadata_dfs, ignore_index=True)
        
        # Remove duplicate isolate_ids (keep first occurrence)
        if metadata_df.duplicated(subset=['isolate_id']).sum() > 0:
            logger.warning(f"Found {metadata_df.duplicated(subset=['isolate_id']).sum()} duplicate isolate_ids in metadata, keeping first occurrence")
            metadata_df = metadata_df.drop_duplicates(subset=['isolate_id'], keep='first')
        
        # Write combined metadata to CSV
        metadata_output_file = os.path.join(args.output_dir, "combined_metadata.csv")
        logger.info(f"Writing combined metadata to {metadata_output_file}")
        metadata_df.to_csv(metadata_output_file, index=False)
        metadata_written = True
    else:
        logger.error("No metadata files were processed")
        metadata_written = False
    
    # Write sequences to output files organized by segment and subtype
    logger.info("Summary of records by segment and subtype:")
    written_counts = {}
    for segment, groups in segment_records.items():
        for group_key, records in groups.items():
            # Create output directory path as segment/group
            segment_output_dir = os.path.join(args.output_dir, segment, group_key)
            if not os.path.isdir(segment_output_dir):
                os.makedirs(segment_output_dir)

            output_file = os.path.join(segment_output_dir, "raw_sequences.fasta.xz")
            with open_sequence_file(output_file, 'wt') as handle:
                SeqIO.write(records, handle, "fasta")

            written_counts[(segment, group_key)] = len(records)
            logger.info(f"{segment}/{group_key}: {len(records)} records written to {output_file}")

    # Aggregate accounting. Every record read is either written or counted under
    # exactly one skip reason, so these totals reconcile.
    total_written = sum(written_counts.values())
    total_skipped = sum(skipped.values())
    logger.info("")
    logger.info("Record accounting:")
    logger.info(f"  read:    {total_records_read:,}")
    logger.info(f"  written: {total_written:,}")
    logger.info(f"  skipped: {total_skipped:,}")
    for reason in sorted(skipped):
        logger.info(f"    - {reason.replace('_', ' ')}: {skipped[reason]:,}")
    # This cannot fire against the loop as currently written -- every `continue`
    # increments exactly one `skipped` counter, and the only other path appends
    # to segment_records -- so it is a guard against future drift, not a runtime
    # condition. Kept because the failure it catches (a new `continue` added
    # without a matching counter) is silent otherwise, and fatal rather than a
    # warning so it matches the exit-nonzero posture of the checks below.
    reconciliation_error = None
    if total_records_read != total_written + total_skipped:
        reconciliation_error = (
            f"accounting does not reconcile: {total_records_read:,} read but "
            f"{total_written:,} written + {total_skipped:,} skipped"
        )

    # Fail loudly. Previously main() had no return, so `sys.exit(main())` was
    # always `sys.exit(None)` -- exit 0 even if every record had been dropped or
    # no metadata written, and Snakemake would record the outputs as good.
    errors = []
    if not metadata_written:
        errors.append("no metadata was written")
    missing = sorted(expected_combinations - set(written_counts))
    if missing:
        errors.append(
            "configured combinations produced no sequences: "
            + ", ".join(f"{seg}/{sub}" for seg, sub in missing)
        )
    # A "wrote zero records" check used to live here. It was unreachable:
    # written_counts is populated only from segment_records, whose entries are
    # created by the .append() above, so no entry can have length zero. A
    # combination that produced nothing is absent entirely, which `missing`
    # already catches.
    if reconciliation_error:
        errors.append(reconciliation_error)
    errors.extend(malformed_dates)
    if errors:
        for err in errors:
            logger.error(err)
        return 1

    logger.info(
        f"Wrote {len(written_counts)} segment/subtype combinations "
        f"({len(expected_combinations)} configured)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())