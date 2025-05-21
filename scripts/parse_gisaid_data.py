"""
Script to parse GISAID data files and organize sequences by segment
"""

import argparse
import os
import sys
import glob
import lzma
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Parse GISAID data files')
    parser.add_argument('--subtype', required=True, help='Flu subtype (e.g., H7N9)')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--segments', nargs='+', help='List of segments to process (e.g., HA NA)')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    print(f"Processing data for {args.subtype}")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Read in sequence data and metadata
    data_dir = os.path.join('data', args.subtype)
    fasta_files = glob.glob(os.path.join(data_dir, "*.fasta"))
    metadata_files = glob.glob(os.path.join(data_dir, "*.xls"))
    
    # Verify that we have at least one fasta and one metadata file
    if len(fasta_files) == 0 or len(metadata_files) == 0:
        print(f"Error: Expected at least one FASTA file and one XLS file in {data_dir}")
        print(f"Found {len(fasta_files)} FASTA files and {len(metadata_files)} XLS files")
        return 1
    
    # Iterate over sequences, group by segment, and make one FASTA file per
    # segment with all sequences for that segment
    segment_records = defaultdict(list)
    epi_isl_ids = set()
    
    # Track which EPI_ISL IDs have been seen for each segment
    segment_epi_isl_ids = defaultdict(set)
    
    # List of segments to keep (from config file or all segments if not specified)
    valid_segments = set(args.segments) if args.segments else None
    
    for fasta_file in fasta_files:
        # Read in records
        records = list(SeqIO.parse(fasta_file, "fasta"))
        print(f"Read {len(records)} records from {fasta_file}")
        
        # Iterate through records and parse metadata from the sequence ID
        for record in records:
            # Parse the sequence ID
            try:
                (epi, segment, name, epi_isl, seq_subtype) = record.id.split('|')
                
                # Update set of EPI_ISL IDs
                epi_isl_ids.add(epi_isl)

                # Skip segments not in the config file
                if valid_segments and segment not in valid_segments:
                    continue
                
                # Skip record if we've already seen this EPI_ISL ID for this segment
                if epi_isl in segment_epi_isl_ids[segment]:
                    print(f"Warning: Duplicate EPI_ISL ID {epi_isl} found for segment {segment}. Skipping record.")
                    continue
                segment_epi_isl_ids[segment].add(epi_isl)
                
                # Check that the subtype matches what we expect
                if not seq_subtype.endswith(args.subtype):
                    print(f"Warning: Record {record.id} has mismatched subtype: {seq_subtype} vs expected {args.subtype}")
                    continue
                
                # Update record ID and description to just use the EPI_ISL
                record.id = epi_isl
                record.description = epi_isl
                
                # Store record by segment
                segment_records[segment].append(record)
                
            except ValueError:
                print(f"Warning: Could not parse ID for record: {record.id}")
                continue
    
    # For each segment, print summary and then write the records to a compressed XZ FASTA file
    print(f"\nTotal unique EPI_ISL IDs: {len(epi_isl_ids)}")
    print("\nSummary of records by segment:")
    for segment, records in segment_records.items():
        assert len(records) == len(segment_epi_isl_ids[segment]), "Mismatch between number of records and EPI_ISL IDs"
        segment_output_dir = os.path.join(args.output_dir, segment)
        if not os.path.isdir(segment_output_dir):
            os.makedirs(segment_output_dir)
        output_file = os.path.join(segment_output_dir, "raw_sequences.fasta.xz")
        with lzma.open(output_file, 'wt') as handle:
            SeqIO.write(records, handle, "fasta")
        print(f"Segment {segment}: {len(records)} records written to {output_file}") 
    
    # Load all metadata and write to an output CSV file
    dfs = []
    for metadata_file in metadata_files:
        
        # Read in data from the first sheet
        df = pd.read_excel(metadata_file, sheet_name=0)
        print(f"Read {len(df)} records from {metadata_file}")
        
        # Subset to a few columns of interest
        cols = [
            'Isolate_Id', 'Isolate_Name', 'Subtype', 'Clade', 'Passage_History',
            'Location', 'Host', 'Collection_Date'
        ]
        df = (
            df[cols]
            .rename(columns={col : col.lower() for col in cols})
        )
        dfs.append(df)
    
    # Concatenate all metadata dataframes
    metadata_df = pd.concat(dfs)
    
    # Make sure that all Isolate_Id values are unique
    assert sum(metadata_df.duplicated(subset=['isolate_id'])) == 0, "Duplicate isolate_id values found in metadata."

    # Make sure that all observed epi_isl_ids are in the metadata
    metadata_epi_isl_ids = set(metadata_df['isolate_id'].unique())
    missing_ids = epi_isl_ids.difference(metadata_epi_isl_ids)
    if missing_ids:
        print(f"{len(missing_ids)} EPI_ISL IDs from sequences are not found in metadata")
        print(f"First few missing IDs: {list(missing_ids)[:5]}")
        return 1
    
    # Write metadata to CSV
    metadata_output_file = os.path.join(args.output_dir, f"{args.subtype}_metadata.csv")
    print(f"\nWriting metadata to {metadata_output_file}")
    metadata_df.to_csv(metadata_output_file, index=False)

if __name__ == "__main__":
    sys.exit(main())

