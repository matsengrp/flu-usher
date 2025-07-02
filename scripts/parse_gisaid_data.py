"""
Script to parse GISAID data files from multiple directories and organize sequences by segment.
For HA and NA segments, further split by subtype (e.g., H1, N1).
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
import re

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Parse GISAID data files from multiple directories')
    parser.add_argument('--input-dirs', nargs='+', required=True, help='Input directories containing GISAID data')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--segments', nargs='+', help='List of segments to process (e.g., HA NA)')
    return parser.parse_args()

def extract_ha_na_subtype(subtype_str):
    """Extract HA and NA subtypes from a full subtype string like H1N1"""
    # Match patterns like H1N1, H3N2, etc.
    match = re.match(r'(H\d+)(N\d+)', subtype_str)
    if match:
        return match.group(1), match.group(2)
    return None, None

def main():
    # Parse command line arguments
    args = parse_args()
    
    print(f"Processing data from {len(args.input_dirs)} input directories")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # List of segments to keep (from config file or all segments if not specified)
    valid_segments = set(args.segments) if args.segments else None
    
    # Iterate over sequences, group by segment (and subtype for HA/NA)
    # For HA/NA: segment_records[segment][subtype] = [records]
    # For others: segment_records[segment]['all'] = [records]
    segment_records = defaultdict(lambda: defaultdict(list))
    
    # Track which EPI_ISL IDs have been seen for each segment-subtype combination
    segment_subtype_epi_isl_ids = defaultdict(lambda: defaultdict(set))
    
    # Collect all metadata dataframes
    all_metadata_dfs = []
    
    # Process each input directory
    for data_dir in args.input_dirs:
        if not os.path.exists(data_dir):
            print(f"Warning: Directory {data_dir} does not exist, skipping")
            continue
            
        print(f"\nProcessing directory: {data_dir}")
        
        fasta_files = glob.glob(os.path.join(data_dir, "*.fasta"))
        metadata_files = glob.glob(os.path.join(data_dir, "*.xls"))
        
        if len(fasta_files) == 0 or len(metadata_files) == 0:
            print(f"Warning: Expected at least one FASTA file and one XLS file in {data_dir}")
            print(f"Found {len(fasta_files)} FASTA files and {len(metadata_files)} XLS files")
            continue
        
        # Process FASTA files
        for fasta_file in fasta_files:
            records = list(SeqIO.parse(fasta_file, "fasta"))
            print(f"  Read {len(records)} records from {fasta_file}")
            
            for record in records:
                try:
                    # Parse the sequence ID
                    (epi, segment, name, epi_isl, seq_subtype) = record.id.split('|')
                    
                    # Skip segments not in the config file
                    if valid_segments and segment not in valid_segments:
                        continue
                    
                    # Determine the grouping key based on segment type
                    if segment == 'HA':
                        # Extract H subtype (e.g., H1 from H1N1)
                        ha_subtype, _ = extract_ha_na_subtype(seq_subtype)
                        if ha_subtype:
                            group_key = ha_subtype
                        else:
                            print(f"Warning: Could not extract HA subtype from {seq_subtype}")
                            continue
                    elif segment == 'NA':
                        # Extract N subtype (e.g., N1 from H1N1)
                        _, na_subtype = extract_ha_na_subtype(seq_subtype)
                        if na_subtype:
                            group_key = na_subtype
                        else:
                            print(f"Warning: Could not extract NA subtype from {seq_subtype}")
                            continue
                    else:
                        # For non-HA/NA segments, group all together
                        group_key = 'all'
                    
                    # Skip record if we've already seen this EPI_ISL ID for this segment-group combination
                    if epi_isl in segment_subtype_epi_isl_ids[segment][group_key]:
                        print(f"Warning: Duplicate EPI_ISL ID {epi_isl} found for segment {segment} group {group_key}. Skipping.")
                        continue
                    
                    segment_subtype_epi_isl_ids[segment][group_key].add(epi_isl)
                    
                    # Update record ID and description to just use the EPI_ISL
                    record.id = epi_isl
                    record.description = epi_isl
                    
                    # Store record by segment and group
                    segment_records[segment][group_key].append(record)
                    
                except ValueError:
                    print(f"Warning: Could not parse ID for record: {record.id}")
                    continue
        
        # Process metadata files
        for metadata_file in metadata_files:
            df = pd.read_excel(metadata_file, sheet_name=0)
            print(f"  Read {len(df)} metadata records from {metadata_file}")
            
            # Subset to columns of interest
            cols = [
                'Isolate_Id', 'Isolate_Name', 'Subtype', 'Clade', 'Passage_History',
                'Location', 'Host', 'Collection_Date'
            ]
            df = (
                df[cols]
                .rename(columns={col : col.lower() for col in cols})
            )
            # Convert the 'collection_date' column to datetime
            df['collection_date'] = pd.to_datetime(df['collection_date'], errors='coerce')
            all_metadata_dfs.append(df)
    
    # Concatenate all metadata dataframes
    if all_metadata_dfs:
        metadata_df = pd.concat(all_metadata_dfs, ignore_index=True)
        
        # Remove duplicate isolate_ids (keep first occurrence)
        if metadata_df.duplicated(subset=['isolate_id']).sum() > 0:
            print(f"Warning: Found {metadata_df.duplicated(subset=['isolate_id']).sum()} duplicate isolate_ids in metadata, keeping first occurrence")
            metadata_df = metadata_df.drop_duplicates(subset=['isolate_id'], keep='first')
        
        # Write combined metadata to CSV
        metadata_output_file = os.path.join(args.output_dir, "combined_metadata.csv")
        print(f"\nWriting combined metadata to {metadata_output_file}")
        metadata_df.to_csv(metadata_output_file, index=False)
    else:
        print("Warning: No metadata files were processed")
    
    # Write sequences to output files organized by segment and subtype
    print("\nSummary of records by segment and subtype:")
    for segment, groups in segment_records.items():
        for group_key, records in groups.items():
            # Create output directory path as segment/group
            segment_output_dir = os.path.join(args.output_dir, segment, group_key)
            if not os.path.isdir(segment_output_dir):
                os.makedirs(segment_output_dir)
            
            output_file = os.path.join(segment_output_dir, "raw_sequences.fasta.xz")
            with lzma.open(output_file, 'wt') as handle:
                SeqIO.write(records, handle, "fasta")
            
            print(f"  {segment}/{group_key}: {len(records)} records written to {output_file}")

if __name__ == "__main__":
    sys.exit(main())