from Bio import Entrez, SeqIO
import os
import json
import argparse
import sys
import time
import yaml

def parse_args():
    parser = argparse.ArgumentParser(description="Download reference sequences for flu-usher")
    parser.add_argument("--config", required=True, help="Path to config.yaml file")
    parser.add_argument("--output-base-dir", required=True, help="Base output directory (e.g., 'results')")
    parser.add_argument("--wait-time", type=int, default=30, help="Seconds to wait between downloads (default: 30)")
    return parser.parse_args()

def download_gene_sequence(accession, output_dir):
    """Download gene sequence in FASTA format using accession number"""
    # Fetch the sequence
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    sequence = handle.read()
    handle.close()

    # Save to file
    output_fasta = os.path.join(output_dir, "reference.fasta")
    with open(output_fasta, "w") as f:
        f.write(sequence)

    print(f"  Sequence saved to {output_fasta}")
    return output_fasta

def download_gene_gff(accession, output_dir):
    """Download GFF file for gene using accession number"""
    # First get the GI number or other identifiers if needed
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Fetch the GFF
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gff3", retmode="text")
    gff_content = handle.read()
    handle.close()

    # Save to file
    output_gff = os.path.join(output_dir, "reference.gff")
    with open(output_gff, "w") as f:
        f.write(gff_content)

    print(f"  GFF saved to {output_gff}")
    return output_gff

def create_pathogen_json(output_dir, fasta_file, gff_file):
    """Create pathogen.json for Nextclade"""
    pathogen_json = {
        "schemaVersion": "3.0.0",
        "files": {
            "genomeAnnotation": os.path.basename(gff_file),
            "pathogenJson": "pathogen.json",
            "reference": os.path.basename(fasta_file),
        },
        "alignmentParams": {
            "minSeedCover": 0.1
        }
    }

    # Write pathogen_json to a JSON file
    json_file_path = os.path.join(output_dir, "pathogen.json")
    with open(json_file_path, 'w') as f:
        json.dump(pathogen_json, f, indent=2)

    print(f"  Pathogen JSON saved to {json_file_path}")
    return json_file_path

def download_reference_set(segment, subtype, accession, output_base_dir):
    """Download FASTA, GFF, and create pathogen.json for a segment/subtype combination"""
    # Create output directory
    output_dir = os.path.join(output_base_dir, segment, subtype, "reference")

    print(f"\nProcessing {segment}/{subtype} (accession: {accession})")

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        print(f"  Created directory: {output_dir}")
    else:
        print(f"  Directory exists: {output_dir}")

    try:
        # Download FASTA
        fasta_file = download_gene_sequence(accession, output_dir)

        # Wait 3 seconds before downloading GFF
        time.sleep(3)

        # Download GFF
        gff_file = download_gene_gff(accession, output_dir)

        # Create pathogen.json
        json_file = create_pathogen_json(output_dir, fasta_file, gff_file)

        print(f"  ✓ Successfully downloaded reference files for {segment}/{subtype}")
        return True

    except Exception as e:
        print(f"  ✗ Error downloading reference for {segment}/{subtype}: {e}", file=sys.stderr)
        return False

def main():
    args = parse_args()

    # Set Entrez email
    Entrez.email = 'test_window238476@gmail.com'

    # Load config file
    print(f"Loading configuration from {args.config}")
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # Build list of segment/subtype/accession combinations
    combinations = []

    # Add HA subtypes
    for subtype in config["ha_subtypes"]:
        ref_key = f"HA_{subtype}"
        if ref_key in config["references"]:
            combinations.append(("HA", subtype, config["references"][ref_key]))

    # Add NA subtypes
    for subtype in config["na_subtypes"]:
        ref_key = f"NA_{subtype}"
        if ref_key in config["references"]:
            combinations.append(("NA", subtype, config["references"][ref_key]))

    # Add internal segments (all subtypes combined)
    for segment in config["segments"]:
        if segment not in ["HA", "NA"]:
            ref_key = f"{segment}_all"
            if ref_key in config["references"]:
                combinations.append((segment, "all", config["references"][ref_key]))

    print(f"\nFound {len(combinations)} segment/subtype combinations to process")
    print(f"Wait time between downloads: {args.wait_time} seconds")
    print(f"Estimated total runtime: ~{len(combinations) * args.wait_time // 60} minutes\n")
    print("=" * 70)

    # Process each combination sequentially
    success_count = 0
    for i, (segment, subtype, accession) in enumerate(combinations, 1):
        # Wait before downloading (except for first download)
        if i > 1:
            print(f"\nWaiting {args.wait_time} seconds before next download...")
            time.sleep(args.wait_time)

        print(f"\n[{i}/{len(combinations)}]", end=" ")

        # Download reference files
        success = download_reference_set(segment, subtype, accession, args.output_base_dir)

        if success:
            success_count += 1
        else:
            print(f"\nFailed to download reference for {segment}/{subtype}", file=sys.stderr)
            sys.exit(1)

    # Final summary
    print("\n" + "=" * 70)
    print(f"\nCompleted: {success_count}/{len(combinations)} downloads successful")
    print(f"All reference data downloaded successfully to {args.output_base_dir}\n")

if __name__ == "__main__":
    main()
