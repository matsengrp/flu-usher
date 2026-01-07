from Bio import Entrez, SeqIO
import os
import json
import argparse
import sys
import time
import random

def parse_args():
    parser = argparse.ArgumentParser(description="Download reference sequences for flu-usher")
    parser.add_argument("--accession", required=True, help="Reference sequence accession number")
    parser.add_argument("--output-dir", required=True, help="Output directory for reference files")
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
    
    print(f"Sequence saved to {output_fasta}")
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
    
    print(f"GFF saved to {output_gff}")
    return output_gff

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    print("Current working directory:", os.getcwd())
    print("Output directory:", args.output_dir)
    if os.path.isdir(args.output_dir):
        print("Output directory already exists:", args.output_dir)
    if not os.path.isdir(args.output_dir):
        print("Making output directory:", args.output_dir)
        os.makedirs(args.output_dir)
    
    # Download reference files
    Entrez.email = 'test_window238476@gmail.com'

    # Wait for a random amount of time before downloading
    wait_time = random.uniform(0, 500)
    print(f"Waiting {wait_time:.1f} seconds before downloading to avoid error related to too many requests.")
    time.sleep(wait_time)
    fasta_file = download_gene_sequence(args.accession, args.output_dir)
    time.sleep(20)  # Wait for 20 seconds before downloading GFF
    gff_file = download_gene_gff(args.accession, args.output_dir)
    
    # Create pathogen.json for Nextclade
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
    json_file_path = os.path.join(args.output_dir, "pathogen.json")
    with open(json_file_path, 'w') as f:
        json.dump(pathogen_json, f, indent=2)
    
    print(f"Pathogen JSON saved to {json_file_path}")
    print(f"Reference data downloaded successfully to {args.output_dir}")

if __name__ == "__main__":
    main()