# Flu-UShER Pipeline

A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER. The pipeline processes influenza sequences by segment and subtype, creating separate trees for HA and NA segments by subtype (e.g., H1, H3, N1, N2) while combining all subtypes for other segments (PB2, PB1, PA, NP, MP, NS).

## Directory Structure

```
flu-usher/
├── Snakefile                # Main pipeline file
├── config/
│   └── config.yaml          # Configuration file
├── data/                    # Input data directories
│   ├── H1N1/                # Example: H1N1 sequences (all segments)
│   ├── H3N2/                # Example: H3N2 sequences (all segments)
│   ├── H5N1/                # Example: H5N1 sequences (all segments)
│   └── H7N9/                # Example: H7N9 sequences (all segments)
│       ├── sequences.fasta  # One or more FASTA files containing sequences
│       └── metadata.xls     # One or more metadata files in Excel format
├── logs/                    # Log files (created by the pipeline)
├── results/                 # Output results (organized by segment/subtype)
│   ├── HA/                  # HA segment results by subtype
│   │   ├── H1/              # H1 subtype tree and files
│   │   ├── H3/              # H3 subtype tree and files
│   │   └── ...
│   ├── NA/                  # NA segment results by subtype
│   │   ├── N1/              # N1 subtype tree and files
│   │   ├── N2/              # N2 subtype tree and files
│   │   └── ...
│   └── PB2/                 # Other segment results
│       └── all/             # All subtypes combined
├── scripts/
│   ├── parse_gisaid_data.py # Script to parse GISAID data by segment
│   ├── curate_msa.py        # Sequence curation script
│   └── download_ref_seq.py  # Reference sequence download script
└── notebooks/               # Jupyter notebooks for development and analysis
```

## Usage

1. **Set up your environment**

   ```bash
   # Create the conda environment
   conda env create -f environment.yml
   
   # Activate the environment
   conda activate usher
   ```

2. **Configure the pipeline**

   Edit `config/config.yaml` to:
   - Specify input directories containing your GISAID data
   - List HA subtypes to analyze (e.g., H1, H3, H5, H7, H9)
   - List NA subtypes to analyze (e.g., N1, N2, N9)
   - Set reference accession numbers for each segment-subtype combination
   - Adjust filtering thresholds for sequence curation (max gaps and ambiguities)
   - Set desired number of threads

3. **Prepare your GISAID data**

   The pipeline expects two types of files with input data downloaded from GISAID:
   - One or more FASTA files containing nucleotide sequences
      - FASTA sequence ID format: `EPI|SEGMENT|NAME|EPI_ISL|SUBTYPE`
   - One or more XLS files with metadata for the above sequences
   
   The pipeline will automatically:
   - Aggregate all FASTA and metadata files from each input directory
   - Parse sequences by segment
   - Group HA and NA segments by subtype (defined by the `subtype` column in the metadata file)
   - Combine all subtypes for internal segments (PB2, PB1, PA, NP, MP, NS)

4. **Run the pipeline**

   ```bash
   # Test run (dry-run)
   snakemake -np
   
   # Run the full pipeline
   snakemake --cores <number_of_cores>
   
   # Run for specific segment-subtype combinations
   snakemake --cores 8 results/HA/H5/opt_tree.pb.gz
   snakemake --cores 8 results/NA/N1/opt_tree.pb.gz
   snakemake --cores 8 results/PB2/all/opt_tree.pb.gz
   ```

5. **Output**

   The pipeline organizes outputs by segment and subtype:
   
   **For HA and NA segments** (e.g., `results/HA/H5/` or `results/NA/N1/`):
   - `raw_sequences.fasta.xz`: Parsed sequences for this segment-subtype
   - `reference/`: Reference sequence and Nextclade dataset files
   - `msa.fasta.xz`: Multiple sequence alignment from Nextclade
   - `curated_msa.fasta.xz`: Quality-filtered alignment (e.g., gaps < 5%, ambiguities < 1%)
   - `curated_reference.fasta/gff/gtf`: Reference files matching the curated alignment
   - `curated_msa.vcf.gz`: Variant call format file for UShER
   - `preopt_tree.pb.gz`: Initial parsimony tree
   - `opt_tree.pb.gz`: Optimized phylogenetic tree
   - `opt_tree.jsonl.gz`: Interactive Taxonium visualization file
   
   **For the other segments** (e.g., `results/PB2/all/` or `results/NP/all/`):
   - Same outputs as above, but combining all influenza subtypes
   
   **Global outputs**:
   - `results/combined_metadata.csv`: Aggregated metadata from all input files

## Pipeline Steps

1. **Parse GISAID data** (`parse_gisaid_data.py`): 
   - Aggregates sequences from multiple input directories
   - Splits by segment and subtype
   - Creates unified metadata file
   
2. **Download reference** (`download_ref_seq.py`): 
   - Fetches appropriate reference sequence for each segment-subtype
   - Creates Nextclade dataset configuration
   
3. **Align sequences** (Nextclade): 
   - Performs codon-aware alignment to reference
   - Includes reference in the output alignment
   
4. **Curate alignment** (`curate_msa.py`): 
   - Filters sequences by quality (gaps, ambiguities)
   - Extracts only coding regions based on GFF annotations
   - Sanitizes sequence IDs for UShER compatibility
   
5. **Create VCF** (faToVcf): 
   - Converts FASTA alignment to variant format
   - Includes reference and handles ambiguous nucleotides
   
6. **Build initial tree** (usher-sampled): 
   - Creates parsimony-based phylogenetic tree
   
7. **Optimize tree** (matOptimize): 
   - Refines tree topology to minimize parsimony score
   
8. **Create visualization** (usher_to_taxonium): 
   - Converts tree to Taxonium format
   - Incorporates metadata for interactive exploration

## Requirements

All dependencies are specified in the environment.yml file.
