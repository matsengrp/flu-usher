# Flu-UShER Pipeline

A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER.

## TODO

* troubleshoot H7N9 NA
* H3N2 HA multiple files
* download all H3N2 data
* summary plots for alignments: number downloaded, length distribution, number retained
* only keep a few metadata columns?
* add metadata to taxonium tree
* rerooting (genome with sequences for all segments)

## Directory Structure

```
flu-usher/
├── Snakefile                # Main pipeline file
├── config/
│   └── config.yaml          # Configuration file
├── data/                    # Input data directory (organized by subtype)
│   └── H7N9/                # Example input directory for H7N9 subtype
│       ├── sequences.fasta  # One or more FASTA files containing sequences
│       └── metadata.xls     # One or more metadata files in Excel format
├── logs/                    # Log files (created by the pipeline)
├── results/                 # Output results (organized by subtype-segment)
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
   - List your subtypes and segments to analyze
   - Set reference accession numbers for each subtype-segment combination
   - Adjust filtering thresholds for sequence curation
   - Set desired number of threads

3. **Prepare your GISAID data**

   For each influenza subtype you want to analyze (e.g., H7N9):
   - Create a directory under `data/` with the name of the subtype (e.g., `data/H7N9/`)
   - Place at least one FASTA file containing sequences in the directory (e.g., `data/H7N9/sequences.fasta`)
   - Place at least one Excel (.xls) metadata file in the same directory (e.g., `data/H7N9/metadata.xls`)
   
   The pipeline will aggregate data from all FASTA and metadata files present in each subtype directory.
   
   The FASTA file should have sequence IDs in the format: `EPI|SEGMENT|NAME|EPI_ISL|SUBTYPE`

4. **Run the pipeline**

   ```bash
   # Test run (dry-run)
   snakemake -np
   
   # Run the pipeline
   snakemake --cores <number_of_cores>
   
   # Run for specific subtype-segment combinations
   snakemake --cores <number_of_cores> results/H7N9-HA/opt_tree.pb.gz
   ```

5. **Output**

   For each subtype/segment combination, the pipeline produces:
   - `results/<subtype>/<segment>/raw_sequences.fasta`: Parsed sequences for this segment
   - `results/<subtype>/<segment>/reference/`: Reference data for Nextclade
   - `results/<subtype>/<segment>/msa.fasta.xz`: Aligned sequences
   - `results/<subtype>/<segment>/curated_msa.fasta.xz`: Curated alignment
   - `results/<subtype>/<segment>/curated_reference.gff`: GFF that matches the curated alignment
   - `results/<subtype>/<segment>/curated_msa.vcf.gz`: VCF format for UShER
   - `results/<subtype>/<segment>/preopt_tree.pb.gz`: Initial UShER tree
   - `results/<subtype>/<segment>/opt_tree.pb.gz`: Optimized tree
   - `results/<subtype>/<segment>/opt_tree.jsonl.gz`: Taxonium visualization file

## Pipeline Steps

1. **Parse GISAID data**: Splits combined FASTA files into segment-specific files and creates metadata file
2. **Download reference**: Fetches reference sequences and creates Nextclade dataset
3. **Align sequences**: Uses Nextclade to align sequences to reference
4. **Curate alignment**: Extracts coding regions and sanitizes sequence IDs
5. **Create VCF**: Converts alignment to VCF format for UShER
6. **Build initial tree**: Uses UShER to create a parsimony-based tree
7. **Optimize tree**: Uses matOptimize to refine the tree
8. **Create visualization**: Converts tree to Taxonium format for interactive viewing

## Requirements

All dependencies are specified in the environment.yml file.