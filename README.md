# Flu-UShER Pipeline

A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER.

## Directory Structure

```
flu-usher/
├── Snakefile                # Main pipeline file
├── config/
│   └── config.yaml          # Configuration file
├── data/                    # Input data (organized by type-segment)
│   ├── H7N9-HA/             # Example input directory
│   │   └── H7N9-HA.fasta    # Input sequences
│   └── ...
├── logs/                    # Log files (created by the pipeline)
├── results/                 # Output results (organized by type-segment)
├── scripts/
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
   - List your types and segments to analyze
   - Set reference accession numbers for each type-segment combination
   - Adjust filtering thresholds for sequence curation
   - Set desired number of threads

3. **Prepare your input data**

   For each type-segment combination you want to analyze (e.g., H7N9-HA), create:
   - A directory under `data/` with the name `{type}-{segment}` (e.g., `data/H7N9-HA/`)
   - A FASTA file inside that directory with the same name as the directory (e.g., `data/H7N9-HA/H7N9-HA.fasta`)

4. **Run the pipeline**

   ```bash
   # Test run (dry-run)
   snakemake -np
   
   # Run the pipeline
   snakemake --cores <number_of_cores>
   
   # Run for specific type-segment combinations
   snakemake --cores <number_of_cores> results/H7N9-HA/opt_tree.pb.gz
   ```

5. **Output**

   For each type-segment combination, the pipeline produces:
   - `results/<type>-<segment>/reference/`: Reference data for Nextclade
   - `results/<type>-<segment>/msa.fasta.xz`: Aligned sequences
   - `results/<type>-<segment>/curated_msa.fasta.xz`: Curated alignment
   - `results/<type>-<segment>/curated_msa.vcf.gz`: VCF format for UShER
   - `results/<type>-<segment>/preopt_tree.pb.gz`: Initial UShER tree
   - `results/<type>-<segment>/opt_tree.pb.gz`: Optimized tree
   - `results/<type>-<segment>/opt_tree.jsonl.gz`: Taxonium visualization file

## Pipeline Steps

1. **Download reference**: Automatically fetches reference sequences and associated GFF files and creates Nextclade dataset
2. **Align sequences**: Uses Nextclade to align sequences to the reference
3. **Curate alignment**: Extracts coding regions and sanitizes sequence IDs
4. **Create VCF**: Converts alignment to VCF format for UShER
5. **Build initial tree**: Uses UShER to create a parsimony-based tree
6. **Optimize tree**: Uses matOptimize to refine the tree
7. **Create visualization**: Converts tree to Taxonium format for interactive viewing

## Requirements

All dependencies are specified in the environment.yml file.