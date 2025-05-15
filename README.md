# Flu-UShER Pipeline

A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER.

## TODO
* View tree using taxonium
* Incorporating metadata
   * When downloading, only include EPI_ISL, subtype, segment
   * Download all segments at once, then add a step that splits them into separate files?
   * In curation step, simplify header to only include EPI_ISL
   * Add CSV of metadata file with a specific subset of columns
* Update pipeline to read in gene start and end from a gff3 file, rather than defining these things in the config
* Apply to other segments and subtypes

## Directory Structure

```
flu-usher/
├── Snakefile                # Main pipeline file
├── config/
│   └── config.yaml          # Configuration file
├── data/                    # Input data (organized by type.segment)
│   ├── h7n9.ha/             # Example input directory
│   │   ├── sequences.fasta  # Input sequences
│   │   └── ...              # Nextclade reference data
│   └── ...
├── logs/                    # Log files (created by the pipeline)
├── results/                 # Output results (organized by type.segment)
└── scripts/
    └── curate_msa.py        # Sequence curation script
```

## Usage

1. **Set up your data**

   For each type-segment combination (e.g., h7n9.ha), create a directory under `data/` containing:
   - `sequences.fasta`: Input sequences
   - Nextclade reference dataset files

2. **Configure the pipeline**

   Edit `config/config.yaml` to:
   - List your types and segments
   - Set gene-specific parameters
   - Adjust filtering thresholds
   - Set desired number of threads

3. **Run the pipeline**

   ```bash
   # Test run (dry-run)
   snakemake -np
   
   # Run the pipeline
   snakemake --cores <number_of_cores>
   
   # Run for specific type-segment combinations
   snakemake --cores <number_of_cores> results/h7n9.ha/opt_tree.pb.gz
   ```

4. **Output**

   For each type-segment combination, the pipeline produces:
   - `results/<type>.<segment>/msa.fasta.xz`: Aligned sequences
   - `results/<type>.<segment>/curated_msa.fasta.xz`: Curated alignment
   - `results/<type>.<segment>/curated_msa.vcf.gz`: VCF format for UShER
   - `results/<type>.<segment>/preopt_tree.pb.gz`: Initial UShER tree
   - `results/<type>.<segment>/opt_tree.pb.gz`: Optimized tree

## Requirements

- Snakemake
- Nextclade
- UShER
- matOptimize
- faToVcf (from UCSC)
- Python with Biopython