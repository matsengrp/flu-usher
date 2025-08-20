# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Flu-UShER Pipeline - A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER. The pipeline processes influenza sequences by segment and subtype, organizing HA and NA segments by their specific subtypes (e.g., H1, H3, N1, N2) while combining all subtypes for internal segments (PB2, PB1, PA, NP, MP, NS).

## Key Commands

### Environment Setup
```bash
# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate usher
```

### Running the Pipeline
```bash
# Test run (dry-run) to see what will be executed
snakemake -np

# Run the full pipeline with specified cores
snakemake --cores 12

# Run for specific segment-subtype combinations
snakemake --cores 8 results/HA/H5/opt_tree.pb.gz
snakemake --cores 8 results/NA/N1/opt_tree.pb.gz
snakemake --cores 8 results/PB2/all/opt_tree.pb.gz
```

### Workflow Control
```bash
# Force rerun from a specific rule
snakemake --forcerun <rule_name> --cores 8

# Generate workflow visualization
snakemake --dag | dot -Tpdf > workflow.pdf
```

## Architecture and Directory Structure

### Reorganized Pipeline Structure (Current)
The pipeline was recently reorganized from a subtype-first to a segment-first structure:

```
results/
├── HA/          # HA segment results by subtype
│   ├── H1/      # Individual subtype results
│   ├── H3/
│   ├── H5/
│   ├── H7/
│   └── H9/
├── NA/          # NA segment results by subtype
│   ├── N1/
│   ├── N2/
│   └── N9/
└── {PB2,PB1,PA,NP,MP,NS}/  # Internal segments
    └── all/     # All subtypes combined
```

### Key Components

1. **Snakefile**: Main workflow definition that orchestrates the entire pipeline
   - Defines rules for each processing step
   - Manages dependencies between steps
   - Handles parallel execution

2. **config/config.yaml**: Configuration file containing:
   - Input directories for GISAID data
   - HA/NA subtypes to analyze
   - Reference sequences for each segment-subtype
   - Quality filtering thresholds
   - Maximum parsimony score for branch pruning (default: 20)

3. **scripts/**: Python scripts for specific tasks:
   - `parse_gisaid_data.py`: Parses GISAID FASTA/metadata files and organizes by segment/subtype
   - `download_ref_seq.py`: Downloads reference sequences from NCBI and creates Nextclade datasets
   - `curate_msa.py`: Filters alignments by quality metrics and extracts coding regions

4. **notebooks/**: Jupyter notebooks for analysis and development
   - `analyze_alignments.ipynb`: Analyzes sequence statistics across segments/subtypes
   - Note: Notebooks may need updates to work with the new directory structure

### Pipeline Workflow

1. **Parse GISAID Data** → Aggregates sequences from multiple input directories, splits by segment/subtype
2. **Download References** → Fetches appropriate reference sequences for each segment-subtype
3. **Align Sequences** → Uses Nextclade for codon-aware alignment
4. **Curate Alignment** → Filters by quality (gaps < 5%, ambiguities < 1%)
5. **Create VCF** → Converts FASTA to variant format for UShER
6. **Build Tree** → Creates initial parsimony tree with usher-sampled
7. **Prune Tree** → Removes leaf nodes with very long branches (>20 parsimony score)
8. **Optimize Tree** → Refines topology with matOptimize
9. **Create Visualization** → Generates Taxonium format for interactive exploration

### Input Data Requirements

The pipeline expects GISAID data in each input directory:
- FASTA files with sequence ID format: `EPI|SEGMENT|NAME|EPI_ISL|SUBTYPE`
- XLS metadata files with corresponding sequence information
- The `subtype` column in metadata determines HA/NA grouping

### Important Notes

- No formal testing infrastructure or linting setup currently exists
- The pipeline uses compressed outputs (.xz, .gz) to save disk space
- All logs are saved in the `logs/` directory organized by segment/subtype
- The pipeline can process multiple influenza subtypes simultaneously
- Reference sequences are specified in config.yaml and can be customized