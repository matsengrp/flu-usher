# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Flu-UShER Pipeline - A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER. The pipeline processes influenza sequences by segment and subtype, organizing HA and NA segments by their specific subtypes (e.g., H1, H3, N1, N2) while combining all subtypes for internal segments (PB2, PB1, PA, NP, MP, NS).

## Key Commands

### Environment Setup
```bash
# Create the main conda environment (provides Snakemake and larch build deps)
conda env create -f environment.yml

# Activate the environment
conda activate flu-usher
```

Per-step dependencies are managed via separate conda environments in `envs/` (fatovcf, historydag, larch, nextclade, python, taxonium, usher). Snakemake creates and activates these automatically when run with `--use-conda`.

### Running the Pipeline
```bash
# Test run (dry-run) to see what will be executed
snakemake -np --use-conda

# Run the full pipeline with specified cores
snakemake --cores 12 --use-conda

# Run for specific segment-subtype combinations
snakemake --cores 8 --use-conda results/HA/H5/final_tree.jsonl.gz
snakemake --cores 8 --use-conda results/NA/N1/final_tree.jsonl.gz
snakemake --cores 8 --use-conda results/PB2/all/final_tree.jsonl.gz

# Run for geographic subtrees
snakemake --cores 8 --use-conda results/HA/H5/geographic_trees/north_america_tree.jsonl.gz

# Run for temporal subtrees
snakemake --cores 8 --use-conda results/HA/H5/temporal_trees/early_tree.jsonl.gz
```

### Workflow Control
```bash
# Force rerun from a specific rule
snakemake --forcerun <rule_name> --cores 8 --use-conda

# Generate workflow visualization
snakemake --dag | dot -Tpdf > workflow.pdf
```

## Architecture and Directory Structure

### Reorganized Pipeline Structure (Current)
The pipeline was recently reorganized from a subtype-first to a segment-first structure:

```
results/
├── combined_metadata.csv                    # Aggregated metadata
├── combined_metadata_augmented.csv          # Metadata with host_group, geographic_group, temporal_group
├── HA/          # HA segment results by subtype
│   ├── H1/      # Individual subtype results
│   │   ├── raw_sequences.fasta.xz
│   │   ├── reference/
│   │   ├── msa.fasta.xz
│   │   ├── curated_msa.fasta.xz
│   │   ├── curated_unaligned_coding_seqs.fasta.xz
│   │   ├── curated_reference.{fasta,txt,gff,gtf}
│   │   ├── curated_root.fasta
│   │   ├── randomized_{0,1,2,...}/  # Multiple randomizations
│   │   │   ├── msa.fasta.xz
│   │   │   ├── msa.vcf
│   │   │   ├── preopt_tree.pb.gz
│   │   │   ├── opt_tree.pb.gz
│   │   │   └── dag.pb
│   │   ├── larch_merged_dag.pb
│   │   ├── trimmed_dag.pb
│   │   ├── sampled_tree.{nh,pb.gz}
│   │   ├── final_tree.{pb.gz,jsonl.gz}
│   │   ├── geographic_trees/
│   │   │   ├── {geo_group}_samples.txt
│   │   │   ├── {geo_group}_tree.pb.gz
│   │   │   └── {geo_group}_tree.jsonl.gz
│   │   └── temporal_trees/
│   │       ├── {temporal_group}_samples.txt
│   │       ├── {temporal_group}_tree.pb.gz
│   │       └── {temporal_group}_tree.jsonl.gz
│   ├── H3/
│   ├── H5/
│   ├── H7/
│   └── H9/
├── NA/          # NA segment results by subtype
│   ├── N1/      # Same structure as HA subtypes
│   ├── N2/
│   └── N9/
└── {PB2,PB1,PA,NP,MP,NS}/  # Internal segments
    └── all/     # All subtypes combined (same structure)
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
   - Quality filtering thresholds (max_frac_gaps, max_frac_ambig)
   - Number of randomizations for tree building (n_randomizations)
   - Number of threads for parallel execution
   - Geographic groups to extract for geographic subtree analysis (geographic_groups_to_extract)
   - Temporal groups to extract for temporal subtree analysis (temporal_groups_to_extract)
   - Optional rerooting specifications for final trees (reroot)

3. **scripts/**: Python scripts for specific tasks:
   - `parse_gisaid_data.py`: Parses GISAID FASTA/metadata files and organizes by segment/subtype
   - `download_ref_seq.py`: Downloads reference sequences from NCBI and creates Nextclade datasets
   - `curate_and_extract_coding_seqs.py`: Curates alignments by quality metrics, extracts coding regions and per-gene unaligned coding sequences
   - `randomize_alignment.py`: Creates randomized versions of alignments for multiple tree builds
   - `trim_dag.py`: Trims suboptimal trees from merged DAGs
   - `convert_DAG_protobuf_to_newick_samples.py`: Samples representative trees from DAGs
   - `create_root_samples_file.py`: Creates sample files for root sequence extraction
   - `extract_root_sequence.py`: Infers root sequences from tree mutations
   - `simplified_host_classifier.py`: Host classification logic (used by augment_metadata.py)
   - `augment_metadata.py`: Adds host_group, geographic_group, and temporal_group columns to metadata
   - `create_samples_file.py`: Creates sample files for subtree extraction by any metadata column
   - `create_temporal_samples_file.py`: Creates sample files for temporal subtree extraction (per-tree median date split)

4. **notebooks/**: Jupyter notebooks for analysis and development
   - `analyze_alignments.ipynb`: Analyzes sequence statistics across segments/subtypes
   - Note: Notebooks may need updates to work with the new directory structure

### Pipeline Workflow

1. **Parse GISAID Data** → Aggregates sequences from multiple input directories, splits by segment/subtype
2. **Download References** → Fetches appropriate reference sequences for each segment-subtype
3. **Align Sequences** → Uses Nextclade for codon-aware alignment
4. **Curate Alignment** → Filters by quality (gaps < 5%, ambiguities < 1%) and extracts coding regions
5. **Create Unaligned Coding Sequences** → Extracts unaligned coding sequences from curated alignments
6. **Randomize Alignments** → Creates multiple randomized versions of alignment (n_randomizations)
7. **Create VCF** → Converts each randomized FASTA to variant format for UShER
8. **Build Initial Trees** → Creates initial parsimony tree for each randomization with usher-sampled
9. **Optimize Trees** → Refines topology for each tree with matOptimize
10. **Convert to DAGs** → Converts each optimized tree to DAG representation (larch-usher)
11. **Merge DAGs** → Combines all DAGs into single merged DAG (larch-dagutil)
12. **Trim DAG** → Removes suboptimal trees from merged DAG
13. **Sample Tree** → Samples a representative tree from trimmed DAG
14. **Create MAT Protobuf** → Converts sampled tree to MAT protobuf format
15. **Reroot Tree** → Optionally reroots tree at specified node (matUtils extract)
16. **Create Root Sequence** → Infers root sequence from tree or uses reference
17. **Augment Metadata** → Adds host_group, geographic_group, and temporal_group columns
18. **Extract Subtrees** → Creates subtrees for each geographic region and temporal period (matUtils extract)
19. **Create Visualizations** → Generates Taxonium format for full tree and all subtrees

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
- The pipeline uses a DAG-based approach (via larch and historydag) to build consensus trees from multiple randomized alignments
- Multiple randomizations help explore tree space and produce more robust phylogenies
- Subtrees are automatically extracted for specified geographic regions (e.g., north_america, europe, asia) and temporal periods (early/late split at per-tree median collection date)
- Trees can be optionally rerooted using the `reroot` configuration parameter
- The final outputs are interactive Taxonium visualization files (.jsonl.gz)