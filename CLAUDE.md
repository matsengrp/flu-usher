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

Per-step dependencies are managed via separate conda environments in `envs/` (ete3, fatovcf, historydag, nextclade, pastml, python, taxonium, usher; larch is not among them -- see the note on `tree_to_dag`). Snakemake creates and activates these automatically when run with `--use-conda`.

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
```

### Workflow Control
```bash
# ALWAYS put file targets BEFORE these flags. --forcerun, --rerun-triggers,
# --until, --omit-from and --allowed-rules all take nargs='+' and swallow every
# following argument up to the next flag. `snakemake --forcerun <rule> <target>`
# silently absorbs <target> as another rule name, leaving no positional target,
# and falls back to building `rule all` -- which on interrupt deletes the
# declared outputs of whatever was running. --rerun-triggers at least rejects
# the target with "invalid choice"; --forcerun fails silently.
snakemake --cores 8 --use-conda <target> --forcerun <rule_name>

# `--rerun-triggers mtime` narrows the default trigger set (mtime, params,
# input, software-env, code) down to mtime alone, disabling params-based
# invalidation. That is safe for config.yaml, contrary to what this file used
# to say: download_all_references passes config.yaml to its script, so it is a
# real input, and it sits at the head of the DAG. Editing any config value
# changes that file's mtime, re-runs it, and cascades through the reference
# files to everything -- measured, one reroot target changed re-runs 1235 of
# 1235 jobs, identically with and without the params trigger.
#
# The flip side is that config.yaml invalidation is maximally blunt: there is
# no such thing as editing one config value and re-running only what depends
# on it. Expect a full re-derivation from any config edit.
snakemake -n --use-conda <target> --rerun-triggers mtime

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
│   │   ├── unaligned_coding_seqs/               # one file per gene
│   │   │   └── curated_unaligned_{GENE}.fasta.xz
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
│   │   ├── rerooted_tree.nh                     # only where `reroot` is configured
│   │   ├── reference_origin_tree.pb.gz          # MAT before rebasing onto its root
│   │   ├── final_tree.{pb.gz,jsonl.gz}
│   │   ├── final_tree_root.fasta                # the root sequence; final_tree's origin
│   │   ├── reroot_sequence_check.txt            # gate report; seqcheck/ inputs are temp()
│   │   ├── geographic_trees/
│   │   │   ├── {geo_group}_samples.txt
│   │   │   ├── {geo_group}_tree.pb.gz
│   │   │   └── {geo_group}_tree.jsonl.gz
│   │   ├── host_ancestral/                     # PastML per-node host_group
│   │   │   ├── combined_ancestral_states.tab
│   │   │   └── host_tree.html
│   │   └── subtype_ancestral/                  # PastML per-node subtype (H*N*)
│   │       ├── combined_ancestral_states.tab
│   │       └── subtype_tree.html
│   ├── H3/
│   ├── H5/
│   ├── H7/
│   └── H9/
├── NA/          # NA segment results by subtype
│   ├── N1/      # Same structure as HA subtypes
│   ├── N2/
│   ├── N6/
│   ├── N8/
│   └── N9/
└── {PB2,PB1,PA,NP,MP,NS}/  # Internal segments
    └── all/     # All subtypes combined (same structure)
```

### Key Components

1. **Snakefile**: Main workflow definition that orchestrates the entire pipeline
   - Defines rules for each processing step
   - Manages dependencies between steps
   - Handles parallel execution

2. **config.yaml** (repo root, loaded by `Snakefile:3`): Configuration file containing:
   - Input directories for GISAID data
   - HA/NA subtypes to analyze
   - Reference sequences for each segment-subtype
   - Quality filtering thresholds (max_frac_gaps, max_frac_ambig)
   - Number of randomizations for tree building (n_randomizations)
   - Number of threads for parallel execution
   - Geographic groups to extract for geographic subtree analysis (geographic_groups_to_extract)
   - Optional rerooting specifications for final trees (reroot)

3. **scripts/**: Python scripts for specific tasks:
   - `parse_gisaid_data.py`: Parses GISAID FASTA/metadata files and organizes by segment/subtype
   - `download_ref_seq.py`: Downloads reference sequences from NCBI and creates Nextclade datasets
   - `curate_and_extract_coding_seqs.py`: Curates alignments by quality metrics, extracts coding regions and per-gene unaligned coding sequences
   - `randomize_alignment.py`: Creates randomized versions of alignments for multiple tree builds
   - `trim_dag.py`: Trims suboptimal trees from merged DAGs
   - `convert_DAG_protobuf_to_newick_samples.py`: Samples representative trees from DAGs
   - `reroot_newick.py`: Reroots a newick at a named leaf via ete3 `set_outgroup`; validates the target exists, is unique, is a leaf, and ends up at the root
   - `check_tree_sequences.py`: Asserts no sample's sequence changed across rerooting, and that the final tree is rooted at the configured target
   - `rebase_mat_root.py`: Moves a MAT's origin from the reference onto its own root, emptying the root's mutation list and writing the root sequence out
   - `create_root_samples_file.py`: Creates sample files for root sequence extraction
   - `extract_root_sequence.py`: Infers root sequences from tree mutations
   - `simplified_host_classifier.py`: Host classification logic (used by augment_metadata.py)
   - `augment_metadata.py`: Adds host_group, geographic_group, and temporal_group columns to metadata
   - `create_samples_file.py`: Creates sample files for subtree extraction by any metadata column
   - `prepare_host_annotation.py`: Builds the global 2-column (isolate_id, host_group) CSV consumed by PastML
   - `prepare_subtype_annotation.py`: Builds the global 2-column (isolate_id, subtype) CSV consumed by PastML, normalizing the raw GISAID `subtype` ("A / H5N1") to `H*N*` form inline
   - `utils.py`: Shared helpers imported by the above — `open_sequence_file()` for plain/gz/xz IO, `setup_logging()`, `sanitize_id()`, and the GFF parsing used to derive coding coordinates
   - `test_curate_and_extract_coding_seqs.py`, `test_parse_gisaid_data.py`, `test_create_samples_file.py`, `test_check_tree_sequences.py`, `test_reroot_newick.py`: unittest suites

4. **notebooks/**: Jupyter notebooks for analysis and development
   - `analyze_alignments.ipynb`: Analyzes sequence statistics across segments/subtypes
   - Note: Notebooks may need updates to work with the new directory structure

### Pipeline Workflow

1. **Parse GISAID Data** → Aggregates sequences from multiple input directories, splits by segment/subtype
2. **Download References** → Fetches appropriate reference sequences for each segment-subtype
3. **Align Sequences** → Uses Nextclade for codon-aware alignment
4. **Curate Alignment** → Filters by quality (`max_frac_gaps: 0.03`, `max_frac_ambig: 0.00` in `config.yaml`) and extracts coding regions
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
15. **Reroot Tree** → Where `reroot` is configured, reroots the *newick* at the named leaf with ete3 `set_outgroup` (`reroot_newick`), then rebuilds the MAT from the alignment with `usher` (`create_final_mat`), which keeps its mutations recorded against the reference. Combinations with no configured reroot symlink `reference_origin_tree.pb.gz` to `sampled_tree.pb.gz`; step 16 then runs on all 16 alike, so `final_tree.pb.gz` is always a real file and always means the same thing. Replaced `matUtils extract -y`, which moved the MAT's origin onto the new root and refused outright when the pre-reroot root carried a mutation (issue #49). `check_tree_sequences` then guards that the MAT is still recorded against the reference, and asserts the tree is rooted where config asked — the sequence comparison alone cannot tell, since both trees are built from the same VCF and so agree for any topology.

    Both MAT-building rules (`create_mat_protobuf`, `create_final_mat`) use `usher`, not `matOptimize`. These steps annotate a fixed topology; they do not search for a better one. `matOptimize` is a parsimony optimiser, and parsimony is an *unrooted* criterion, so it normalises the tree at load even under `-N 0` and does not preserve the input rooting. That silently broke NA/N1: its reroot target `EPI_ISL_5878` is the only one of the 14 whose terminal branch carries no mutations (the others run 6–87 in `sampled_tree.pb.gz`, 5–62 after rerooting), and it was the only one that lost its rooting. Empirically the rooting survives `matOptimize` only when the target's terminal branch carries mutations; the exact normalisation step responsible has not been traced to matOptimize's source, and the observed failure is not a simple collapse — NA/N1's root landed 54 mutations away under a 3-way polytomy, 56 substitutions from the intended root sequence, not merged into its neighbour. Switching both rules to `usher` fixed it and moved `curated_root.fasta` by one position in 2 of 16 combinations. `matOptimize` remains correct for `optimize_tree`, which really is topology optimisation.
16. **Rebase MAT onto its Root** → `rebase_final_mat` moves the tree's origin from `curated_reference.fasta` onto the tree's own root, so the root carries no mutations and its sequence ships as `final_tree_root.fasta`. Reconstructing a sample from `final_tree.pb.gz` therefore starts from that file, not from the curated reference. This is pure bookkeeping: every branch below the root already records a parent-to-child change and is unaffected, so only the root's own mutation list (emptied) and each mutation's `ref_nuc` annotation (repointed) change. It is the same operation `matUtils extract -y` performed, minus its two defects — refusing whenever the root carried mutations, and discarding the root sequence instead of emitting it, which is what made the shift silent in issue #49.
17. **Create Root Sequence** → Infers root sequence from tree or uses reference
18. **Augment Metadata** → Adds host_group, geographic_group, and temporal_group columns
19. **Extract Subtrees** → Creates subtrees for each configured geographic region (matUtils extract)
20. **Infer Per-Node Host States** → Runs PastML / DOWNPASS on each `final_tree.nwk` to reconstruct `host_group` at every node; outputs `{segment}/{subtype}/host_ancestral/combined_ancestral_states.tab`
21. **Infer Per-Node Subtype States** → Runs PastML / DOWNPASS on each `final_tree.nwk` to reconstruct `subtype` (`H*N*`, normalized from the raw GISAID `subtype` inside `prepare_subtype_annotation.py`) at every node; outputs `{segment}/{subtype}/subtype_ancestral/combined_ancestral_states.tab`. On HA per-subtype trees the H is fixed and only the inferred N partner varies (analogously for NA trees); on internal-segment trees neither letter is constrained, so the full `H*N*` can vary along the tree.
22. **Create Visualizations** → Generates Taxonium format for full tree and geographic subtrees
23. **Execute Notebooks** → Runs the notebooks under `notebooks/` once all 16 `final_tree.jsonl.gz` files exist; writes `results/notebooks/{notebook}.done` sentinels
24. **Record Input md5sums** → Writes `results/input_data_md5sums.txt`, a provenance manifest over every FASTA/XLS file discovered under the configured `input_dirs`

### Input Data Requirements

The pipeline expects GISAID data in each input directory:
- FASTA files with sequence ID format: `EPI|SEGMENT|NAME|EPI_ISL|SUBTYPE`
- XLS metadata files with corresponding sequence information
- The `subtype` column in metadata determines HA/NA grouping

### Important Notes

- Tests live in seven `scripts/test_*.py` modules (202 unittest tests: curate 91, check_tree_sequences 31, rebase_mat_root 30, utils 18, parse_gisaid 15, create_samples 9, reroot_newick 8). Run them with `python -m unittest discover` from within `scripts/`. Under `envs/python.yaml` that reports `OK (skipped=8)`: `test_reroot_newick.py` is `skipIf`-guarded because ete3 lives in its own env, so run it separately under `envs/ete3.yaml` to actually exercise those 8 — `OK (skipped=N)` otherwise reads like a pass. No linting setup currently exists, and 10 of the 17 script modules are untested.
- The pipeline uses compressed outputs (.xz, .gz) to save disk space
- All logs are saved in the `logs/` directory organized by segment/subtype
- The pipeline can process multiple influenza subtypes simultaneously
- Reference sequences are specified in config.yaml and can be customized
- The pipeline uses a DAG-based approach (via larch and historydag) to build consensus trees from multiple randomized alignments
- Multiple randomizations help explore tree space and produce more robust phylogenies
- Subtrees are automatically extracted for specified geographic regions (e.g., north_america, europe, asia)
- Trees can be optionally rerooted using the `reroot` configuration parameter
- The final outputs are interactive Taxonium visualization files (.jsonl.gz)