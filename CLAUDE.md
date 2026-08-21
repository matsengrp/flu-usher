# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Flu-UShER Pipeline - A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER. The pipeline processes influenza sequences by segment and subtype, organizing HA and NA segments by their specific subtypes (e.g., H1, H3, N1, N2) while combining all subtypes for internal segments (PB2, PB1, PA, NP, MP, NS).

Design rationale — why the pipeline is built the way it is, with the measurements
behind each decision — lives in the "Design Notes" section of `README.md`. The
Snakefile's own comments say only what each rule does and point there. Keep it
that way: when you change a rule's reasoning, update `README.md`, not the
Snakefile comment.

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
# input, software-env, code) down to mtime alone, which disables params-based
# invalidation. Only execute_analyze_alignments and execute_analyze_dags still
# declare config.yaml as an input -- their notebooks open it directly, which
# nothing but an mtime dependency can track -- so under mtime-only a config
# edit re-runs those two and nothing else, and every substantive output goes
# quietly stale. Measured: a reroot edit under mtime-only schedules 3 jobs.
# Use it for a quick dry run, not to decide whether real work is up to date.
#
# Under the default triggers, config invalidation is precise: each rule names
# the config values it consumes in `params:`, so editing a reroot target
# re-runs the reroot rules and leaves the alignments alone. It did not used to
# be -- download_all_references took config.yaml as an input at the head of the
# DAG, so any edit re-downloaded the references and cascaded through
# everything. Measured before that change, on the 1240-job DAG of the time: one
# reroot target re-ran 1235 of 1240 jobs, identically with and without the
# params trigger. (The five spared were parse_gisaid_data, augment_metadata, the
# two annotation rules and input_data_md5sums -- everything expensive
# re-derived.) The DAG is 1720 jobs since the scaffold rules added three
# 160-job rules; count it with `snakemake -n --forceall`, as a plain `-n`
# reports only the jobs that are out of date.
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
│   │   │   ├── scaffold_msa.fasta.xz            # time-spread subset
│   │   │   ├── scaffold_samples.txt             # ids drawn into it
│   │   │   ├── scaffold_tree.pb.gz              # backbone the full placement seeds from
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
   - Number of taxa drawn into each randomization's scaffold tree (scaffold_n_taxa), and the most any single collection year may contribute to it (scaffold_max_year_fraction)
   - Seed for sampling a tree out of the trimmed DAG (tree_sample_seed)
   - Number of threads for parallel execution
   - Geographic groups to extract for geographic subtree analysis (geographic_groups_to_extract)
   - Optional rerooting specifications for final trees (reroot)

3. **scripts/**: Python scripts for specific tasks:
   - `parse_gisaid_data.py`: Parses GISAID FASTA/metadata files and organizes by segment/subtype. Keeps `collection_date` as the string GISAID supplied, at whichever of its three precisions that is, and fails the run on a value that is not a real date at one of them (issue #55; see "Collection dates" in the Design Notes of `README.md`)
   - `download_ref_seq.py`: Downloads reference sequences from NCBI and creates Nextclade datasets
   - `curate_and_extract_coding_seqs.py`: Curates alignments by quality metrics, extracts coding regions and per-gene unaligned coding sequences
   - `randomize_alignment.py`: Creates randomized versions of alignments for multiple tree builds
   - `create_scaffold_alignment.py`: Subsets an alignment to `scaffold_n_taxa` sequences spread evenly over collection year, capping any one year's contribution at `scaffold_max_year_fraction` of what it holds, for the backbone tree each randomization is seeded with
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
   - `test_*.py`: ten unittest suites, one per tested module (see Important Notes for counts)

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
8. **Build Initial Trees** → Each randomization first draws a `scaffold_n_taxa` subset spread evenly over collection year, with no year giving more than `scaffold_max_year_fraction` of what it holds (`create_scaffold_alignment`), builds and optimizes a backbone tree from it alone (`build_scaffold_tree`), then places every remaining sequence onto that backbone with usher-sampled (`create_initial_tree`, `-i` not `-t`). Placing all 86,232 HA/H3 sequences onto an *empty* tree instead settled on a measurably worse high-level shape — issue #53. Why the backbone exists, how the draw works, and the two measurement traps involved: see the "Design Notes" section of `README.md`.
9. **Optimize Trees** → Refines topology for each tree with matOptimize
10. **Convert to DAGs** → Converts each optimized tree to DAG representation (larch-usher)
11. **Merge DAGs** → Combines all DAGs into single merged DAG (larch-dagutil)
12. **Trim DAG** → Removes suboptimal trees from merged DAG
13. **Sample Tree** → Samples a representative tree from trimmed DAG
14. **Create MAT Protobuf** → Converts sampled tree to MAT protobuf format
15. **Reroot Tree** → Where `reroot` is configured, reroots the *newick* at the named leaf with ete3 `set_outgroup` (`reroot_newick`), then rebuilds the MAT from the alignment with `usher` (`create_final_mat`), which keeps its mutations recorded against the reference. Combinations with no configured reroot symlink `reference_origin_tree.pb.gz` to `sampled_tree.pb.gz`; step 16 then runs on all 16 alike, so `final_tree.pb.gz` is always a real file and always means the same thing. Replaced `matUtils extract -y` (issue #49). `check_tree_sequences` then guards that the MAT is still recorded against the reference, and asserts the tree is rooted where config asked. See "Rerooting" and "The sequence-identity gate" in the Design Notes of `README.md`.

    Both MAT-building rules (`create_mat_protobuf`, `create_final_mat`) use `usher`, not `matOptimize`, because they annotate a fixed topology and must keep it, rooting included. `matOptimize` remains correct for `optimize_tree`, which really is topology optimisation. Why, and the NA/N1 failure that established it: see "Why `usher` and not `matOptimize`" in the Design Notes of `README.md`.
16. **Rebase MAT onto its Root** → `rebase_final_mat` moves the tree's origin from `curated_reference.fasta` onto the tree's own root, so the root carries no mutations and its sequence ships as `final_tree_root.fasta`. Reconstructing a sample from `final_tree.pb.gz` therefore starts from that file, not from the curated reference. This is pure bookkeeping: only the root's own mutation list (emptied) and each mutation's `ref_nuc` annotation (repointed) change. It is the same operation `matUtils extract -y` performed, minus its two defects — see "Rerooting" in the Design Notes of `README.md`.
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

- Tests live in ten `scripts/test_*.py` modules (275 unittest tests: curate 91, create_scaffold_alignment 44, check_tree_sequences 31, rebase_mat_root 30, parse_gisaid_data 25, utils 21, extract_root_sequence 13, create_samples 9, reroot_newick 8, config_params_sync 3). Run them with `python -m unittest discover` from within `scripts/`. Under `envs/python.yaml` that reports `OK (skipped=8)`: `test_reroot_newick.py` is `skipIf`-guarded because ete3 lives in its own env, so run it separately under `envs/ete3.yaml` to actually exercise those 8 — `OK (skipped=N)` otherwise reads like a pass. No linting setup currently exists, and 9 of the 18 script modules are untested (augment_metadata, convert_DAG_protobuf_to_newick_samples, create_root_samples_file, download_ref_seq, prepare_host_annotation, prepare_subtype_annotation, randomize_alignment, simplified_host_classifier, trim_dag).
- The pipeline uses compressed outputs (.xz, .gz) to save disk space
- All logs are saved in the `logs/` directory organized by segment/subtype
- The pipeline can process multiple influenza subtypes simultaneously
- Reference sequences are specified in config.yaml and can be customized
- The pipeline uses a DAG-based approach (via larch and historydag) to build consensus trees from multiple randomized alignments
- Multiple randomizations help explore tree space and produce more robust phylogenies
- Subtrees are automatically extracted for specified geographic regions (e.g., north_america, europe, asia)
- Trees can be optionally rerooted using the `reroot` configuration parameter
- The final outputs are interactive Taxonium visualization files (.jsonl.gz)