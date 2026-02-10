# Flu-UShER Pipeline

A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER. The pipeline processes influenza sequences by segment and subtype, creating separate trees for HA and NA segments by subtype (e.g., H1, H3, N1, N2) while combining all subtypes for other segments (PB2, PB1, PA, NP, MP, NS).

## Directory Structure

```
flu-usher/
├── Snakefile                # Main pipeline file
├── config.yaml              # Configuration file
├── data/                    # Input data directories
│   ├── H1N1/                # Example: H1N1 sequences (all segments)
│   ├── H3N2/                # Example: H3N2 sequences (all segments)
│   ├── H5N1/                # Example: H5N1 sequences (all segments)
│   └── H7N9/                # Example: H7N9 sequences (all segments)
│       ├── sequences.fasta  # One or more FASTA files containing sequences
│       └── metadata.xls     # One or more metadata files from GISAID
├── logs/                    # Log files (created by the pipeline)
├── results/                 # Output results (organized by segment/subtype)
│   ├── HA/                  # HA segment results by subtype
│   │   ├── H1/              # H1 subtype tree and files
│   │   ├── H3/              # H3 subtype tree and files
│   │   ├── H5/              # H5 subtype tree and files
│   │   ├── H7/              # H7 subtype tree and files
│   │   └── H9/              # H9 subtype tree and files
│   ├── NA/                  # NA segment results by subtype
│   │   ├── N1/              # N1 subtype tree and files
│   │   ├── N2/              # N2 subtype tree and files
│   │   └── N9/              # N9 subtype tree and files
│   └── PB2/                 # Other segment results (PB1, PA, NP, MP, NS)
│       └── all/             # All subtypes combined
├── scripts/
│   ├── parse_gisaid_data.py                        # Parse GISAID data by segment
│   ├── download_ref_seq.py                         # Download reference sequences
│   ├── curate_and_extract_coding_seqs.py           # Curate alignments and extract coding sequences
│   ├── randomize_alignment.py                      # Randomize alignment order
│   ├── trim_dag.py                                 # Trim suboptimal trees from DAG
│   ├── convert_DAG_protobuf_to_newick_samples.py   # Sample tree from DAG
│   ├── create_root_samples_file.py                 # Create samples file for root extraction
│   ├── extract_root_sequence.py                    # Infer root sequence from tree
│   ├── simplified_host_classifier.py               # Classify hosts into groups
│   └── create_host_samples_file.py                 # Create samples file for host extraction
└── notebooks/               # Jupyter notebooks for development and analysis
```

## Usage

1. **Set up your environment**

   ```
   # Create the conda environment
   conda env create -f environment.yml
   
   # Activate the environment
   conda activate usher
   ```

2. **Configure the pipeline**

   Edit `config.yaml` to:
   - Specify input directories containing your GISAID data
   - List HA subtypes to analyze (e.g., H1, H3, H5, H7, H9)
   - List NA subtypes to analyze (e.g., N1, N2, N9)
   - Set reference accession numbers for each segment-subtype combination
   - Adjust filtering thresholds for sequence curation (max_frac_gaps, max_frac_ambig)
   - Set number of randomizations for tree building (n_randomizations, default: 10)
   - Set desired number of threads
   - Specify host groups to extract for host-specific subtree analysis (host_groups_to_extract)
   - Optionally specify rerooting nodes for final trees (reroot)

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

   ```
   snakemake --cores <number_of_cores> --use-conda
   ```

5. **Output**

   The pipeline organizes outputs by segment and subtype:
   
   **For HA and NA segments** (e.g., `results/HA/H5/` or `results/NA/N1/`):
   - `raw_sequences.fasta.xz`: Parsed sequences for this segment-subtype
   - `reference/`: Reference sequence and Nextclade dataset files
   - `msa.fasta.xz`: Multiple sequence alignment from Nextclade
   - `msa.tsv.xz`: Nextclade alignment metadata
   - `curated_msa.fasta.xz`: Quality-filtered alignment (e.g., gaps < 3%, ambiguities < 0%)
   - `unaligned_coding_seqs/`: Directory containing per-gene unaligned coding sequences extracted from curated sequences
   - `curated_reference.fasta/txt/gff/gtf`: Reference files matching the curated alignment
   - `curated_root.fasta`: Root sequence (reference or inferred from rerooting)
   - `randomized_{n}/`: Directory for each randomized alignment (n = 0, 1, 2, ...)
     - `msa.fasta.xz`: Randomized alignment
     - `msa.vcf`: Variant call format file
     - `preopt_tree.pb.gz`: Initial parsimony tree
     - `opt_tree.pb.gz`: Optimized tree
     - `dag.pb`: DAG representation of the tree
   - `larch_merged_dag.pb`: Merged DAG from all randomizations
   - `trimmed_dag.pb`: Trimmed DAG with suboptimal trees removed
   - `sampled_tree.nh`: Newick format tree sampled from trimmed DAG
   - `sampled_tree.pb.gz`: MAT protobuf of sampled tree
   - `final_tree.pb.gz`: Final tree (rerooted if specified in config)
   - `final_tree.jsonl.gz`: Interactive Taxonium visualization file
   - `host_specific_trees/`: Host-specific subtree visualizations
     - `{host_group}_samples.txt`: Sample list for each host group
     - `{host_group}_tree.pb.gz`: Extracted subtree for each host group
     - `{host_group}_tree.jsonl.gz`: Taxonium visualization for each host group
   
   **For the other segments** (e.g., `results/PB2/all/` or `results/NP/all/`):
   - Same outputs as above, but combining all influenza subtypes

   **Global outputs**:
   - `results/combined_metadata.csv`: Aggregated metadata from all input files
   - `results/combined_metadata_with_host_groups.csv`: Metadata with host group classifications added
   - `results/notebooks/`: Executed analysis notebooks
     - `analyze_metadata.html`: Metadata analysis report
     - `analyze_alignments.html`: Alignment statistics report

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
   - Generates alignment metadata TSV

4. **Curate alignment and extract coding sequences** (`curate_and_extract_coding_seqs.py`):
   - Filters sequences by quality (gaps, ambiguities, terminal gaps)
   - Extracts only coding regions based on GFF annotations
   - Sanitizes sequence IDs for UShER compatibility
   - Filters duplicate sequences
   - Validates CDS for all genes (sequences must pass ALL validations)
   - Extracts per-gene unaligned coding sequences from curated alignment
   - Preserves original nucleotides removed during alignment

5. **Randomize alignment** (`randomize_alignment.py`):
   - Creates multiple randomized versions of the alignment (keeping reference at top)
   - Each randomization uses a different seed for variation in tree building

6. **Create VCF** (faToVcf):
   - Converts each randomized FASTA alignment to variant format
   - Includes reference and handles ambiguous nucleotides

7. **Build initial tree** (usher-sampled):
   - Creates parsimony-based phylogenetic tree for each randomization
   - Uses empty starting tree and builds incrementally

8. **Optimize tree** (matOptimize):
   - Refines tree topology to minimize parsimony score for each randomization
   - Multiple optimization rounds to improve tree quality

9. **Convert to DAG** (larch-usher):
   - Converts each optimized tree to a DAG (Directed Acyclic Graph) representation
   - Allows for representing multiple equally parsimonious tree topologies

10. **Merge DAGs** (larch-dagutil):
    - Combines DAGs from all randomizations into a single merged DAG
    - Trims redundant structures during merge

11. **Trim DAG** (`trim_dag.py`):
    - Removes suboptimal trees from the merged DAG
    - Retains only the most parsimonious tree topologies

12. **Sample tree from DAG** (`convert_DAG_protobuf_to_newick_samples.py`):
    - Samples a representative tree from the trimmed DAG
    - Outputs in Newick format

13. **Create MAT protobuf** (matOptimize):
    - Converts the sampled Newick tree to MAT protobuf format
    - Required for downstream matUtils operations

14. **Reroot tree** (matUtils extract):
    - Reroots the tree at a specified node if configured
    - Otherwise creates a symlink to the sampled tree

15. **Create root sequence** (`extract_root_sequence.py` or symlink):
    - If rerooted: Infers root sequence from tree mutations
    - If not rerooted: Uses reference sequence as root

16. **Add host groups** (`simplified_host_classifier.py`):
    - Classifies hosts into simplified groups
    - Adds host_group column to metadata

17. **Extract host-specific subtrees** (matUtils extract):
    - Creates separate subtrees for each host group
    - Includes all samples from the specified host group plus the root

18. **Create visualizations** (usher_to_taxonium):
    - Converts final tree and host-specific subtrees to Taxonium format
    - Incorporates metadata for interactive exploration

19. **Execute analysis notebooks** (jupyter nbconvert):
    - Runs analysis notebooks after all pipeline outputs are complete
    - Generates HTML reports in `results/notebooks/`
    - Includes metadata analysis and alignment statistics

## Requirements

All dependencies are specified in the environment.yml file.
