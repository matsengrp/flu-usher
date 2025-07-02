# Flu-UShER Pipeline

A Snakemake pipeline for building phylogenetic trees of influenza virus sequences using UShER.

## TODO

* start downloading additional H1N1 data
* workflow:
   * build single tree with sequences from all hosts
   * generate counts from big tree
   * also generate counts from tree subset to sequences from a given host or subtype or whatever
* For HA and NA, build one tree for each subtype regardless of the subtype of other segments
   * H1, H3, H5, H7, H9
   * N1, N2, N9
* For other genes, build one tree with all genes across all subtypes
   * PB2
   * PB1
   * etc.
* Could specify which subtypes to combine for each segment, like:
   * N1: H1N1, H5N1
   * PB2: H1N1, N3N2, etc.
   * could use H3N2 reference when combining across all subtypes
* Get a handle on the overall topologies of the trees

* summary plots
   * summary plots for alignments
      * number downloaded, length distribution, number retained, frac gaps at each site
   * summary plots for final sequences
      * number of sequences per host
      * number of sequences sampled over time
   * summary of trees
      * number of mutations per branch

* using workaround for H1N1 HA

* how identify mutations in test set that go against DASM expectations?
* how identify mutations that go against DMS expectations? Fit selection factor * DMS? How assess significance?
* add tests?
* filter out nonhuman sequences for H3N2?
* allow more mutations per branch? For SARS2 was for whole genome; for flu is only for part; so already more than before

* H5N1 NS: why did the alignment step filter out so many sequences?

* rerooting (genome with sequences for all segments)
* add code for computing counts to pipeline

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
