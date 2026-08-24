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
│   ├── create_scaffold_alignment.py                # Subset an alignment to a time-spread sample
│   ├── trim_dag.py                                 # Trim suboptimal trees from DAG
│   ├── convert_DAG_protobuf_to_newick_samples.py   # Sample tree from DAG
│   ├── create_root_samples_file.py                 # Create samples file for root extraction
│   ├── extract_root_sequence.py                    # Infer root sequence from tree
│   ├── simplified_host_classifier.py               # Classify hosts into groups (used by augment_metadata.py)
│   ├── augment_metadata.py                         # Add host and geographic group columns to metadata
│   ├── create_samples_file.py                      # Create samples file for subtree extraction (generic column filter)
│   ├── prepare_host_annotation.py                  # Build 2-col CSV (isolate_id, host_group) for PastML
│   └── prepare_subtype_annotation.py               # Build 2-col CSV (isolate_id, subtype) for PastML, normalizing "A / H5N1" → "H5N1"
└── notebooks/               # Jupyter notebooks for development and analysis
```

## Usage

1. **Set up your environment**

   Create and activate the conda environment:
   ```
   conda env create -f environment.yml
   conda activate flu-usher
   ```

   Clone [larch](https://github.com/matsengrp/larch), then build and install it from source into the environment, where `$CONDA_PREFIX` is the path to the active conda environment (e.g., `/home/hhaddox/miniforge3/envs/flu-usher/`):
   ```
   git clone https://github.com/matsengrp/larch.git
   cd larch
   git submodule update --init --recursive
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
   make -j16
   make install
   ```

   Other dependencies are managed per-step via conda environments in `envs/`, which Snakemake invokes automatically.

2. **Configure the pipeline**

   Edit `config.yaml` to:
   - Specify input directories containing your GISAID data
   - List HA subtypes to analyze (e.g., H1, H3, H5, H7, H9)
   - List NA subtypes to analyze (e.g., N1, N2, N9)
   - Set reference accession numbers for each segment-subtype combination
   - Adjust filtering thresholds for sequence curation (max_frac_gaps, max_frac_ambig)
   - Set number of randomizations for tree building (n_randomizations, default: 10)
   - Set desired number of threads
   - Specify geographic groups to extract (geographic_groups_to_extract)
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
     - `scaffold_msa.fasta.xz`: Time-spread subset of the alignment (see [Design Notes](#design-notes))
     - `scaffold_samples.txt`: Ids that went into the scaffold, one per line
     - `scaffold_tree.pb.gz`: Optimized backbone tree built from that subset
     - `preopt_tree.pb.gz`: Initial parsimony tree
     - `opt_tree.pb.gz`: Optimized tree
     - `dag.pb`: DAG representation of the tree
   - `larch_merged_dag.pb`: Merged DAG from all randomizations
   - `trimmed_dag.pb`: Trimmed DAG with suboptimal trees removed
   - `sampled_tree.nh`: Newick format tree sampled from trimmed DAG
   - `sampled_tree.pb.gz`: MAT protobuf of sampled tree
   - `final_tree.pb.gz`: Final tree (rerooted if specified in config)
   - `final_tree.jsonl.gz`: Interactive Taxonium visualization file
   - `geographic_trees/`: Geographic subtree visualizations
     - `{geo_group}_samples.txt`: Sample list for each geographic group
     - `{geo_group}_tree.pb.gz`: Extracted subtree for each geographic group
     - `{geo_group}_tree.jsonl.gz`: Taxonium visualization for each geographic group
   - `host_ancestral/`: Per-node host inference (PastML / DOWNPASS)
     - `combined_ancestral_states.tab`: Tab-separated file of inferred `host_group` for every node (leaves + internals); the `node` column joins to internal node IDs in `final_tree.pb.gz`. Ambiguous internals may appear on multiple rows (one per equally-parsimonious state).
     - `host_tree.html`: Interactive PastML visualization of the ancestral reconstruction.
     - `named.tree_final_tree.nwk`: PastML's labeled-tree output (identical names to `final_tree.nwk` here).
   - `subtype_ancestral/`: Per-node subtype inference (PastML / DOWNPASS)
     - `combined_ancestral_states.tab`: Tab-separated file of inferred `subtype` (`H*N*` form, normalized from the raw GISAID `subtype` column) for every node (leaves + internals); same join + ambiguity semantics as `host_ancestral/`. On HA per-subtype trees the H letter is fixed across all tips and only the inferred N partner varies (analogously for NA trees); on internal-segment trees neither letter is constrained, so the full `H*N*` can vary along the tree.
     - `subtype_tree.html`: Interactive PastML visualization of the ancestral reconstruction.
   
   **For the other segments** (e.g., `results/PB2/all/` or `results/NP/all/`):
   - Same outputs as above, but combining all influenza subtypes

   **Global outputs**:
   - `results/combined_metadata.csv`: Aggregated metadata from all input files
   - `results/combined_metadata_augmented.csv`: Metadata with host_group and geographic_group columns added
   - `results/host_annotation.csv`: 2-column (isolate_id, host_group) annotation table for PastML
   - `results/subtype_annotation.csv`: 2-column (isolate_id, subtype) annotation table for PastML (subtype normalized to `H*N*` from the raw `subtype` column)
   - `results/input_data_md5sums.txt`: Provenance manifest with md5 hashes for every input FASTA / XLS file under the configured `input_dirs`, one file per line in standard `md5sum` format (`<hash>  <path>`)
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

7. **Build initial tree** (`create_scaffold_alignment.py`, faToVcf, usher-sampled, matOptimize):
   - Draws `scaffold_n_taxa` sequences spread evenly across collection year into a scaffold alignment
   - Builds and optimizes a backbone tree from that subset alone
   - Places every remaining sequence onto the backbone with usher-sampled
   - Why the backbone exists, and how the draw works: see [Design Notes](#design-notes)

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

13. **Create MAT protobuf** (usher):
    - Converts the sampled Newick tree to MAT protobuf format
    - Required for downstream matUtils operations

14. **Reroot tree** (`reroot_newick.py`, usher, `rebase_mat_root.py`, `check_tree_sequences.py`):
    - Where `reroot` is configured, reroots the *newick* at the named leaf, then rebuilds the MAT from the alignment; other combinations symlink the sampled MAT through unchanged
    - Moves the MAT's origin from the reference onto the tree's own root, writing that root out as `final_tree_root.fasta`
    - Checks that no sample's sequence changed across rerooting and that the tree is rooted where the config asked, writing `reroot_sequence_check.txt`
    - Why rerooting happens in newick space rather than with `matUtils extract -y`: see [Design Notes](#design-notes)

15. **Create root sequence** (`extract_root_sequence.py` or symlink):
    - If rerooted: Infers root sequence from tree mutations
    - If not rerooted: Uses reference sequence as root

16. **Augment metadata** (`augment_metadata.py`):
    - Adds host_group column (human, avian, swine, bovine, other)
    - Adds geographic_group column (north_america, europe, asia, other)

17. **Extract geographic subtrees** (matUtils extract):
    - Creates separate subtrees for each configured geographic region (filter on the augmented `geographic_group` column)
    - Each subtree includes the root sequence plus matching samples

18. **Infer per-node host states** (PastML, `prepare_host_annotation.py`):
    - Builds a 2-column annotation table from `combined_metadata_augmented.csv` mapping `isolate_id → host_group`
    - Runs PastML with the DOWNPASS parsimony method on `final_tree.nwk` to reconstruct `host_group` for every internal node
    - **Inputs:** `final_tree.nwk`, `combined_metadata_augmented.csv`
    - **Outputs:** `host_ancestral/combined_ancestral_states.tab` (per-node host_group; node IDs match `final_tree.pb.gz`), `host_ancestral/host_tree.html` (interactive viz), `host_ancestral/named.tree_final_tree.nwk`
    - Ambiguous internals appear on multiple rows of `combined_ancestral_states.tab` (one per equally-parsimonious state)

19. **Infer per-node subtype states** (PastML, `prepare_subtype_annotation.py`):
    - Builds a 2-column annotation table from `combined_metadata_augmented.csv` mapping `isolate_id → subtype`, normalizing the raw GISAID `subtype` ("A / H5N1") to `H*N*` form ("H5N1") inside the prep script (`unknown` for non-matches)
    - Runs PastML with DOWNPASS on `final_tree.nwk` to reconstruct `subtype` for every internal node
    - **Inputs:** `final_tree.nwk`, `combined_metadata_augmented.csv`
    - **Outputs:** `subtype_ancestral/combined_ancestral_states.tab`, `subtype_ancestral/subtype_tree.html`
    - On HA per-subtype trees the H letter is fixed across all tips and only the inferred N partner varies (analogously for NA trees); on internal-segment trees (`PB2/all`, `PB1/all`, `PA/all`, `NP/all`, `MP/all`, `NS/all`) neither letter is constrained, so the full `H*N*` can vary along the tree and ancestral subtype tracks reassortment events

20. **Create visualizations** (usher_to_taxonium):
    - Converts final tree and all subtrees to Taxonium format
    - Incorporates metadata (including host and geographic groups, and the raw `collection_date`) for interactive exploration

21. **Execute analysis notebooks** (jupyter nbconvert):
    - Runs analysis notebooks after all pipeline outputs are complete
    - Generates HTML reports in `results/notebooks/`
    - Includes metadata analysis and alignment statistics

22. **Record input data md5 sums** (md5sum):
    - Walks every `.fasta` and `.xls` file under the configured `input_dirs` and writes their md5 hashes to `results/input_data_md5sums.txt` as a provenance manifest
    - Files are listed as explicit Snakemake inputs, so the manifest is regenerated whenever any input data file changes

## Design Notes

Why the pipeline is built the way it is. The Snakefile carries only short
comments saying what each rule does; the reasoning lives here, so that these
notes can be read as a whole and can carry evidence at length.

### Why each tree search starts from a time-spread scaffold

Issue #53. Placing all 86,232 HA/H3 sequences onto an empty tree settled on a
high-level shape that is measurably not the most parsimonious one available:
sequences collected from 2006 on sat ~70 branches deeper on the trunk than
2005's, with no matching rise in divergence. Building a ~1000-sequence backbone
first and placing everything else onto it fixed that *and* improved the score,
which is what says the search rather than the data was at fault. On the same
86,232 sequences, HA/H3's `final_tree` went from 185,398 to 185,301 parsimony
and 2006's median depth from 161 to 93.

The clearest evidence is per-randomization rather than on the merged result.
Rerooted and measured by cohort, 4 of the 10 baseline `opt_tree`s had 2006
sitting 18-28 branches *deeper* than 2013; with a scaffold, 0 of 10 do, and 9 of
10 are clean at 2/2571 leaves inverted. So the per-randomization search was
itself part of the problem — an earlier write-up claimed the search was never at
fault, which was wrong, and rested on the two randomizations that happened to be
clean.

Two measurement traps, both of which earlier write-ups fell into:

- **Do not compare a scaffolded subset against the same leaves extracted from a
  full tree.** Restricting an 86k-leaf tree to 1k leaves inflates its parsimony
  there, because a globally optimal tree does not induce an optimal subtree. The
  scaffolded tree is globally better (185,301 vs 185,398) and still scores 8365
  on such a subset against a de novo 8298.
- **Do not quote the share of pre-2009 leaves deeper than the median 2013 leaf
  on its own.** It reads 35% for the pathological baseline tree while the
  aggregate median gap simultaneously reads -46, because pre-2009 is dominated
  by shallow pre-2005 leaves. Use the 2006 / 2007 cohort medians.

### How the scaffold is sampled

Quotas are per collection *year*, not per sequence: HA/H3 has eight sequences
from 1963 and 11,786 from 2024, and a plain random draw of 1000 would take ~150
from 2024 alone and miss most years before 1991 entirely. Years that cannot fill
their share hand the remainder back to be redistributed — repeatedly, since a
redistribution can itself overfill a small year. Sequences with no collection
date are not candidates, since they carry no information about where on the time
axis they belong; they are placed later, in the pass that adds every
non-scaffold sequence, so nothing is dropped.

Each draw uses that randomization's seed, so the backbones differ: on HA/H3 any
two share ~19% of their ids, which is what lets the merged DAG keep
exploring topology space. Note the draw depends on the seed alone and *not* on
the alignment's order — the sort in `select_scaffold_ids` is deliberate — so the
shuffle decides only the order the chosen records are written in. On small
combinations the backbones are necessarily more alike: NA/N9 draws 1000 of just
2048 dated sequences, and any two of its scaffolds share ~49%.

A per-year quota on its own is not enough to keep them apart, which is what
`scaffold_max_year_fraction` is for. A year holding 21 sequences against a
quota of 17 hands over almost the same 17 under every seed: it spends its slots
and buys no diversity. The parameter caps what any one year may give at that
fraction of what it holds, floored at one sequence so no year is dropped for
being small — `int()` alone would delete HA/H7's 1902 and NA/N9's 1968, single
sequences that are exactly what the scaffold exists to anchor. Mean pairwise
overlap between two backbones, at 0.5 versus unbounded, measured over the 10
randomizations: HA/H3 19% vs 38%, NA/N8 30% vs 49%, HA/H7 41% vs 58%, NA/N9 49%
vs 78%, with the distinct-sequence union rising in every case (HA/H3 5366 →
6209 over ten draws).

What the bound is really reclaiming is wasted slots, and this is the clearest
way to see why it is not a trade against historical coverage. Unbounded, on
HA/H3, 255 of each backbone's 1000 draws come from years holding no more
sequences than their quota — every one of them 1990 or earlier. Those years have
nothing to vary, so all ten backbones contain the *identical* early set, and
nine of the ten copies tell the merged DAG nothing the first had not already
given it. A quarter of every backbone is spent re-deriving the same sequences
ten times. (This is also why unbounded whole-backbone overlap reads so high:
across the years that are genuinely sampled it is 17%.)

Bounding it does not cost that coverage, because what reaches the merged DAG is
the union of the ten draws, not any one of them. Each backbone takes about half
of a small year, but a different half, and a sequence in a year drawn at 0.5 is
missed by all ten draws only with probability of order 2⁻¹⁰. So the collection
still sees essentially every early sequence while the freed slots go to varied
recent ones, and the union rises rather than falls: NA/N9 from 1885 to 2046 of
its 2048 candidates, HA/H3 from 5366 to 6209. What genuinely changes is
per-backbone early-taxon density — each individual backbone resolves its trunk
from fewer early taxa — and the floor of one keeps every year present in every
backbone regardless.

The constraint to respect when changing it is capacity, not coverage: at 0.5 a
combination can field about half its dated sequences, so NA/N9 offers 1016
against the 1000 requested — a margin of 16, the tightest of the 16
combinations, where every other has at least 1599. Below about 0.5 the small
combinations cannot fill the draw at all (at 0.25 NA/N9 manages 509 and HA/H7
800), which `select_scaffold_ids` raises on rather than shipping a short
scaffold and letting `scaffold_n_taxa` become a fiction.

`create_scaffold_alignment` reads `combined_metadata.csv` rather than the
augmented file. It needs only `collection_date`, and depending on the augmented
file would put host and geography classification upstream of every tree in the
pipeline.

`scaffold.vcf` is `temp()` while the full `msa.vcf` is not: it has a single
consumer, and faToVcf regenerates it in seconds from the scaffold alignment,
which is kept precisely so that it can.

`build_scaffold_tree` runs usher-sampled and matOptimize in one rule, because
both take seconds at this size and there is nothing to gain from a checkpoint
between them. Their full-size equivalents (`create_initial_tree`,
`optimize_tree`) stay separate because they are the two most expensive steps in
the pipeline.

`create_initial_tree` takes the scaffold as a MAT (`-i`) rather than as a newick
(`-t`), because `matUtils extract -t` segfaults in this usher build even on a
~1600-node tree (reproduced, core dumped), so the routine way of handing the
optimized scaffold over as a newick is unavailable. Feeding the MAT straight
through costs nothing either way, since usher re-derives the backbone's states
from the full VCF. Samples already in the scaffold are skipped rather than
duplicated: HA/H3's `preopt_tree` holds exactly 86,232 samples and
`placement_stats.tsv` has 85,231 rows, i.e. 86,232 - 1,001.

### Collection dates

GISAID supplies a `Collection_Date` for every isolate in the corpus as it
stands today, at one of three precisions: 582,016 as `YYYY-MM-DD`, 15,488 as `YYYY-MM` and 31,602 as `YYYY`,
over the 629,106 isolates in the 21 configured input dirs. `parse_gisaid_data`
writes the string through to `combined_metadata.csv` exactly as given, and
consumers that need a full date parse it themselves.

It did not always. The column used to go through
`pd.to_datetime(df['collection_date'], errors='coerce')`, and pandas >= 2.0
infers a *single* format from the column's first non-null element, then coerces
everything not matching it to `NaT`. The call sat inside the per-file loop, so
one row decided the fate of that file's whole date column. Five of the 56
`data/*/*.xls` files open on a partial date, and in exactly those five every
fully-resolved date was destroyed: `H3N2-window-8.xls` (`2022-11`) 19,811,
`H1N1-window-1.xls` (`2009`) 17,908, `H3NX` (`2001`) 4,653, `H8NX` (`2006`) 237
and `H15NX` (`1983`) 10 — 42,619 valid dates, silently, because `errors='coerce'`
has nothing to complain about. `combined_metadata.csv` came out with 87,854 of
629,106 dates blank: those 42,619 plus the 45,235 genuine partials the same call
discarded. Issue #55.

The narrow fix is `format='mixed'`, which recovers all 42,619. It is the wrong
one: it also maps `2009` to January 1 and `2022-11` to the 1st, inventing a day
for 47,090 rows. The bug already did this to the 1,855 partials that happened to
match their file's inferred format — `EPI_ISL_65694` is `2009` in the `.xls` and
`2009-01-01` in today's CSV — and normalising the rest would grow the fabrication
by 25x while destroying the one thing a consumer needs in order to drop partial
dates: the ability to tell them apart. Hence the raw string. There is no
precision column either, because precision is the string's length and no
consumer needs it named.

The one in-pipeline consumer already reads the column as text.
`create_scaffold_alignment` takes `date[:4]`, so partials are usable year
buckets and the fix strictly enlarges its candidate pool — every scaffold
candidate in all 16 combinations is now dated, and the tightest combination,
NA/N9, goes from 1869 dated of 2048 to all 2048, so the `scaffold_n_taxa`
headroom described in `config.yaml` strictly increases and no combination moves
toward usher-sampled's "No samples to place". Everything downstream of that only
displays the column: `usher_to_taxonium` passes `collection_date` into the tree
view at whatever precision GISAID supplied. `augment_metadata` was the second
consumer until issue #58 — it split the dates into `early`/`late` halves against
a global median, on a strict `%Y-%m-%d` parse that left all 47,090 partials
`unknown` — but that column was display-only as well, and the raw date shows
the same thing at full resolution.

What the raw string buys is that a consumer needing exact dates can reject
partials instead of being handed invented ones. What it costs is that such a
consumer must now do so deliberately: a bare `pd.to_datetime(..., errors='coerce')`
downstream of this file will infer `%Y` from a leading year and quietly `NaT`
every full date — the same trap, one repo over. Parse with `format='ISO8601'`,
never `'mixed'`: on `['03/04/2020', '2020-05-06']`, `ISO8601` returns `NaT` for
the ambiguous value while `'mixed'` silently guesses March 4.

`parse_gisaid_data` now fails the run on any date matching none of the three
precisions, listing the file and up to five offending isolate ids. Blanks are
logged but tolerated, since a blank is missing data rather than a change in what
GISAID ships. All 56 files pass today. The strict `%Y-%m-%d` consumer that makes
this check worth having is no longer in this repo: it is
`scripts/date_filters.py:parse_iso_date` in `flu-dasm-antigenic-evo`, which
returns `None` on a date it cannot parse and so drops the sequence from
chronumental's node dating and the root-to-tip regression without a word. That
is the "one repo over" named above, and it is why an impossible date has to die
here rather than downstream.

### Rerooting: newick space, not `matUtils extract -y`

Issue #49. `matUtils extract -y` rebases the MAT into the new root's coordinate
frame: the new root is written with zero mutations, so the file stops recording
how it differs from the reference. That is self-consistent and readable if you
know the new root's sequence, but it refuses outright when the existing root
carries a mutation (which broke HA/H3), and the frame shift is implicit —
nothing records it. Rerooting the topology in newick space instead leaves
mutation assignment to the MAT-building rule, which reads it off the alignment
and so keeps the MAT in the reference's frame.

`rebase_final_mat` then performs that frame shift deliberately, for all 16
combinations alike so that `final_tree.pb.gz` means the same thing everywhere.
It is pure bookkeeping: every branch below the root already records a
parent-to-child change independent of what the origin is, so only the root's own
mutation list (emptied) and each mutation's `ref_nuc` annotation (repointed)
change. Unlike `-y`, it does not refuse when the root carries mutations, and it
emits the root sequence rather than discarding it — that discard is what made
the shift silent in #49.

### Why `usher` and not `matOptimize` builds the final MATs

Both `create_mat_protobuf` and `create_final_mat` annotate a *fixed* topology;
they do not search for a better one. `matOptimize` is a parsimony optimizer, and
parsimony is an *unrooted* criterion, so it normalizes the tree at load —
collapsing zero-mutation branches — even under `-N 0`, and does not preserve the
input rooting.

That silently broke NA/N1. Its reroot target `EPI_ISL_5878` is the only one of
the 14 whose terminal branch carries no mutations (the others run 6-87 in
`sampled_tree.pb.gz`, 5-62 after rerooting), and it was the only one that lost
its rooting: the root landed 54 mutations away under a 3-way polytomy, 56
substitutions from the intended root sequence — not a simple collapse into its
neighbour. Empirically the rooting survives `matOptimize` only where the
target's terminal branch carries mutations; the exact normalization step
responsible has not been traced to matOptimize's source. Switching both rules to
`usher`, which takes the topology as authoritative, fixed it and moved
`curated_root.fasta` by one position in 2 of 16 combinations. `matOptimize`
remains correct for `optimize_tree`, which really is topology optimization.

### The sequence-identity gate

`check_tree_sequences` is a standing guard against the class of bug in #49:
rerooting is a pure re-representation, so every sequence the tree implies must
survive it unchanged. It runs for every combination, including the two that are
not rerooted, where it is a cheap tautology.

It also asserts that the tree is rooted where the config asked, because the
sequence comparison alone cannot tell: both trees are built from the same VCF
and so agree for any topology.

Everything derived from the final tree — the Taxonium conversions, the
geographic subtrees, `final_tree.nwk` and the two PastML rules that read it —
depends on the gate's report rather than merely being listed alongside it in
`rule all`. That dependency is what makes it a gate. Building one combination by
targeting its `final_tree.jsonl.gz` produces a DAG that never reaches `rule
all`, so without those edges the check would be skipped for exactly the workflow
people use most and a bad tree would be published anyway.

### larch is built from source, not installed from conda

`tree_to_dag` and `larch_merge` have no `conda:` directive, and `snakemake
--lint` flags both for it. That is accepted, not an oversight: larch is built
from source into the main flu-usher environment (`environment.yml` carries its
build and runtime deps, not larch itself), from commit `0ac4146`, built
2025-10-27.

`envs/larch.yaml` used to point at the packaged `larch-phylo`, and was deleted
rather than wired up. The package's "0.1.0" version string is a stale
VERSION-file fallback — conda-build strips `.git`, so its `git describe` never
runs — not a lower release: it is built from `c2e75a2`, which is the `v0.1.3`
tag. The real gap is the 5 commits from there to `0ac4146`, and adopting it
would not have degraded anything quietly, it would have broken outright:
`fb1bb78` is what added VCF support to `larch-dagutil`, so the packaged binary
has no `-v` option and `larch_merge` passes one. It also pins `python_abi` 3.8
and `boost-cpp` 1.76. Packaging larch properly is a separate change that has to
be validated against the trees it produces.

### Snakemake bookkeeping

**`script_deps`.** Only rules using Snakemake's `script:` directive are
code-tracked. The scripts invoked as `shell: python scripts/x.py` appear in no
`input:` block, so editing one leaves Snakemake reporting every downstream
output up to date. Declaring them as inputs closes that gap. Every script
imports `utils`, so it is always included.

**Config invalidation.** `download_all_references` deliberately does *not* take
`config.yaml` as an input. It sits at the head of the DAG, so depending on the
file's mtime made every config edit — a reroot target, a seed, a filtering
threshold — re-download the references and re-derive the whole pipeline;
measured before the change, one altered reroot target re-ran 1235 of 1235 jobs.
Naming the consumed values in `params:` instead keeps the default `params`
rerun-trigger precise, and the script still reads the file for the rest.

**Wildcard constraints.** Without them the default wildcard regex is `.+`, which
matches `/`. A typo in a target path then resolves against the wrong rule with a
wildcard spanning directory separators, instead of failing.

**Empty input list.** `input_data_md5sums` fails at parse time when no GISAID
input files are found, rather than degrading silently. With an empty list,
`md5sum {input} > {output}` becomes a bare `md5sum` reading stdin, which writes
`d41d8cd98f00b204e9800998ecf8427e  -` and exits 0 — a provenance manifest
recording nothing. `set -euo pipefail` would not catch it, because
zero-argument `md5sum` succeeds.

**Thread declarations.** usher detects hardware concurrency when `-T` is
omitted, so a rule that did not declare `threads` claimed 1 core from the
scheduler and then spawned 32. Both MAT-building rules declare it and pass it
through.

**Private scratch directories.** `usher` and `usher-sampled` write
`final-tree.nh` and other fixed-name files beside their `-d` target, so rules
that would otherwise share a directory each get a `mktemp -d`. `matUtils`
prepends `-d` to the `-S` path even when it is absolute and then silently writes
nothing, which is why `-S` must stay a bare filename.

**One rule per notebook.** The notebooks produce figures unevenly —
`analyze_alignments` 3, `analyze_dags` 1, `analyze_metadata` 0 — so a single
`{notebook}` wildcard rule cannot declare them, because an output list cannot
depend on a wildcard value. Each writes its executed copy under `results/`
instead of running `--inplace` over the git-tracked source, which rewrote
tracked files and churned output-cell diffs on every run while declaring only a
`.done` sentinel.

## Requirements

The main conda environment (`environment.yml`) provides Snakemake and the build dependencies needed to compile larch. Per-step dependencies are managed via separate conda environments in `envs/`, which Snakemake invokes automatically using the `--use-conda` flag.
