# Methods Outline: Flu-UShER Phylogenetic Pipeline

## Data Collection and Preprocessing

- Downloaded influenza virus sequence data and associated metadata from GISAID
- Grouped sequences by genome segment
- For internal segments (PB2, PB1, PA, NP, MP, NS), combined all subtypes into one group for analysis
- For HA and NA segments, we further divided the data by subtype, focusing on the H1, H3, H5, H7, and H9 subtypes for HA and the N1, N2, and N9 subtypes for NA

## Sequence Alignment

- For each group of sequences from above (PB2, HA-H1, HA-H3, etc.), we aligned the sequences to a group-specific reference sequence we downloaded from NCBI (see Table X)
- We did so using nextclade which, for a given group of sequences, iterates over each input sequence and does a codon-aware pairwise alignment to the group's reference, then strips all insertions relative to the reference, and combines the resulting sequences into a single multiple-sequence alignment where all sequences are the same length (including gap characters).

## Quality Control and Sequence Curation

- Applied stringent quality filters to remove low-quality sequences:
  - Excluded sequences with >3% gaps in alignment
  - Excluded sequences with any ambiguous nucleotides (N's)
- Extracted only protein-coding regions based on reference genome annotations (GFF format)
- Removed duplicate sequences
- Filtered out sequences if their coding sequences (unaligned) were not a multiple of three and did not start with a start codon and end with a stop codon.

## Phylogenetic Tree Construction

- Employed a robust tree-building strategy using multiple randomized alignments (n=10) to thoroughly explore tree space
- Converted alignments to variant call format (VCF) relative to reference sequence
- Constructed initial maximum parsimony trees using UShER with incremental sample placement
- Refined tree topologies through iterative optimization (5 rounds) to minimize parsimony score
- Represented phylogenetic uncertainty using Directed Acyclic Graphs (DAGs) to capture multiple equally parsimonious tree topologies
- Merged DAGs from all randomized alignments into consensus representation
- Removed suboptimal tree topologies to retain only most parsimonious solutions
- Sampled representative tree from final trimmed DAG

## Tree Rooting

- For select segment-subtype combinations, rerooted trees at biologically meaningful nodes representing early divergence points
- Inferred ancestral root sequences from mutation patterns along branches leading to root node
- For non-rerooted trees, used reference sequence as root

## Host Classification and Subtree Analysis

- Classified all sequences by host taxonomy into major groups:
  - Avian species
  - Human
  - Swine
  - Other mammalian hosts (bovine, equine, canine, feline, marine mammals)
  - Laboratory/unknown
- Extracted host-specific phylogenetic subtrees to analyze transmission patterns within host groups
- Retained root node in all subtrees to preserve ancestral context

## Data Visualization

- Generated interactive phylogenetic trees in Taxonium format for web-based exploration
- Integrated metadata including:
  - Isolation location and date
  - Host species and host group classification
  - Viral subtype
  - Passage history
  - Genetic clade assignments
- Created separate visualizations for full trees and host-specific subtrees
