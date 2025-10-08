"""
Compute consensus sequences for each leaf across multiple UShER trees.
"""

import argparse
import bte
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import lzma
import sys


def extract_leaf_sequences(tree_file, ref_seq):
    """
    Extract sequences for all leaf nodes in a tree.
    
    Args:
        tree_file: Path to the UShER tree (.pb.gz file)
        ref_seq: Reference sequence (SeqRecord object)
    
    Returns:
        Dict mapping leaf IDs to their sequences
    """
    # Get reference sequence string
    ref_seq_str = str(ref_seq.seq)
    
    # Load tree
    tree = bte.MATree(tree_file)
    
    # Get all leaf nodes
    leaf_nodes = tree.get_leaves()
    
    # Extract sequence for each leaf
    sequences = {}
    for leaf in leaf_nodes:
        # Get mutations relative to reference
        muts = tree.get_haplotype(leaf.id)
        
        # Apply mutations to reference sequence
        seq = list(ref_seq_str)
        for mut in muts:
            wt_nt = mut[0]
            site = int(mut[1:-1])
            mut_nt = mut[-1]
            
            # Verify mutation is valid
            if ref_seq_str[site - 1] != wt_nt:
                print(f"Warning: Reference mismatch at site {site} for leaf {leaf.id}: "
                      f"expected {wt_nt}, found {ref_seq_str[site - 1]}", file=sys.stderr)
            
            seq[site - 1] = mut_nt
        
        sequences[leaf.id] = ''.join(seq)
    
    return sequences


def compute_consensus(sequences_list):
    """
    Compute consensus sequence from multiple sequences.
    
    Args:
        sequences_list: List of sequences (strings)
    
    Returns:
        Consensus sequence (string)
    """
    if not sequences_list:
        return ""
    
    seq_length = len(sequences_list[0])
    consensus = []
    
    for pos in range(seq_length):
        # Get all nucleotides at this position
        nucs = [seq[pos] for seq in sequences_list]
        
        # Find most common nucleotide (majority vote)
        counter = Counter(nucs)
        most_common = counter.most_common(1)[0][0]
        consensus.append(most_common)
    
    return ''.join(consensus)


def main():
    parser = argparse.ArgumentParser(description="Compute consensus sequences across multiple UShER trees")
    parser.add_argument("--trees", nargs='+', required=True, 
                        help="Paths to UShER trees (.pb.gz)")
    parser.add_argument("--reference", required=True, 
                        help="Path to reference sequence (FASTA)")
    parser.add_argument("--output", required=True, 
                        help="Output FASTA file with consensus sequences (will be xz-compressed)")
    
    args = parser.parse_args()
    
    # Read reference sequence once
    ref_seq = SeqIO.read(args.reference, "fasta")
    ref_id = ref_seq.id
    print(f"Reference sequence ID: {ref_id}", file=sys.stderr)
    
    # Extract sequences from all trees
    all_tree_sequences = []
    all_leaf_ids = None
    
    for i, tree_file in enumerate(args.trees):
        print(f"Processing tree {i+1}/{len(args.trees)}: {tree_file}...", file=sys.stderr)
        sequences = extract_leaf_sequences(tree_file, ref_seq)
        
        # Get leaf IDs and ensure consistency across trees
        leaf_ids = set(sequences.keys())
        if all_leaf_ids is None:
            all_leaf_ids = leaf_ids
            print(f"  Found {len(all_leaf_ids)} leaves", file=sys.stderr)
        else:
            if leaf_ids != all_leaf_ids:
                print(f"ERROR: Tree {i+1} has different leaf IDs than previous trees!", file=sys.stderr)
                print(f"  Expected {len(all_leaf_ids)} leaves, found {len(leaf_ids)}", file=sys.stderr)
                
                # Show differences
                missing = all_leaf_ids - leaf_ids
                extra = leaf_ids - all_leaf_ids
                if missing:
                    print(f"  Missing leaves (first 10): {list(missing)[:10]}", file=sys.stderr)
                if extra:
                    print(f"  Extra leaves (first 10): {list(extra)[:10]}", file=sys.stderr)
                
                sys.exit(1)
        
        all_tree_sequences.append(sequences)
    
    # Compute consensus sequences for each leaf
    print(f"\nComputing consensus sequences for {len(all_leaf_ids)} leaves...", file=sys.stderr)
    consensus_sequences = {}
    
    for leaf_id in sorted(all_leaf_ids):
        # Collect sequences for this leaf from all trees
        leaf_sequences = [tree_seqs[leaf_id] for tree_seqs in all_tree_sequences]
        
        # Compute consensus
        consensus = compute_consensus(leaf_sequences)
        consensus_sequences[leaf_id] = consensus
    
    # Verify we have the expected number of output sequences
    assert len(consensus_sequences) == len(all_leaf_ids), \
        f"Output sequence count mismatch: {len(consensus_sequences)} != {len(all_leaf_ids)}"
    
    # Check that reference sequence exists in consensus sequences
    if ref_id not in consensus_sequences:
        print(f"ERROR: Reference sequence '{ref_id}' not found in consensus sequences!", file=sys.stderr)
        print(f"Available sequence IDs (first 10): {sorted(list(consensus_sequences.keys()))[:10]}", file=sys.stderr)
        sys.exit(1)
    
    # Write to output FASTA with reference sequence first
    output_file = args.output if args.output.endswith('.xz') else f"{args.output}.xz"
    print(f"Writing {len(consensus_sequences)} consensus sequences to {output_file}...", file=sys.stderr)
    records = []
    
    # Add reference sequence first
    record = SeqRecord(Seq(consensus_sequences[ref_id]), id=ref_id, description="")
    records.append(record)
    
    # Add remaining sequences in sorted order (excluding reference)
    for seq_id in sorted(consensus_sequences.keys()):
        if seq_id != ref_id:
            record = SeqRecord(Seq(consensus_sequences[seq_id]), id=seq_id, description="")
            records.append(record)
    
    # Write to xz-compressed file
    with lzma.open(output_file, 'wt') as handle:
        SeqIO.write(records, handle, "fasta")
    
    # Verify we wrote the correct number of sequences
    assert len(records) == len(consensus_sequences), \
        f"Output record count mismatch: {len(records)} records != {len(consensus_sequences)} consensus sequences"
    
    print(f"Done! Written {len(records)} sequences to {output_file} (reference '{ref_id}' first).", file=sys.stderr)


if __name__ == "__main__":
    main()