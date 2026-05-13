"""
Drop terminals whose branch length exceeds --max-branch-length from a newick
tree, and write the filtered tree.

Motivation: chronumental's SVI fitting can place extreme date predictions on
internal nodes when the tree contains terminals with very long branches that
lack date anchors (e.g. synthetic reverse-genetics reassortants like
A/Puerto Rico/8-RGcH5-1/1934). Those predictions can exceed pandas' ~292-year
timedelta limit and break the dates-TSV writer. Pruning the long-branch tips
before dating avoids the issue without affecting the topology of the rest of
the tree.
"""
import argparse
from Bio import Phylo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input newick tree")
    parser.add_argument("--output", required=True, help="Output filtered newick tree")
    parser.add_argument("--dropped", required=True,
                        help="Output TSV listing dropped tips (columns: name, branch_length)")
    parser.add_argument("--max-branch-length", type=float, required=True,
                        help="Drop terminals whose branch length exceeds this value")
    args = parser.parse_args()

    tree = Phylo.read(args.input, "newick")
    n_before = sum(1 for _ in tree.get_terminals())

    long_tips = [t for t in tree.get_terminals()
                 if (t.branch_length or 0) > args.max_branch_length]
    long_tips.sort(key=lambda t: -(t.branch_length or 0))

    with open(args.dropped, "w") as f:
        f.write("name\tbranch_length\n")
        for t in long_tips:
            f.write(f"{t.name}\t{t.branch_length}\n")

    print(f"Tree has {n_before} terminals; "
          f"{len(long_tips)} have branch_length > {args.max_branch_length}")
    for t in long_tips:
        print(f"  dropping {t.name}\tbranch_length={t.branch_length}")

    for t in long_tips:
        tree.prune(t)

    n_after = sum(1 for _ in tree.get_terminals())
    Phylo.write(tree, args.output, "newick")
    print(f"Wrote filtered tree with {n_after} terminals to {args.output}")


if __name__ == "__main__":
    main()
