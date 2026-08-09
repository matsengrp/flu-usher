"""
Reroot a newick tree at a named leaf, using it as an outgroup.

Replaces `matUtils extract -y`, which reroots by moving the MAT's origin onto
the new root: that root is written with zero mutations, so the file no longer
records how it differs from the reference. The result is self-consistent --
every branch's parent-to-child mutations are correct, and it reads correctly if
you know the new root's sequence -- but it carries two costs. It refuses
outright when the pre-reroot root already carries mutations, which is what
broke HA/H3, and the move of origin is implicit: nothing in the artifact
records it, so any consumer pairing final_tree.pb.gz with
curated_reference.fasta silently misreads it. See issue #49.

Rerooting the newick instead leaves mutation assignment to `usher`, which reads
it off the alignment, so the MAT stays recorded against the reference and needs
no reinterpretation. (Not matOptimize: create_final_mat must use usher, because
matOptimize would collapse the very branch this rerooting depends on -- see the
comment on that rule.) matUtils' own log confirms the topology semantics match:
"New root was a leaf node; retaining it as leaf node on new root internal node"
is what set_outgroup does with a leaf.
"""
import argparse
import sys

from ete3 import Tree

from utils import setup_logging

logger = setup_logging(__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input newick file")
    parser.add_argument("--output", required=True, help="Output newick file")
    parser.add_argument("--root", required=True,
                        help="Leaf name to reroot at, used as the outgroup")
    args = parser.parse_args()

    # format=1 keeps the internal node names historydag writes (internal_1234).
    tree = Tree(args.input, format=1)

    # `tree & name` matches internal nodes too and raises a bare TreeError on a
    # miss, so check explicitly rather than letting either case through quietly.
    matches = tree.search_nodes(name=args.root)
    if not matches:
        logger.error(f"reroot target '{args.root}' is not present in {args.input}")
        sys.exit(1)
    if len(matches) > 1:
        logger.error(
            f"reroot target '{args.root}' appears {len(matches)} times in "
            f"{args.input}; it must be unique"
        )
        sys.exit(1)

    node = matches[0]
    if not node.is_leaf():
        logger.error(
            f"reroot target '{args.root}' is an internal node in {args.input}; "
            "only leaves are supported"
        )
        sys.exit(1)

    n_leaves = len(tree)
    tree.set_outgroup(node)

    # set_outgroup rewrites the root, so a lost or duplicated leaf here means a
    # silently truncated tree downstream.
    if len(tree) != n_leaves:
        logger.error(
            f"rerooting changed the leaf count from {n_leaves} to {len(tree)}"
        )
        sys.exit(1)

    # The point of the rule. A no-op set_outgroup would leave every downstream
    # check happy -- the sequence gate passes whether or not rerooting happened,
    # because both trees are matOptimize products of the same alignment -- so
    # this is the only place that asserts the tree is rooted where it was asked.
    children = [child.name for child in tree.children]
    if args.root not in children:
        logger.error(
            f"reroot target '{args.root}' is not a child of the new root; "
            f"root children are {children}"
        )
        sys.exit(1)

    # sampled_tree.nh carries no branch lengths; ete3 invents dist=1.0 and writes
    # them out. Inert -- usher recomputes lengths as mutation counts from the
    # VCF, and the current pipeline already feeds it a length-free newick.
    tree.write(outfile=args.output, format=1)
    logger.info(
        f"Rerooted {n_leaves} leaves at '{args.root}' -> {args.output}"
    )


if __name__ == "__main__":
    main()
