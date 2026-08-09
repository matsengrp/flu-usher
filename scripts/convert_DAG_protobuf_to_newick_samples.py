"""
Sample a single history out of a trimmed DAG and write it as a newick tree.

Reproducibility needs two things, and neither alone is sufficient:

1. `random.seed(...)`. fast_sample() draws edges with random.choices() from
   Python's global `random`; historydag documents seeding via random.seed.
2. A fixed PYTHONHASHSEED, set by the create_newick rule. Clades are keyed on
   frozensets of string labels, so the traversal order fast_sample() walks
   depends on string hash randomisation. With the seed alone, four runs on one
   byte-identical trimmed_dag.pb still produced four different trees.
"""

import argparse
import random

import historydag as hdag


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sample a history from a trimmed DAG and write it as newick"
    )
    parser.add_argument("--input", required=True, help="Trimmed DAG protobuf")
    parser.add_argument("--output", required=True, help="Output newick file")
    parser.add_argument(
        "--seed", type=int, required=True, help="Seed for the history sample"
    )
    return parser.parse_args()


def name_internal_nodes_distinctly(n):
    """
    Name a DAG node for the newick output.

    Internal node ids are prefixed so they cannot collide with leaf names,
    which some downstream tools treat as the same taxon. The UA node is
    unnamed.
    """
    if n.is_leaf():
        return n.label.node_id
    elif not n.is_ua_node():
        return "internal_" + str(n.label.node_id)
    return ""


def main():
    args = parse_args()
    random.seed(args.seed)

    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(args.input)
    newick = dag.fast_sample().to_newick(
        name_func=name_internal_nodes_distinctly, features=[], feature_funcs={}
    )
    with open(args.output, "w") as f:
        f.write(newick)


if __name__ == "__main__":
    main()
