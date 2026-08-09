"""
Trim a merged history DAG down to its optimal-weight histories.

Invoked with PYTHONHASHSEED fixed (see the trim_dag rule). historydag keys its
clades on frozensets of string labels, so set iteration order -- and therefore
the order nodes are written to the output protobuf -- depends on Python's string
hash randomisation. Without a fixed PYTHONHASHSEED the same input DAG produces a
different trimmed_dag.pb on every run.
"""

import argparse

import historydag as hdag


def parse_args():
    parser = argparse.ArgumentParser(
        description="Trim a merged DAG to its optimal-weight histories"
    )
    parser.add_argument("--input", required=True, help="Merged DAG protobuf")
    parser.add_argument("--output", required=True, help="Trimmed DAG protobuf")
    return parser.parse_args()


def main():
    args = parse_args()
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(
        args.input, compact_genomes=True
    )
    dag.trim_optimal_weight()
    dag.to_protobuf_file(args.output)


if __name__ == "__main__":
    main()
