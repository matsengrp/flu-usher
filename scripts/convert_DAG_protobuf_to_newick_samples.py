import historydag as hdag

input_dag_path = snakemake.input.dag_protobuf
output_newick_path = snakemake.output.newick

dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_dag_path)


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


with open(output_newick_path, "w") as f:
    f.write(
        dag.fast_sample().to_newick(
            name_func=name_internal_nodes_distinctly, features=[], feature_funcs={}
        )
    )
