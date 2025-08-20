import historydag as hdag

input_pb_name = snakemake.input.dag_protobuf
output_tree_name = snakemake.output.newick

dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_pb_name)
name_func = lambda n: n.label.node_id if (not n.is_ua_node()) else ""

# MB alt name_func to remove internal node names that might mirror leaves:
def alt_name_func(n):
    if n.is_leaf():
        return n.label.node_id
    elif not n.is_ua_node():
        return "internal_"+str(n.label.node_id)
    return ""

with open(output_tree_name, "w") as f:
    f.write(dag.fast_sample().to_newick(name_func=alt_name_func, features=[], feature_funcs={}))
