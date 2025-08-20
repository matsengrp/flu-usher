import historydag as hdag

input_pb_name = snakemake.input.dag_protobuf
output_tree_name = snakemake.output.trimmed_dag_protobuf

dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_pb_name, compact_genomes=True)
dag.trim_optimal_weight()
dag.to_protobuf_file(output_tree_name)

