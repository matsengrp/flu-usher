import historydag as hdag

input_dag_path = snakemake.input.dag_protobuf
output_dag_path = snakemake.output.trimmed_dag_protobuf

dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_dag_path, compact_genomes=True)
dag.trim_optimal_weight()
dag.to_protobuf_file(output_dag_path)
