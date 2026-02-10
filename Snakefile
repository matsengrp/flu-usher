# Snakemake pipeline for flu-usher
import glob
configfile: "config.yaml"


# Define the final outputs that should be created for each segment-subtype combination
rule all:
    input:
        # Taxonium visualization trees for HA segments by subtype
        expand("results/HA/{subtype}/final_tree.jsonl.gz",
               subtype=config["ha_subtypes"]),
        # Unaligned coding sequences for HA segments by subtype (per-gene output directories)
        expand("results/HA/{subtype}/unaligned_coding_seqs/",
               subtype=config["ha_subtypes"]),
        # Host-specific Taxonium visualizations for HA segments by subtype
        expand("results/HA/{subtype}/host_specific_trees/{host_group}_tree.jsonl.gz",
               subtype=config["ha_subtypes"],
               host_group=config["host_groups_to_extract"]),
        # Taxonium visualization trees for NA segments by subtype
        expand("results/NA/{subtype}/final_tree.jsonl.gz",
               subtype=config["na_subtypes"]),
        # Unaligned coding sequences for NA segments by subtype (per-gene output directories)
        expand("results/NA/{subtype}/unaligned_coding_seqs/",
               subtype=config["na_subtypes"]),
        # Host-specific Taxonium visualizations for NA segments by subtype
        expand("results/NA/{subtype}/host_specific_trees/{host_group}_tree.jsonl.gz",
               subtype=config["na_subtypes"],
               host_group=config["host_groups_to_extract"]),
        # Taxonium visualization trees for other segments (all subtypes combined)
        expand("results/{segment}/all/final_tree.jsonl.gz",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        # Unaligned coding sequences for other segments (all subtypes combined, per-gene output directories)
        expand("results/{segment}/all/unaligned_coding_seqs/",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        # Host-specific Taxonium visualizations for other segments (all subtypes combined)
        expand("results/{segment}/all/host_specific_trees/{host_group}_tree.jsonl.gz",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]],
               host_group=config["host_groups_to_extract"]),
        # Executed analysis notebooks
        "results/notebooks/analyze_metadata.html",
        "results/notebooks/analyze_alignments.html"

# Parse GISAID data files from all input directories at once
rule parse_gisaid_data:
    output:
        metadata="results/combined_metadata.csv",
        # Generate output files for each segment-subtype combination
        fastas=expand("results/HA/{subtype}/raw_sequences.fasta.xz",
                     subtype=config["ha_subtypes"]) + \
               expand("results/NA/{subtype}/raw_sequences.fasta.xz", 
                     subtype=config["na_subtypes"]) + \
               expand("results/{segment}/all/raw_sequences.fasta.xz",
                     segment=[s for s in config["segments"] if s not in ["HA", "NA"]])
    params:
        input_dirs=config["input_dirs"],
        segments=" ".join(config["segments"])
    log:
        "logs/parse_gisaid_data.log"
    shell:
        """
        python scripts/parse_gisaid_data.py \
            --input-dirs {params.input_dirs} \
            --output-dir results \
            --segments {params.segments} \
            &> {log}
        """

# Download all reference sequences sequentially with fixed wait time between downloads
# This prevents NCBI rate limiting and provides predictable runtime
rule download_all_references:
    output:
        # Explicitly list all expected outputs for Snakemake dependency tracking
        fastas=expand("results/HA/{subtype}/reference/reference.fasta",
                     subtype=config["ha_subtypes"]) + \
              expand("results/NA/{subtype}/reference/reference.fasta",
                     subtype=config["na_subtypes"]) + \
              expand("results/{segment}/all/reference/reference.fasta",
                     segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        gffs=expand("results/HA/{subtype}/reference/reference.gff",
                   subtype=config["ha_subtypes"]) + \
             expand("results/NA/{subtype}/reference/reference.gff",
                    subtype=config["na_subtypes"]) + \
             expand("results/{segment}/all/reference/reference.gff",
                    segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        jsons=expand("results/HA/{subtype}/reference/pathogen.json",
                    subtype=config["ha_subtypes"]) + \
              expand("results/NA/{subtype}/reference/pathogen.json",
                     subtype=config["na_subtypes"]) + \
              expand("results/{segment}/all/reference/pathogen.json",
                     segment=[s for s in config["segments"] if s not in ["HA", "NA"]])
    params:
        config_file="config.yaml",
        wait_time=30
    log:
        "logs/download_all_references.log"
    shell:
        """
        python scripts/download_ref_seq.py \
            --config {params.config_file} \
            --output-base-dir results \
            --wait-time {params.wait_time} \
            &> {log}
        """

# Use Nextclade to align sequences to the reference
rule align_sequences:
    input:
        sequences="results/{segment}/{subtype}/raw_sequences.fasta.xz",
        reference_gff="results/{segment}/{subtype}/reference/reference.gff",
        reference_fasta="results/{segment}/{subtype}/reference/reference.fasta",
        reference_json="results/{segment}/{subtype}/reference/pathogen.json"
    output:
        alignment="results/{segment}/{subtype}/msa.fasta.xz",
        tsv="results/{segment}/{subtype}/msa.tsv.xz"
    threads: config["threads"]
    log:
        "logs/{segment}/{subtype}/nextclade.log"
    params:
        dataset_dir="results/{segment}/{subtype}/reference/"
    shell:
        """
        nextclade run {input.sequences} \
            --input-dataset {params.dataset_dir} \
            --include-reference \
            --penalty-gap-open-out-of-frame 16 \
            --jobs {threads} \
            --output-fasta {output.alignment} \
            --output-tsv {output.tsv} \
            >& {log}
        """

# Curate alignment and extract unaligned coding sequences in a single pass.
# Filters sequences by quality metrics (gaps, ambiguous nucleotides, terminal gaps, duplicates)
# and validates CDS for all genes. Sequences must pass ALL validations to be included.
# Outputs curated MSA, reference files, and per-gene unaligned coding sequences.
rule curate_and_extract_coding_seqs:
    input:
        alignment="results/{segment}/{subtype}/msa.fasta.xz",
        gff="results/{segment}/{subtype}/reference/reference.gff",
        tsv="results/{segment}/{subtype}/msa.tsv.xz",
        raw_sequences="results/{segment}/{subtype}/raw_sequences.fasta.xz"
    output:
        # Curated MSA and reference files
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz",
        curated_ref_fasta="results/{segment}/{subtype}/curated_reference.fasta",
        curated_ref_txt="results/{segment}/{subtype}/curated_reference.txt",
        curated_ref_gff="results/{segment}/{subtype}/curated_reference.gff",
        curated_ref_gtf="results/{segment}/{subtype}/curated_reference.gtf",
        # Per-gene unaligned coding sequences
        coding_seqs_dir=directory("results/{segment}/{subtype}/unaligned_coding_seqs/")
    params:
        output_dir="results/{segment}/{subtype}",
        output_coding_dir="results/{segment}/{subtype}/unaligned_coding_seqs",
        max_frac_gaps=config["max_frac_gaps"],
        max_frac_ambig=config["max_frac_ambig"]
    log:
        "logs/{segment}/{subtype}/curate_and_extract_coding_seqs.log"
    shell:
        """
        python scripts/curate_and_extract_coding_seqs.py \
            --input {input.alignment} \
            --gff {input.gff} \
            --tsv {input.tsv} \
            --raw-sequences {input.raw_sequences} \
            --output-dir {params.output_dir} \
            --output-coding-dir {params.output_coding_dir} \
            --max_frac_gaps {params.max_frac_gaps} \
            --max_frac_ambig {params.max_frac_ambig} \
            --filter_duplicates \
            &> {log}
        """

# Randomize the alignment order while keeping the reference at the top
rule randomize_alignment:
    input:
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz"
    output:
        "results/{segment}/{subtype}/randomized_{n}/msa.fasta.xz"
    params:
        seed=lambda wildcards: int(wildcards.n)
    log:
        "logs/{segment}/{subtype}/randomize_{n}.log"
    shell:
        """
        python scripts/randomize_alignment.py \
            --input {input.curated_msa} \
            --output {output} \
            --seed {params.seed} \
            &> {log}
        """

# Convert the MSA to VCF format for UShER
rule create_vcf:
    input:
        msa="results/{segment}/{subtype}/randomized_{n}/msa.fasta.xz"
    output:
        "results/{segment}/{subtype}/randomized_{n}/msa.vcf"
    log:
        "logs/{segment}/{subtype}/randomized_{n}/create_vcf.log"
    shell:
        """
        xzcat {input.msa} \
        | faToVcf -includeRef -ambiguousToN -includeNoAltN stdin {output} 2> {log}
        """

# Create initial tree with usher-sampled
# TODO add this arg back in? -e 5
rule create_initial_tree:
    input:
        vcf="results/{segment}/{subtype}/randomized_{n}/msa.vcf"
    output:
        tree="results/{segment}/{subtype}/randomized_{n}/preopt_tree.pb.gz",
        empty_tree=temp("results/{segment}/{subtype}/randomized_{n}/emptyTree.nwk")
    params:
        outdir="results/{segment}/{subtype}/randomized_{n}"
    threads: 2
    log:
        stdout="logs/{segment}/{subtype}/randomized_{n}/usher-sampled.log",
        stderr="logs/{segment}/{subtype}/randomized_{n}/usher-sampled.stderr"
    shell:
        """
        echo '()' > {output.empty_tree}
        usher-sampled -T {threads} -A \
            -t {output.empty_tree} \
            -v {input.vcf} \
            -d {params.outdir}/ \
            -o {output.tree} \
            --optimization_radius 0 --batch_size_per_process 5 \
            > {log.stdout} 2> {log.stderr}
        """

# Optimize the tree with matOptimize
rule optimize_tree:
    input:
        tree="results/{segment}/{subtype}/randomized_{n}/preopt_tree.pb.gz",
        vcf="results/{segment}/{subtype}/randomized_{n}/msa.vcf"
    output:
        "results/{segment}/{subtype}/randomized_{n}/opt_tree.pb.gz"
    threads: config["threads"]
    log:
        "logs/{segment}/{subtype}/randomized_{n}/optimize.log"
    shell:
        """
        matOptimize -T {threads} -m 0.00000001 -M 5 \
            -i {input.tree} \
            -v {input.vcf} \
            -o {output} \
            &> {log}
        """

# Convert optimized trees to DAGs
rule tree_to_dag:
    input:
        tree="results/{segment}/{subtype}/randomized_{n}/opt_tree.pb.gz",
        ref="results/{segment}/{subtype}/curated_reference.txt",
        vcf="results/{segment}/{subtype}/randomized_{n}/msa.vcf"
    output:
        "results/{segment}/{subtype}/randomized_{n}/dag.pb"
    conda:
        "larch"
    log:
        "logs/{segment}/{subtype}/randomized_{n}/tree_to_dag.log"
    shell:
        """
        larch-usher -i {input.tree} \
            -r {input.ref} \
            -v {input.vcf} \
            -o {output} \
            -c 0 \
            &> {log}
        """

# Use larch to merge multiple DAGs into a single DAG
rule larch_merge:
    input:
        dags=expand("results/{{segment}}/{{subtype}}/randomized_{n}/dag.pb",
                    n=range(config["n_randomizations"])),
        vcf="results/{segment}/{subtype}/randomized_0/msa.vcf"
    output:
        "results/{segment}/{subtype}/larch_merged_dag.pb"
    conda:
        "larch"
    log:
        "logs/{segment}/{subtype}/larch_merge.log"
    params:
        dag_args=lambda wildcards, input: " ".join([f"-i {dag}" for dag in input.dags])
    shell:
        """
        larch-dagutil \
            {params.dag_args} \
            -v {input.vcf} \
            --trim \
            -o {output} \
            &> {log}
        """

# Use larch-usher to optimize the trees
# rule larch_optimize:
#     input:
#         dag="results/{segment}/{subtype}/larch_merged_dag.pb",
#     output:
#         opt_dag="results/{segment}/{subtype}/larch_optimized_dag.pb"
#     conda:
#         "larch"
#     threads: config["threads"]
#     log:
#         "logs/{segment}/{subtype}/larch_optimize.log"
#     shell:
#         """
#         larch-usher -i {input.dag} \
#             -o {output.opt_dag} \
#             -c 1 \
#             -s 0 \
#             --max-subtree-clade-size 2000 \
#             --trim \
#             --quiet \
#             &> {log}
#         """

# Trim the DAG to remove suboptimal trees
rule trim_dag:
    input:
        dag_protobuf="results/{segment}/{subtype}/larch_merged_dag.pb"
    output:
        trimmed_dag_protobuf="results/{segment}/{subtype}/trimmed_dag.pb"
    conda:
        "historydag"
    log:
        "logs/{segment}/{subtype}/trim_dag.log"
    script:
        "scripts/trim_dag.py"

# Create a newick tree from the trimmed DAG
rule create_newick:
    input:
        dag_protobuf="results/{segment}/{subtype}/trimmed_dag.pb"
    output:
        newick="results/{segment}/{subtype}/sampled_tree.nh"
    conda:
        "historydag"
    log:
        "logs/{segment}/{subtype}/create_newick.log"
    script:
        "scripts/convert_DAG_protobuf_to_newick_samples.py"

# Create MAT protobuf from newick tree
rule create_mat_protobuf:
    input:
        nh_file="results/{segment}/{subtype}/sampled_tree.nh",
        vcf_file="results/{segment}/{subtype}/randomized_0/msa.vcf"
    output:
        protobuf_name="results/{segment}/{subtype}/sampled_tree.pb.gz"
    log:
        "logs/{segment}/{subtype}/create_mat_protobuf.log"
    shell:
        """
        matOptimize -t {input.nh_file} -v {input.vcf_file} -o {output.protobuf_name} -N 0 &> {log}
        """

# Reroot the tree if specified in config, otherwise create symlink
rule reroot_tree:
    input:
        "results/{segment}/{subtype}/sampled_tree.pb.gz"
    output:
        "results/{segment}/{subtype}/final_tree.pb.gz"
    log:
        "logs/{segment}/{subtype}/reroot.log"
    run:
        key = f"{wildcards.segment}_{wildcards.subtype}"
        if key in config.get("reroot", {}):
            new_root = config["reroot"][key]
            shell("""
                matUtils extract -i {input} \
                    -o {output} \
                    -y {new_root} \
                    &> {log}
            """)
        else:
            shell("""
                ln -sf $(basename {input}) {output} \
                && echo "Created symlink for {key} (no rerooting specified)" > {log}
            """)

# Extract root sequence from tree mutations or create symlink to reference
rule create_root_fasta:
    input:
        tree="results/{segment}/{subtype}/sampled_tree.pb.gz",
        msa="results/{segment}/{subtype}/curated_msa.fasta.xz",
        ref="results/{segment}/{subtype}/curated_reference.fasta"
    output:
        "results/{segment}/{subtype}/curated_root.fasta"
    params:
        new_root=lambda wildcards: config.get("reroot", {}).get(f"{wildcards.segment}_{wildcards.subtype}", None)
    log:
        "logs/{segment}/{subtype}/create_root_fasta.log"
    run:
        if params.new_root:
            # Rerooting specified - infer root sequence from tree
            samples_file = f"results/{wildcards.segment}/{wildcards.subtype}/root_samples.txt"
            paths_file = f"results/{wildcards.segment}/{wildcards.subtype}/root_paths.txt"

            # Create samples file (reference + new root)
            shell("""
                python scripts/create_root_samples_file.py \
                    --reference {input.ref} \
                    --new-root {params.new_root} \
                    --output {samples_file} \
                    &>> {log}
            """)

            # Extract mutation paths using matUtils
            shell("""
                matUtils extract -i {input.tree} \
                    -s {samples_file} \
                    -S {paths_file} \
                    &>> {log}
            """)

            # Infer root sequence from mutations and validate against MSA
            shell("""
                python scripts/extract_root_sequence.py \
                    --reference {input.ref} \
                    --msa {input.msa} \
                    --paths {paths_file} \
                    --new-root-name {params.new_root} \
                    --output {output} \
                    &>> {log}
            """)
        else:
            # No rerooting specified, use reference
            shell("""
                ln -sf $(basename {input.ref}) {output} \
                && echo "Created symlink (no rerooting specified)" > {log}
            """)


# Add host group classifications to combined metadata
rule add_host_groups:
    input:
        metadata="results/combined_metadata.csv"
    output:
        "results/combined_metadata_with_host_groups.csv"
    log:
        "logs/add_host_groups.log"
    shell:
        """
        python scripts/simplified_host_classifier.py {input.metadata} {output} 2> {log}
        """

# Create samples file for host-specific subtree extraction
rule create_host_samples_file:
    input:
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz",
        metadata="results/combined_metadata_with_host_groups.csv",
        root="results/{segment}/{subtype}/curated_root.fasta"
    output:
        "results/{segment}/{subtype}/host_specific_trees/{host_group}_samples.txt"
    log:
        "logs/{segment}/{subtype}/create_host_samples_{host_group}.log"
    shell:
        """
        python scripts/create_host_samples_file.py \
            --curated-msa {input.curated_msa} \
            --metadata {input.metadata} \
            --host-group {wildcards.host_group} \
            --root {input.root} \
            --output {output} \
            &> {log}
        """

# Extract host-specific subtree using matUtils, collapsing trees before outputting
rule extract_host_subtree:
    input:
        tree="results/{segment}/{subtype}/final_tree.pb.gz",
        samples="results/{segment}/{subtype}/host_specific_trees/{host_group}_samples.txt"
    output:
        "results/{segment}/{subtype}/host_specific_trees/{host_group}_tree.pb.gz"
    log:
        "logs/{segment}/{subtype}/extract_host_subtree_{host_group}.log"
    shell:
        """
        matUtils extract \
            -i {input.tree} \
            -s {input.samples} \
            -O \
            -o {output} \
            &> {log}
        """

# Convert the final tree to Taxonium format for visualization
rule convert_to_taxonium:
    input:
        final_tree="results/{segment}/{subtype}/final_tree.pb.gz",
        metadata="results/combined_metadata_with_host_groups.csv"
    output:
        "results/{segment}/{subtype}/final_tree.jsonl.gz"
    log:
        "logs/{segment}/{subtype}/taxonium.log"
    shell:
        """
        usher_to_taxonium \
            --input {input.final_tree} \
            --metadata {input.metadata} \
            --key_column isolate_id \
            --columns isolate_name,subtype,clade,passage_history,location,host,host_group,collection_date \
            --output {output} \
            &> {log}
        """

# Convert host-specific subtrees to Taxonium format for visualization
rule convert_host_subtree_to_taxonium:
    input:
        tree="results/{segment}/{subtype}/host_specific_trees/{host_group}_tree.pb.gz",
        metadata="results/combined_metadata_with_host_groups.csv"
    output:
        "results/{segment}/{subtype}/host_specific_trees/{host_group}_tree.jsonl.gz"
    log:
        "logs/{segment}/{subtype}/taxonium_host_{host_group}.log"
    shell:
        """
        usher_to_taxonium \
            --input {input.tree} \
            --metadata {input.metadata} \
            --key_column isolate_id \
            --columns isolate_name,subtype,clade,passage_history,location,host,host_group,collection_date \
            --output {output} \
            &> {log}
        """

# Execute analysis notebooks and generate HTML reports
rule execute_notebooks:
    input:
        notebook="notebooks/{notebook}.ipynb",
        # Ensure all main pipeline outputs are complete before running notebooks
        metadata="results/combined_metadata_with_host_groups.csv",
        ha_trees=expand("results/HA/{subtype}/final_tree.jsonl.gz",
                       subtype=config["ha_subtypes"]),
        na_trees=expand("results/NA/{subtype}/final_tree.jsonl.gz",
                       subtype=config["na_subtypes"]),
        other_trees=expand("results/{segment}/all/final_tree.jsonl.gz",
                          segment=[s for s in config["segments"] if s not in ["HA", "NA"]])
    output:
        "results/notebooks/{notebook}.html"
    log:
        "logs/notebooks/{notebook}.log"
    shell:
        """
        mkdir -p results/notebooks
        jupyter nbconvert --to html \
            --execute \
            --ExecutePreprocessor.timeout=600 \
            --output-dir results/notebooks \
            --output {wildcards.notebook}.html \
            {input.notebook} \
            &> {log}
        """
