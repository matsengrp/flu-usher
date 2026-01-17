# Snakemake pipeline for flu-usher
import glob
configfile: "config/config.yaml"


# Define the final outputs that should be created for each segment-subtype combination
rule all:
    input:
        # Trees for HA segments by subtype
        expand("results/HA/{subtype}/final_tree.pb.gz",
               subtype=config["ha_subtypes"]),
        expand("results/HA/{subtype}/final_tree.jsonl.gz",
               subtype=config["ha_subtypes"]),
        # Unaligned coding sequences for HA segments by subtype
        expand("results/HA/{subtype}/curated_unaligned_coding_seqs.fasta.xz",
               subtype=config["ha_subtypes"]),
        # Trees for NA segments by subtype
        expand("results/NA/{subtype}/final_tree.pb.gz",
               subtype=config["na_subtypes"]),
        expand("results/NA/{subtype}/final_tree.jsonl.gz",
               subtype=config["na_subtypes"]),
        # Unaligned coding sequences for NA segments by subtype
        expand("results/NA/{subtype}/curated_unaligned_coding_seqs.fasta.xz",
               subtype=config["na_subtypes"]),
        # Trees for other segments (all subtypes combined)
        expand("results/{segment}/all/final_tree.pb.gz",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        expand("results/{segment}/all/final_tree.jsonl.gz",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        # Unaligned coding sequences for other segments (all subtypes combined)
        expand("results/{segment}/all/curated_unaligned_coding_seqs.fasta.xz",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]])

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

# Download the reference sequence and corresponding GFF, create a Nextclade pathogen.json file,
# and save these files in the same directory as input to Nextclade
rule download_reference:
    output:
        directory("results/{segment}/{subtype}/reference/"),
        gff="results/{segment}/{subtype}/reference/reference.gff",
        fasta="results/{segment}/{subtype}/reference/reference.fasta",
        json="results/{segment}/{subtype}/reference/pathogen.json"
    params:
        accession=lambda wildcards: config["references"][f"{wildcards.segment}_{wildcards.subtype}"]
    log:
        "logs/{segment}/{subtype}/download_reference.log"
    shell:
        """
        python scripts/download_ref_seq.py \
            --accession {params.accession} \
            --output-dir {output[0]} \
            &> {log}
        """

# Use Nextclade to align sequences to the reference
rule align_sequences:
    input:
        sequences="results/{segment}/{subtype}/raw_sequences.fasta.xz",
        dataset_dir="results/{segment}/{subtype}/reference/",
        reference_gff="results/{segment}/{subtype}/reference/reference.gff",
        reference_fasta="results/{segment}/{subtype}/reference/reference.fasta",
        reference_json="results/{segment}/{subtype}/reference/pathogen.json"
    output:
        alignment="results/{segment}/{subtype}/msa.fasta.xz",
        tsv="results/{segment}/{subtype}/msa.tsv.xz"
    threads: config["threads"]
    log:
        "logs/{segment}/{subtype}/nextclade.log"
    shell:
        """
        nextclade run {input.sequences} \
            --input-dataset {input.dataset_dir} \
            --include-reference \
            --jobs {threads} \
            --output-fasta {output.alignment} \
            --output-tsv {output.tsv} \
            >& {log}
        """

# Curate the alignment to only include sites in the CDS of the reference, and sanitize IDs
# so that they do not include special characters that lead to errors in UShER jobs. Also
# output new FASTA, GFF, and GTF files for the curated reference sequence.
# --replace_gaps_with_ref
rule curate_alignment:
    input:
        alignment="results/{segment}/{subtype}/msa.fasta.xz",
        gff="results/{segment}/{subtype}/reference/reference.gff"
    output:
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz",
        curated_ref_fasta="results/{segment}/{subtype}/curated_reference.fasta",
        curated_ref_txt="results/{segment}/{subtype}/curated_reference.txt",
        curated_ref_gff="results/{segment}/{subtype}/curated_reference.gff",
        curated_ref_gtf="results/{segment}/{subtype}/curated_reference.gtf"
    params:
        output_dir="results/{segment}/{subtype}",
        max_frac_gaps=config["max_frac_gaps"],
        max_frac_ambig=config["max_frac_ambig"]
    log:
        "logs/{segment}/{subtype}/curate.log"
    shell:
        """
        python scripts/curate_msa.py \
            --input {input.alignment} \
            --output-dir {params.output_dir} \
            --gff {input.gff} \
            --max_frac_gaps {params.max_frac_gaps} \
            --max_frac_ambig {params.max_frac_ambig} \
            --filter_duplicates \
            &> {log}
        """

# Create unaligned coding sequences from curated aligned sequences
# by removing gaps and re-inserting nucleotides that were removed during alignment
rule create_unaligned_coding_sequences:
    input:
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz",
        tsv="results/{segment}/{subtype}/msa.tsv.xz",
        gff="results/{segment}/{subtype}/reference/reference.gff",
        raw_sequences="results/{segment}/{subtype}/raw_sequences.fasta.xz"
    output:
        "results/{segment}/{subtype}/curated_unaligned_coding_seqs.fasta.xz"
    log:
        "logs/{segment}/{subtype}/create_unaligned_coding_seqs.log"
    shell:
        """
        python scripts/create_unaligned_coding_seqs.py \
            --curated-msa {input.curated_msa} \
            --tsv {input.tsv} \
            --gff {input.gff} \
            --raw-sequences {input.raw_sequences} \
            --output {output} \
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


# TODO condense trees before feeding them to taxonium

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
        python scripts/add_host_groups.py {input.metadata} {output} 2> {log}
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
