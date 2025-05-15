# Snakemake pipeline for flu-usher

configfile: "config/config.yaml"

# Define the final outputs that should be created for each type-segment combination
rule all:
    input:
        expand("results/{type}-{segment}/opt_tree.pb.gz", 
               type=config["types"], 
               segment=config["segments"]),
        expand("results/{type}-{segment}/opt_tree.jsonl.gz", 
               type=config["types"], 
               segment=config["segments"])

# Download reference sequence and prepare Nextclade dataset
rule download_reference:
    output:
        directory("results/{type}-{segment}/reference/"),
        gff="results/{type}-{segment}/reference/reference.gff",
        fasta="results/{type}-{segment}/reference/reference.fasta",
        json="results/{type}-{segment}/reference/pathogen.json"
    params:
        accession=lambda wildcards: config["references"][f"{wildcards.type}_{wildcards.segment}"]
    log:
        "logs/{type}-{segment}/download_reference.log"
    shell:
        """
        python scripts/download_ref_seq.py \
            --type {wildcards.type} \
            --segment {wildcards.segment} \
            --accession {params.accession} \
            --output-dir {output[0]} \
            &> {log}
        """

# Align sequences to the reference using nextclade
rule align_sequences:
    input:
        sequences="data/{type}-{segment}/{type}-{segment}.fasta",
        dataset_dir="results/{type}-{segment}/reference/"
    output:
        alignment="results/{type}-{segment}/msa.fasta.xz"
    threads: config["threads"]
    log:
        "logs/{type}-{segment}/nextclade.log"
    shell:
        """
        nextclade run {input.sequences} \
            --input-dataset {input.dataset_dir} \
            --include-reference \
            --output-fasta {output.alignment} \
            >& {log}
        """

# Curate the alignment to only include CDSs and to sanitize IDs
rule curate_alignment:
    input:
        alignment="results/{type}-{segment}/msa.fasta.xz",
        gff="results/{type}-{segment}/reference/reference.gff"
    output:
        curated="results/{type}-{segment}/curated_msa.fasta.xz"
    params:
        gene_name=lambda wildcards: wildcards.segment
    log:
        "logs/{type}-{segment}/curate.log"
    shell:
        """
        python scripts/curate_msa.py \
            --input {input.alignment} \
            --output {output.curated} \
            --gff {input.gff} \
            --gene-name {params.gene_name} \
            --max_frac_gaps {config[max_frac_gaps]} \
            --max_frac_ambig {config[max_frac_ambig]} \
            &> {log}
        """

# Convert the curated MSA to VCF format for UShER
rule create_vcf:
    input:
        "results/{type}-{segment}/curated_msa.fasta.xz"
    output:
        "results/{type}-{segment}/curated_msa.vcf.gz"
    shell:
        """
        xzcat {input} \
        | faToVcf -includeRef -includeNoAltN stdin stdout \
        | gzip > {output}
        """

# Create initial tree with usher-sampled
rule create_initial_tree:
    input:
        vcf="results/{type}-{segment}/curated_msa.vcf.gz"
    output:
        tree="results/{type}-{segment}/preopt_tree.pb.gz",
        empty_tree=temp("results/{type}-{segment}/emptyTree.nwk")
    params:
        outdir="results/{type}-{segment}"
    threads: config["threads"]
    log:
        stdout="logs/{type}-{segment}/usher-sampled.log",
        stderr="logs/{type}-{segment}/usher-sampled.stderr"
    shell:
        """
        echo '()' > {output.empty_tree}
        usher-sampled -T {threads} -A -e 5 \
            -t {output.empty_tree} \
            -v {input.vcf} \
            -d {params.outdir}/ \
            -o {output.tree} \
            --optimization_radius 0 --batch_size_per_process 100 \
            > {log.stdout} 2> {log.stderr}
        """

# Optimize the tree with matOptimize
rule optimize_tree:
    input:
        tree="results/{type}-{segment}/preopt_tree.pb.gz",
        vcf="results/{type}-{segment}/curated_msa.vcf.gz"
    output:
        "results/{type}-{segment}/opt_tree.pb.gz"
    threads: config["threads"]
    log:
        "logs/{type}-{segment}/optimize.log"
    shell:
        """
        matOptimize -T {threads} -m 0.00000001 -M 1 \
            -i {input.tree} \
            -v {input.vcf} \
            -o {output} \
            &> {log}
        """

# Convert the optimized tree to Taxonium format for visualization
rule convert_to_taxonium:
    input:
        "results/{type}-{segment}/opt_tree.pb.gz"
    output:
        "results/{type}-{segment}/opt_tree.jsonl.gz"
    log:
        "logs/{type}-{segment}/taxonium.log"
    shell:
        """
        usher_to_taxonium --input {input} --output {output} &> {log}
        """