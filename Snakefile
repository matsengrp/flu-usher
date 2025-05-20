# Snakemake pipeline for flu-usher
import glob
configfile: "config/config.yaml"

# Define the final outputs that should be created for each subtype-segment combination
rule all:
    input:
        expand("results/{subtype}/{segment}/opt_tree.pb.gz", 
               subtype=config["subtypes"], 
               segment=config["segments"]),
        expand("results/{subtype}/{segment}/opt_tree.jsonl.gz", 
               subtype=config["subtypes"], 
               segment=config["segments"])

# Parse GISAID data files by segment
rule parse_gisaid_data:
    output:
        metadata="results/{subtype}/{subtype}_metadata.csv",
        fastas=expand("results/{{subtype}}/{segment}/raw_sequences.fasta.xz", 
                      segment=config["segments"])
    log:
        "logs/{subtype}/parse_gisaid_data.log"
    shell:
        """
        python scripts/parse_gisaid_data.py \
            --subtype {wildcards.subtype} \
            --output-dir results/{wildcards.subtype} \
            --segments {config[segments]} \
            &> {log}
        """

# Download the reference sequence and corresponding GFF, create a Nextclade pathogen.json file,
# and save these files in the same directory as input to Nextclade
rule download_reference:
    output:
        directory("results/{subtype}/{segment}/reference/"),
        gff="results/{subtype}/{segment}/reference/reference.gff",
        fasta="results/{subtype}/{segment}/reference/reference.fasta",
        json="results/{subtype}/{segment}/reference/pathogen.json"
    params:
        accession=lambda wildcards: config["references"][f"{wildcards.subtype}_{wildcards.segment}"]
    log:
        "logs/{subtype}/{segment}/download_reference.log"
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
        sequences="results/{subtype}/{segment}/raw_sequences.fasta.xz",
        dataset_dir="results/{subtype}/{segment}/reference/",
        reference_gff="results/{subtype}/{segment}/reference/reference.gff",
        reference_fasta="results/{subtype}/{segment}/reference/reference.fasta",
        reference_json="results/{subtype}/{segment}/reference/pathogen.json"
    output:
        alignment="results/{subtype}/{segment}/msa.fasta.xz"
    threads: config["threads"]
    log:
        "logs/{subtype}/{segment}/nextclade.log"
    shell:
        """
        nextclade run {input.sequences} \
            --input-dataset {input.dataset_dir} \
            --include-reference \
            --output-fasta {output.alignment} \
            >& {log}
        """

# Curate the alignment to only include sites in the CDS of the reference, and sanitize IDs
# so that they do not include special characters that lead to errors in UShER jobs
rule curate_alignment:
    input:
        alignment="results/{subtype}/{segment}/msa.fasta.xz",
        gff="results/{subtype}/{segment}/reference/reference.gff"
    output:
        curated="results/{subtype}/{segment}/curated_msa.fasta.xz",
        curated_gff="results/{subtype}/{segment}/curated_reference.gff"
    params:
        gene_name=lambda wildcards: wildcards.segment
    log:
        "logs/{subtype}/{segment}/curate.log"
    shell:
        """
        python scripts/curate_msa.py \
            --input {input.alignment} \
            --output {output.curated} \
            --gff {input.gff} \
            --gene-name {params.gene_name} \
            --output-gff {output.curated_gff} \
            --max_frac_gaps {config[max_frac_gaps]} \
            --max_frac_ambig {config[max_frac_ambig]} \
            &> {log}
        """

# Convert the curated MSA to VCF format for UShER
rule create_vcf:
    input:
        curated_msa="results/{subtype}/{segment}/curated_msa.fasta.xz"
    output:
        "results/{subtype}/{segment}/curated_msa.vcf.gz"
    shell:
        """
        xzcat {input.curated_msa} \
        | faToVcf -includeRef -includeNoAltN stdin stdout \
        | gzip > {output}
        """

# Create initial tree with usher-sampled
rule create_initial_tree:
    input:
        vcf="results/{subtype}/{segment}/curated_msa.vcf.gz"
    output:
        tree="results/{subtype}/{segment}/preopt_tree.pb.gz",
        empty_tree=temp("results/{subtype}/{segment}/emptyTree.nwk")
    params:
        outdir="results/{subtype}/{segment}"
    threads: config["threads"]
    log:
        stdout="logs/{subtype}/{segment}/usher-sampled.log",
        stderr="logs/{subtype}/{segment}/usher-sampled.stderr"
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
        tree="results/{subtype}/{segment}/preopt_tree.pb.gz",
        vcf="results/{subtype}/{segment}/curated_msa.vcf.gz"
    output:
        "results/{subtype}/{segment}/opt_tree.pb.gz"
    threads: config["threads"]
    log:
        "logs/{subtype}/{segment}/optimize.log"
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
        opt_tree="results/{subtype}/{segment}/opt_tree.pb.gz"
    output:
        "results/{subtype}/{segment}/opt_tree.jsonl.gz"
    log:
        "logs/{subtype}/{segment}/taxonium.log"
    shell:
        """
        usher_to_taxonium --input {input.opt_tree} --output {output} &> {log}
        """