# Snakemake pipeline for flu-usher
import glob
configfile: "config/config.yaml"

# Define the final outputs that should be created for each segment-subtype combination
rule all:
    input:
        # Trees for HA segments by subtype
        expand("results/HA/{subtype}/opt_tree.pb.gz", 
               subtype=config["ha_subtypes"]),
        expand("results/HA/{subtype}/opt_tree.jsonl.gz", 
               subtype=config["ha_subtypes"]),
        # Trees for NA segments by subtype
        expand("results/NA/{subtype}/opt_tree.pb.gz", 
               subtype=config["na_subtypes"]),
        expand("results/NA/{subtype}/opt_tree.jsonl.gz", 
               subtype=config["na_subtypes"]),
        # Trees for other segments (all subtypes combined)
        expand("results/{segment}/all/opt_tree.pb.gz",
               segment=[s for s in config["segments"] if s not in ["HA", "NA"]]),
        expand("results/{segment}/all/opt_tree.jsonl.gz",
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
        input_dirs=config["input_dirs"]
    log:
        "logs/parse_gisaid_data.log"
    shell:
        """
        python scripts/parse_gisaid_data.py \
            --input-dirs {params.input_dirs} \
            --output-dir results \
            --segments {config[segments]} \
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
        alignment="results/{segment}/{subtype}/msa.fasta.xz"
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
            >& {log}
        """

# Curate the alignment to only include sites in the CDS of the reference, and sanitize IDs
# so that they do not include special characters that lead to errors in UShER jobs. Also
# output new FASTA, GFF, and GTF files for the curated reference sequence.
rule curate_alignment:
    input:
        alignment="results/{segment}/{subtype}/msa.fasta.xz",
        gff="results/{segment}/{subtype}/reference/reference.gff"
    output:
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz",
        curated_ref_fasta="results/{segment}/{subtype}/curated_reference.fasta",
        curated_ref_gff="results/{segment}/{subtype}/curated_reference.gff",
        curated_ref_gtf="results/{segment}/{subtype}/curated_reference.gtf"
    params:
        output_dir="results/{segment}/{subtype}"
    log:
        "logs/{segment}/{subtype}/curate.log"
    shell:
        """
        python scripts/curate_msa.py \
            --input {input.alignment} \
            --output-dir {params.output_dir} \
            --gff {input.gff} \
            --max_frac_gaps {config[max_frac_gaps]} \
            --max_frac_ambig {config[max_frac_ambig]} \
            &> {log}
        """

# Convert the curated MSA to VCF format for UShER
rule create_vcf:
    input:
        curated_msa="results/{segment}/{subtype}/curated_msa.fasta.xz"
    output:
        "results/{segment}/{subtype}/curated_msa.vcf.gz"
    shell:
        """
        xzcat {input.curated_msa} \
        | faToVcf -includeRef -includeNoAltN stdin stdout \
        | gzip > {output}
        """

# Create initial tree with usher-sampled
rule create_initial_tree:
    input:
        vcf="results/{segment}/{subtype}/curated_msa.vcf.gz"
    output:
        tree="results/{segment}/{subtype}/preopt_tree.pb.gz",
        empty_tree=temp("results/{segment}/{subtype}/emptyTree.nwk")
    params:
        outdir="results/{segment}/{subtype}"
    threads: config["threads"]
    log:
        stdout="logs/{segment}/{subtype}/usher-sampled.log",
        stderr="logs/{segment}/{subtype}/usher-sampled.stderr"
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
        tree="results/{segment}/{subtype}/preopt_tree.pb.gz",
        vcf="results/{segment}/{subtype}/curated_msa.vcf.gz"
    output:
        "results/{segment}/{subtype}/opt_tree.pb.gz"
    threads: config["threads"]
    log:
        "logs/{segment}/{subtype}/optimize.log"
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
        opt_tree="results/{segment}/{subtype}/opt_tree.pb.gz",
        metadata="results/combined_metadata.csv"
    output:
        "results/{segment}/{subtype}/opt_tree.jsonl.gz"
    log:
        "logs/{segment}/{subtype}/taxonium.log"
    shell:
        """
        usher_to_taxonium \
            --input {input.opt_tree} \
            --metadata {input.metadata} \
            --key_column isolate_id \
            --columns isolate_name,subtype,clade,passage_history,location,host,collection_date \
            --output {output} \
            &> {log}
        """

