# Notes from Angie on how to build pipeline

## TODO

* Incorporating metadata
   * When downloading, only include EPI_ISL, subtype, segment
   * Download all segments at once, then add a step that splits them into separate files?
   * In curation step, simplify header to only include EPI_ISL
   * Add CSV of metadata file with a specific subset of columns
   

## Align sequences to reference

See: `run_nextclade.sh`

Run nextclade like:

nextclade run --input-dataset $nextcladeTypeSegDir \
                --include-reference \
                --jobs $threads \
                --output-fasta aligned.$type.$segment.fa.xz \
                >& nextalign.log

-- where $nextcladeTypeSegDir is a directory containing the reference sequence and nextclade config settings.

`data/NC_026425.1/` is what she uses for H7N9 HA

## Curate the alignment

See: `notebooks/curate_alignment.ipynb`

## Convert MSA to VCF

See: `create_vcf.sh`

Then aligned.$type.$segment.fa.xz is an MSA fasta file with all of the sequences aligned to reference.  UShER doesn't read MSA fasta, it needs the input converted to VCF or the input format for the Maple aligner.  I wrote a little util that converts an MSA to VCF, faToVcf, available from bioconda (conda install ucsc-fatovcf).  Run faToVcf like this:
xzcat aligned.$type.$segment.fa.xz \
| faToVcf -includeRef -includeNoAltN stdin stdout \
| gzip > aligned.$type.$segment.vcf.gz

# Build tree

Then (assuming you have already done 'conda install usher') run usher-sampled like this:
echo '()' > emptyTree.nwk
usher-sampled -T $threads -A -e 5 \
    -t emptyTree.nwk \
    -v aligned.$type.$segment.vcf.gz \
    -o fluA.$type.$segment.preFilter.pb.gz \
    --optimization_radius 0 --batch_size_per_process 100 \
    > usher-sampled.$type.$segment.log 2> usher-sampled.$type.$segment.stderr

This step may not be necessary since you have already limited your input sequences to H7N9, but I apply a branch length filter and max number of private mutations filter:
matUtils extract -i fluA.$type.$segment.preFilter.pb.gz \
    --max-branch-length 250 \
    --max-parsimony 100 \
    -O -o fluA.$type.$segment.preOpt.pb.gz
Then I run matOptimize:
matOptimize -T $threads -m 0.00000001 -M 1 \
    -i fluA.$type.$segment.preOpt.pb.gz \
    -v aligned.$type.$segment.vcf.gz \
    -o fluA.$type.$segment.pb.gz

# Taxomium

Then I run usher_to_taxonium so I can view the tree in taxonium.  Let me know if you want to do that too.  The options to pass to usher_to_taxonium depend on what your metadata file looks like.