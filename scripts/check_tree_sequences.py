"""
Assert that rerooting did not change any sequence stored in the tree.

Rerooting moves the root and redistributes mutations along the backbone, but the
sequence the tree implies for each sample, read against the reference, must not
change. `matUtils extract -y` broke that: it rebased the MAT into the new root's
coordinate frame, so sequences read against curated_reference.fasta came out
shifted by the target's divergence from it (issue #49). The rebased tree was
self-consistent and recoverable if you knew the new root's sequence, but nothing
in the artifact recorded the shift. This check is the standing guard that the
MAT stays in the reference's frame.

Compares the sequences implied by sampled_tree.pb.gz (pre-reroot) and
final_tree.pb.gz (post-reroot), via the root-to-tip mutation paths that
`matUtils extract -S` writes, and exits nonzero if any sample differs.

What this does and does not prove. Both trees are matOptimize products of the
same VCF, and matOptimize assigns mutations that reproduce the input genotypes
exactly, so composing a root-to-tip path returns the alignment's genotype for
any topology. This check therefore passes on a randomly scrambled tree; it is a
coordinate-frame guard, not a rerooting-correctness guard. It catches a MAT that
is not in the reference's frame -- the #49 bug -- and nothing else. That the
tree is rooted where config.yaml asked is asserted separately, here via
--expect-root and in reroot_newick.py; do not read a pass here as confirming it.

Sites the input alignment leaves uncalled are excluded. faToVcf encodes a gap
as missing data rather than as an allele, because a MAT stores substitutions
only; matOptimize then fills those positions in by parsimony from the parent,
so the imputed base legitimately moves when the root moves. Those bases were
never observed, and asserting that a guess is stable under rerooting is not a
property worth having. The count of excluded sites is reported so a combo with
a large gappy region cannot quietly hollow out the check.
"""
import argparse
import sys

from utils import setup_logging

logger = setup_logging(__name__)


def read_reference(path):
    """Return {1-based position: base} for the single record in a FASTA file."""
    seq = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    return {i: base.upper() for i, base in enumerate("".join(seq), start=1)}


def compose(path_field, ref):
    """Collapse a root-to-tip mutation path into {position: allele != reference}.

    `matUtils extract -S` writes space-separated `node:MUT,MUT` chunks in
    root-to-tip order, where each mutation is <parent base><position><new base>.
    Later chunks override earlier ones at the same position, so a mutation and
    its back-mutation cancel and drop out.
    """
    alleles = {}
    for chunk in path_field.split():
        _, _, muts = chunk.partition(":")
        for mut in muts.split(","):
            if not mut:
                continue
            pos = int(mut[1:-1])
            alleles[pos] = mut[-1]
    return {p: a for p, a in alleles.items() if ref.get(p) != a}


def read_paths(path, ref, keep=None):
    """Yield (sample, genotype dict) from a matUtils -S file.

    `keep` restricts to a subset of sample names, so the second pass can hold
    full genotypes for the handful of samples that disagreed without ever
    materialising them for the whole tree.
    """
    with open(path) as handle:
        for line in handle:
            sample, _, path_field = line.rstrip("\n").partition("\t")
            if not sample or (keep is not None and sample not in keep):
                continue
            yield sample, compose(path_field, ref)


def genotype_key(genotype):
    return ";".join(f"{p}{genotype[p]}" for p in sorted(genotype))


def _node_label(token):
    """Name of the node a top-level newick token describes.

    "EPI_ISL_1:55" -> "EPI_ISL_1";  "(a,b)node_2:3" -> "node_2".
    """
    token = token.strip()
    close = token.rfind(")")
    tail = token[close + 1:] if close != -1 else token
    return tail.split(":")[0]


def root_children(newick):
    """Names of the root's direct children in a newick string.

    Hand-rolled rather than via ete3 so this stays in envs/python.yaml: loading
    a 130k-leaf tree just to read two names is not worth a second conda env in
    the rule.
    """
    text = newick.strip()
    end = text.rfind(")")
    if not text.startswith("(") or end == -1:
        raise ValueError("not a newick with an internal root node")
    inner = text[1:end]
    tokens, depth, start = [], 0, 0
    for i, char in enumerate(inner):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            tokens.append(inner[start:i])
            start = i + 1
    tokens.append(inner[start:])
    return [_node_label(t) for t in tokens]


def uncalled_positions(vcf_path, samples, positions):
    """Return {sample: {position, ...}} for genotypes the input left uncalled.

    Restricted to the candidate samples and positions that actually disagreed,
    so this stays cheap on the 700 MB internal-segment VCFs.
    """
    uncalled = {s: set() for s in samples}
    with open(vcf_path) as handle:
        columns = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                columns = [(i, n) for i, n in enumerate(fields) if i >= 9 and n in uncalled]
                continue
            # Reject on position before paying for a full split.
            tab = line.find("\t")
            end = line.find("\t", tab + 1)
            if int(line[tab + 1:end]) not in positions:
                continue
            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            for index, name in columns:
                if fields[index] == ".":
                    uncalled[name].add(pos)
    return uncalled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampled-paths", required=True,
                        help="matUtils -S output for sampled_tree.pb.gz")
    parser.add_argument("--final-paths", required=True,
                        help="matUtils -S output for final_tree.pb.gz")
    parser.add_argument("--reference", required=True,
                        help="curated_reference.fasta the VCF was built against")
    parser.add_argument("--input-vcf", required=True,
                        help="randomized_0/msa.vcf, for the called-site mask")
    parser.add_argument("--output", required=True, help="Report file")
    parser.add_argument("--final-newick", default=None,
                        help="newick of final_tree.pb.gz, for the root check")
    parser.add_argument("--expect-root", default="",
                        help="configured reroot target; empty if not rerooted")
    args = parser.parse_args()

    # Fail before the expensive comparison. This is the invariant the pipeline
    # actually depends on and the one the sequence comparison cannot see.
    if args.expect_root:
        if not args.final_newick:
            logger.error("--expect-root given without --final-newick")
            sys.exit(1)
        with open(args.final_newick) as handle:
            children = root_children(handle.read())
        if args.expect_root not in children:
            logger.error(
                f"final tree is not rooted at '{args.expect_root}'; "
                f"root children are {children}"
            )
            sys.exit(1)
        logger.info(f"Final tree is rooted at '{args.expect_root}'")

    ref = read_reference(args.reference)
    logger.info(f"Reference spans {len(ref)} positions")

    # First pass keeps one hash per sample rather than a genotype map, so the
    # 126k-sample internal segments stay in memory.
    before = {s: hash(genotype_key(g)) for s, g in read_paths(args.sampled_paths, ref)}
    candidates = set()
    seen = set()
    for sample, genotype in read_paths(args.final_paths, ref):
        seen.add(sample)
        if sample not in before or before[sample] != hash(genotype_key(genotype)):
            candidates.add(sample)

    missing = set(before) - seen
    extra = seen - set(before)
    if missing or extra:
        logger.error(
            f"sample sets differ between trees: {len(missing)} only in sampled_tree, "
            f"{len(extra)} only in final_tree"
        )
        sys.exit(1)

    n_samples = len(before)
    # An empty paths file is not a pass. matUtils extract -S writes nothing at
    # all when it dislikes its -S argument (see the comment on the rule that
    # produces these), and Snakemake accepts the empty file as a satisfied
    # output, so without this the whole gate reports success having compared
    # nothing.
    if n_samples == 0:
        logger.error(
            f"no samples in {args.sampled_paths}; refusing to pass a check "
            "that compared nothing"
        )
        sys.exit(1)

    logger.info(f"{n_samples} samples compared, {len(candidates)} differ before masking")

    # Second pass materialises genotypes only for the samples that disagreed.
    differing = {}
    excluded_sites = 0
    if candidates:
        sampled = dict(read_paths(args.sampled_paths, ref, keep=candidates))
        final = dict(read_paths(args.final_paths, ref, keep=candidates))
        positions = set()
        for sample in candidates:
            a, b = sampled[sample], final[sample]
            positions |= set(a) ^ set(b)
            positions |= {p for p in set(a) & set(b) if a[p] != b[p]}
        uncalled = uncalled_positions(args.input_vcf, candidates, positions)
        for sample in sorted(candidates):
            skip = uncalled[sample]
            a = {p: v for p, v in sampled[sample].items() if p not in skip}
            b = {p: v for p, v in final[sample].items() if p not in skip}
            excluded_sites += len(skip)
            if a != b:
                sites = sorted(set(a) ^ set(b)) or sorted(p for p in a if a[p] != b.get(p))
                differing[sample] = sites

    with open(args.output, "w") as handle:
        handle.write(f"samples compared: {n_samples}\n")
        handle.write(f"identical sequences: {n_samples - len(differing)}\n")
        handle.write(f"differing sequences: {len(differing)}\n")
        handle.write(f"samples differing only at uncalled sites: "
                     f"{len(candidates) - len(differing)}\n")
        handle.write(f"uncalled sites excluded: {excluded_sites}\n")
        for sample, sites in list(differing.items())[:20]:
            handle.write(f"  {sample}: {sites[:10]}\n")

    if differing:
        # Snakemake deletes a failed job's outputs, so args.output is gone by
        # the time anyone reads the error. The log survives; put the detail
        # there too rather than pointing at a file that no longer exists.
        logger.error(f"{len(differing)} of {n_samples} sequences changed across rerooting")
        for sample, sites in list(differing.items())[:20]:
            logger.error(f"  {sample}: {sites[:10]}")
        if len(differing) > 20:
            logger.error(f"  ... and {len(differing) - 20} more samples")
        sys.exit(1)

    logger.info(
        f"All {n_samples} sequences identical across rerooting "
        f"({excluded_sites} uncalled sites excluded)"
    )


if __name__ == "__main__":
    main()
