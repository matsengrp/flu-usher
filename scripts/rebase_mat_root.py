"""
Rebase a MAT so its root node is the origin instead of the reference.

usher records the root's divergence from curated_reference.fasta as a mutation
list on the root node, so reconstructing a sample means starting from the
reference and applying the whole root-to-tip path, the root's own mutations
included. That is correct but leaves the tree's origin implicit: the root has a
sequence, and nothing ships it.

This rewrites the tree so the root carries no mutations and writes its sequence
out alongside. Nothing else about the tree changes, because nothing else
depends on the choice: every branch below the root records a parent-to-child
change (par_nuc -> mut_nuc), which is what it is regardless of what the origin
is. Only two things move -- the root's own mutation list, which is emptied, and
each mutation's ref_nuc annotation, which is repointed at the new origin so it
does not keep describing a sequence the file no longer references.

This is the bookkeeping `matUtils extract -y` used to do, minus its two
problems: it refused whenever the root carried mutations (which is now always,
see issue #49), and it discarded the root's sequence rather than emitting it,
so the shift went unrecorded and anything pairing the tree with the reference
misread it silently.

Imports the MAT schema from taxoniumtools, which vendors the compiled
parsimony.proto. UShER ships the .proto but no Python bindings, and taxonium is
already a pipeline dependency, so this avoids vendoring a generated module that
would drift from the version usher actually writes.
"""
import argparse
import gzip
import sys

from utils import setup_logging

logger = setup_logging(__name__)

# usher encodes nucleotides as indices into this string, verified against the
# mutation paths matUtils writes (position 33 ref_nuc=0 mut_nuc=3 is "A33T").
NUCLEOTIDES = "ACGT"


def read_reference(path):
    """Return the single sequence in a FASTA file, uppercased."""
    parts = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if parts:
                    break
                continue
            parts.append(line.strip())
    return "".join(parts).upper()


def _schema():
    """Imported lazily so the pure helpers stay testable outside envs/taxonium."""
    from taxoniumtools import parsimony_pb2
    return parsimony_pb2


def load_mat(path):
    mat = _schema().data()
    with gzip.open(path, "rb") as handle:
        mat.ParseFromString(handle.read())
    return mat


def save_mat(mat, path):
    with gzip.open(path, "wb") as handle:
        handle.write(mat.SerializeToString())


def zero_root_branch_length(newick):
    """Set the root's own branch length to 0, leaving the rest untouched.

    Branch lengths in a MAT newick are mutation counts, so emptying the root's
    mutation list without this leaves the tree asserting a root branch of, say,
    113 mutations that the file no longer contains.
    """
    text = newick.rstrip()
    terminator = ""
    if text.endswith(";"):
        text, terminator = text[:-1], ";"
    close = text.rfind(")")
    if close == -1:
        return newick
    head, tail = text[:close + 1], text[close + 1:]
    if ":" in tail:
        tail = tail[:tail.rindex(":")] + ":0"
    return head + tail + terminator


def root_branch_length(newick):
    """The root's own branch length, or None if it carries none.

    Branch lengths in a MAT newick are mutation counts, so this is a second,
    independent encoding of len(node_mutations[0].mutation) -- which is what
    makes it useful for checking that node_mutations[0] really is the root.
    """
    text = newick.rstrip().rstrip(";")
    close = text.rfind(")")
    if close == -1:
        return None
    tail = text[close + 1:]
    if ":" not in tail:
        return None
    try:
        return int(tail.rsplit(":", 1)[1])
    except ValueError:
        return None


def root_sequence(reference, root_mutations):
    """Apply the root's mutation list to the reference, 1-based positions."""
    seq = list(reference)
    for mutation in root_mutations:
        index = mutation.position - 1
        if not 0 <= index < len(seq):
            logger.error(
                f"root mutation at position {mutation.position} is outside the "
                f"{len(seq)}-base reference"
            )
            sys.exit(1)
        # par_nuc is what usher believes sits at this position in the parent,
        # and the root's parent is the reference. Disagreement means the tree
        # and the reference are not the pair we think they are.
        expected = NUCLEOTIDES[mutation.par_nuc]
        if seq[index] != expected:
            logger.error(
                f"root mutation at position {mutation.position} expects parent "
                f"base {expected}, but the reference has {seq[index]}"
            )
            sys.exit(1)
        # mut_nuc is a repeated field so the schema can carry an ambiguous
        # state. Everything below reads mut_nuc[0], which would silently pick
        # one of several possibilities. usher has not written an ambiguous root
        # mutation for any of this pipeline's 16 combinations, so refuse rather
        # than guess if that ever changes.
        if len(mutation.mut_nuc) != 1:
            logger.error(
                f"root mutation at position {mutation.position} carries "
                f"{len(mutation.mut_nuc)} alternate bases; this script assumes "
                "an unambiguous state and will not choose between them"
            )
            sys.exit(1)
        seq[index] = NUCLEOTIDES[mutation.mut_nuc[0]]
    return "".join(seq)


def rebase_onto_root(node_mutations):
    """Move the origin onto node_mutations[0], in place. Returns (moved, count).

    `moved` is {position: allele} wherever the new origin differs from the old
    one, taken from the root's own mutation list. Every mutation in the tree at
    one of those positions has its ref_nuc repointed, because ref_nuc names the
    origin's base and the origin just changed. par_nuc and mut_nuc are left
    alone: they record a parent-to-child change, which is what it is regardless
    of what the origin is. That asymmetry is the whole reason this is
    bookkeeping rather than a change of content.

    Kept separate from main() so it can be tested without the compiled schema,
    which lives in envs/taxonium; it needs only objects with .mutation lists.
    """
    root = node_mutations[0]
    moved = {m.position: m.mut_nuc[0] for m in root.mutation}
    repointed = 0
    for node in node_mutations:
        for mutation in node.mutation:
            if mutation.position in moved:
                mutation.ref_nuc = moved[mutation.position]
                repointed += 1
    del root.mutation[:]
    return moved, repointed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-mat", required=True)
    parser.add_argument("--reference", required=True,
                        help="curated_reference.fasta the MAT was built against")
    parser.add_argument("--output-mat", required=True)
    parser.add_argument("--output-fasta", required=True,
                        help="the root sequence, the tree's new origin")
    parser.add_argument("--name", default="root",
                        help="FASTA record name for the root sequence")
    args = parser.parse_args()

    reference = read_reference(args.reference)
    mat = load_mat(args.input_mat)
    if not mat.node_mutations:
        logger.error(f"{args.input_mat} has no node mutation lists")
        sys.exit(1)

    root = mat.node_mutations[0]

    # node_mutations is written in newick DFS order, so index 0 is the root.
    # That is usher's convention, not something the schema enforces, and every
    # position this script writes depends on it. The newick states the root's
    # mutation count independently, as its branch length, so the two encodings
    # must agree -- verified to hold across all 16 combinations.
    declared = root_branch_length(mat.newick)
    if declared is not None and declared != len(root.mutation):
        logger.error(
            f"node_mutations[0] carries {len(root.mutation)} mutations but the "
            f"newick gives the root a branch length of {declared}; it may not "
            "be the root"
        )
        sys.exit(1)

    logger.info(f"Root carries {len(root.mutation)} mutations vs the reference")

    # Validates the root's mutations against the reference before anything is
    # mutated in place, so a mismatched pair fails without a partial rewrite.
    sequence = root_sequence(reference, root.mutation)
    moved, repointed = rebase_onto_root(mat.node_mutations)
    mat.newick = zero_root_branch_length(mat.newick)

    save_mat(mat, args.output_mat)
    with open(args.output_fasta, "w") as handle:
        handle.write(f">{args.name}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i:i + 60] + "\n")

    logger.info(
        f"Rebased onto the root: emptied its mutation list, repointed "
        f"{repointed} ref_nuc annotations at {len(moved)} positions"
    )
    logger.info(f"Wrote {len(sequence)}-base root sequence to {args.output_fasta}")


if __name__ == "__main__":
    main()
