#!/bin/bash
# Archive a snapshot of key pipeline outputs under
# /fh/fast/matsen_e/shared/flu-usher-archive/YYMMDD/, preserving the
# in-repo path of each file. Records the current HEAD commit in the
# archive's README.md so each snapshot is tied back to the code that
# produced it.
# Usage: bash scripts/copy_to_shared.sh

set -euo pipefail

SRC="results"
ARCHIVE_ROOT="/fh/fast/matsen_e/shared/flu-usher-archive"
DATE="$(date +%y%m%d)"
DEST="$ARCHIVE_ROOT/$DATE/results"
README="$ARCHIVE_ROOT/README.md"
COMMIT="$(git rev-parse HEAD)"

if [ -e "$ARCHIVE_ROOT/$DATE" ]; then
    echo "Error: $ARCHIVE_ROOT/$DATE already exists. Remove it first or wait until tomorrow." >&2
    exit 1
fi

# All segment/subtype combinations from the pipeline
SEGMENT_SUBTYPES=(
    HA/H1 HA/H3 HA/H5 HA/H7 HA/H9
    NA/N1 NA/N2 NA/N6 NA/N8 NA/N9
    PB2/all PB1/all PA/all NP/all MP/all NS/all
)

for ss in "${SEGMENT_SUBTYPES[@]}"; do
    # Copy newick tree, taxonium jsonl, and the curated MSA
    mkdir -p "$DEST/$ss"
    cp "$SRC/$ss/final_tree.nwk" "$DEST/$ss/"
    cp "$SRC/$ss/final_tree.jsonl.gz" "$DEST/$ss/"
    cp "$SRC/$ss/msa.fasta.xz" "$DEST/$ss/"

    # Copy all unaligned coding sequence files
    mkdir -p "$DEST/$ss/unaligned_coding_seqs"
    cp "$SRC/$ss/unaligned_coding_seqs/"*.fasta.xz "$DEST/$ss/unaligned_coding_seqs/"

    echo "Copied $ss"
done

# Global metadata and the input-data provenance manifest
cp "$SRC/combined_metadata_augmented.csv" "$DEST/"
cp "$SRC/host_annotation.csv" "$DEST/"
cp "$SRC/input_data_md5sums.txt" "$DEST/"

# Make everything group- and world-readable
chmod -R o+rX "$ARCHIVE_ROOT/$DATE"

# Record the commit for this snapshot in the archive's README
echo "- \`$DATE/\` — commit \`$COMMIT\`" >> "$README"

echo "Done. Files copied to $DEST"
echo "Recorded commit $COMMIT under $DATE/ in $README"
