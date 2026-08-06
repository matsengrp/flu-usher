#!/bin/bash
# Archive a snapshot of key pipeline outputs under
# /fh/fast/matsen_e/shared/flu-usher-archive/YYMMDD/, preserving the
# in-repo path of each file. Records the current HEAD commit in the
# archive's README.md so each snapshot is tied back to the code that
# produced it.
# Usage: bash scripts/copy_to_shared.sh

set -euo pipefail
# An unmatched glob must expand to nothing rather than to a literal asterisk,
# so a missing directory fails here instead of handing cp a bogus path.
shopt -s nullglob

SRC="results"
CONFIG="config.yaml"
ARCHIVE_ROOT="/fh/fast/matsen_e/shared/flu-usher-archive"
DATE="$(date +%y%m%d)"
DEST="$ARCHIVE_ROOT/$DATE/results"
README="$ARCHIVE_ROOT/README.md"
COMMIT="$(git rev-parse HEAD)"

if [ -e "$ARCHIVE_ROOT/$DATE" ]; then
    echo "Error: $ARCHIVE_ROOT/$DATE already exists. Remove it first or wait until tomorrow." >&2
    exit 1
fi

# Derive the segment/subtype combinations from config.yaml rather than
# hardcoding them, so adding a subtype to the config also adds it to the
# archive. This mirrors how download_ref_seq.py builds the same list.
mapfile -t SEGMENT_SUBTYPES < <(python3 - "$CONFIG" <<'PY'
import sys
import yaml

with open(sys.argv[1]) as f:
    config = yaml.safe_load(f)

for subtype in config["ha_subtypes"]:
    print(f"HA/{subtype}")
for subtype in config["na_subtypes"]:
    print(f"NA/{subtype}")
for segment in config["segments"]:
    if segment not in ("HA", "NA"):
        print(f"{segment}/all")
PY
)

if [ "${#SEGMENT_SUBTYPES[@]}" -eq 0 ]; then
    echo "Error: no segment/subtype combinations derived from $CONFIG" >&2
    exit 1
fi

for ss in "${SEGMENT_SUBTYPES[@]}"; do
    # Copy newick tree, taxonium jsonl, and the curated MSA
    mkdir -p "$DEST/$ss"
    cp "$SRC/$ss/final_tree.nwk" "$DEST/$ss/"
    cp "$SRC/$ss/final_tree.jsonl.gz" "$DEST/$ss/"
    cp "$SRC/$ss/msa.fasta.xz" "$DEST/$ss/"

    # Copy all unaligned coding sequence files
    coding_seqs=("$SRC/$ss/unaligned_coding_seqs/"*.fasta.xz)
    if [ "${#coding_seqs[@]}" -eq 0 ]; then
        echo "Error: no *.fasta.xz under $SRC/$ss/unaligned_coding_seqs/" >&2
        exit 1
    fi
    mkdir -p "$DEST/$ss/unaligned_coding_seqs"
    cp "${coding_seqs[@]}" "$DEST/$ss/unaligned_coding_seqs/"

    # Copy the PastML per-node ancestral state reconstructions: the states
    # table and the rendered tree, which are the two declared rule outputs.
    for ancestral in host_ancestral subtype_ancestral; do
        mkdir -p "$DEST/$ss/$ancestral"
        cp "$SRC/$ss/$ancestral/combined_ancestral_states.tab" "$DEST/$ss/$ancestral/"
        cp "$SRC/$ss/$ancestral/${ancestral%_ancestral}_tree.html" "$DEST/$ss/$ancestral/"
    done

    echo "Copied $ss"
done

# Global metadata, the PastML annotation tables, and the input-data manifest
cp "$SRC/combined_metadata_augmented.csv" "$DEST/"
cp "$SRC/host_annotation.csv" "$DEST/"
cp "$SRC/subtype_annotation.csv" "$DEST/"
cp "$SRC/input_data_md5sums.txt" "$DEST/"

# Make everything group- and world-readable
chmod -R o+rX "$ARCHIVE_ROOT/$DATE"

# Record the commit for this snapshot in the archive's README
echo "- \`$DATE/\` — commit \`$COMMIT\`" >> "$README"

echo "Done. Files copied to $DEST"
echo "Recorded commit $COMMIT under $DATE/ in $README"
