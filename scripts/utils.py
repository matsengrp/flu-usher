"""
Utility functions shared across flu-usher pipeline scripts.
"""

import gzip
import logging
import lzma
import re

DIGITS = "0123456789"

# IUPAC nucleotide codes plus the gap character, in both cases. matUtils only
# ever writes ACGT for this pipeline's data -- 2,628,520 tokens checked across
# all 16 combinations, no other character appears -- but a MAT can carry an
# ambiguous state, so accepting the codes costs nothing and rejecting anything
# outside the set is what makes the check worth having.
BASES = frozenset("ACGTUNRYSWKMBDHV" + "acgtunryswkmbdhv" + "-")


def setup_logging(name=None):
    """
    Set up logging configuration

    Args:
        name: Logger name, conventionally the caller's __name__. Defaults to
              this module's name, which is what callers that predate the
              argument get.

    Returns:
        logging.Logger: Configured logger
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(name if name is not None else __name__)


def parse_mutation(token):
    """Split a matUtils mutation token into (parent base, position, new base).

    Tokens look like "A33T": the base in the parent, a 1-based position, and
    the base in this node. Note the first character is the *parent's* base, not
    the reference's -- they coincide only for mutations on the root.

    Raises ValueError rather than returning something wrong, because the two
    callers used to slice this format by hand and disagreed about how much to
    check; that format has already produced one production bug (issue #49).
    """
    if len(token) < 3:
        raise ValueError(f"malformed mutation token {token!r}: too short")
    parent, position, mutant = token[0], token[1:-1], token[-1]
    # Both checks are deliberately charset-explicit rather than str.isdigit() /
    # str.isalpha(), which accept any Unicode digit or letter: '2'.isdigit() is
    # True for the superscript form that int() then rejects, and 'A33O'.isalpha()
    # passes with a Greek omega that would be written straight into a FASTA.
    # Failing at the parse boundary is the whole point of this function.
    if not position or position.strip(DIGITS) != "":
        raise ValueError(f"malformed mutation token {token!r}: bad position")
    for base in (parent, mutant):
        if base not in BASES:
            raise ValueError(f"malformed mutation token {token!r}: bad base {base!r}")
    return parent, int(position), mutant


def iter_path_mutations(path_field):
    """Yield (parent base, position, new base) from a `matUtils extract -S` path.

    The field is space-separated `node:MUT,MUT` chunks in root-to-tip order.
    Mutations come out in that order, so a later one at the same position
    supersedes an earlier one, and a mutation followed by its back-mutation
    cancels.
    """
    for chunk in path_field.split():
        _, _, mutations = chunk.partition(":")
        for token in mutations.split(","):
            if token:
                yield parse_mutation(token)


def read_reference(path):
    """Return the first sequence in a FASTA file as one uppercased string.

    Deliberately not Bio.SeqIO: this is called from rules whose environments
    carry no biopython, and the inputs are single-record reference FASTAs.
    """
    parts = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if parts:
                    break
                continue
            parts.append(line.strip())
    return "".join(parts).upper()


def open_sequence_file(path, mode='rt'):
    """
    Open a possibly-compressed sequence file in text mode.

    Compression is selected from the file extension, so the same call site
    serves plain, gzip- and xz-compressed paths. Pipeline sequence files are
    xz-compressed; gzip is supported because some inputs arrive that way.

    Args:
        path: Path to the file
        mode: Text-mode flag, 'rt' to read or 'wt' to write

    Returns:
        An open text-mode file handle. Use as a context manager.
    """
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    if path.endswith('.xz'):
        return lzma.open(path, mode)
    # Plain files take the single-letter text mode: 'rt' -> 'r', 'wt' -> 'w'.
    return open(path, mode[0])


def sanitize_id(seq_id):
    """
    Sanitize sequence ID by removing problematic characters.
    This is used to ensure sequence IDs are compatible with downstream tools.

    Args:
        seq_id: Original sequence ID

    Returns:
        Sanitized sequence ID
    """
    sanitized = seq_id
    for char in ['[', ']', '(', ')', ':', ';', ',', "'", '.']:
        sanitized = sanitized.replace(char, '')
    return sanitized


def extract_attribute_value(attributes, attr_name):
    """
    Extract a specific attribute value from GFF attributes string

    Args:
        attributes: GFF attributes string
        attr_name: Name of the attribute to extract (e.g., 'ID', 'Name')

    Returns:
        str: Attribute value or None if not found
    """
    attr_pattern = f"{attr_name}="
    if attr_pattern in attributes:
        start_idx = attributes.find(attr_pattern) + len(attr_pattern)
        end_idx = attributes.find(';', start_idx)
        if end_idx == -1:
            end_idx = len(attributes)
        return attributes[start_idx:end_idx]
    return None


def get_coding_region_coords(gff_file):
    """
    Extract coding region coordinates from GFF file.
    Finds all gene/CDS features and returns the min start and max end.

    Args:
        gff_file: Path to GFF file

    Returns:
        tuple: (min_start, max_end) of coding region
    """
    logger = logging.getLogger(__name__)

    # Regular expression for parsing GFF
    re_gff = re.compile(r'([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)')

    logger.info(f"Extracting coding region coordinates from {gff_file}")

    starts = []
    ends = []

    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            match = re_gff.match(line)
            if not match:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = match.groups()

            # Look for gene or CDS features
            if feature_type.lower() in ('gene', 'cds'):
                starts.append(int(start))
                ends.append(int(end))

    if not starts:
        raise ValueError(f"Could not find any gene or CDS features in {gff_file}")

    min_start = min(starts)
    max_end = max(ends)

    logger.info(f"Coding region spans {min_start}-{max_end}")
    return min_start, max_end


def extract_all_genes_and_cds(gff_file):
    """
    Extract information for all genes and CDS features from a GFF file

    Args:
        gff_file: Path to GFF file

    Returns:
        list: List of dictionaries containing gene/CDS information
        tuple: (min_start, max_end) coordinates for all features
    """
    logger = logging.getLogger(__name__)

    # Regular expressions for parsing GFF
    re_gff = re.compile(r'([^\t]+)\t([^\t]+)\t([^\t]+)\t(\d+)\t(\d+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)')

    logger.info(f"Extracting all genes and CDS features from {gff_file}")

    features = []

    with open(gff_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            match = re_gff.match(line)
            if not match:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = match.groups()

            # Convert to integers
            start, end = int(start), int(end)

            # Strip newline character from attributes
            attributes = attributes.strip()

            # Look for gene or CDS features
            if feature_type.lower() in ('gene', 'cds'):
                # Extract ID and Name from attributes
                feature_id = extract_attribute_value(attributes, 'ID')
                feature_name = extract_attribute_value(attributes, 'Name')

                # Use ID as fallback for name if Name is not present
                if not feature_name:
                    feature_name = feature_id

                # Use feature_type as fallback if neither ID nor Name is present
                if not feature_name:
                    feature_name = f"{feature_type}_feature"

                feature_info = {
                    'id': feature_id,
                    'name': feature_name,
                    'type': feature_type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'phase': phase,
                    'attributes': attributes,
                    'original_attributes': attributes  # Keep original for reference
                }
                features.append(feature_info)
                logger.info(f"Found {feature_type} ID='{feature_id}' Name='{feature_name}' at position {start}-{end}")

    if not features:
        raise ValueError(f"Could not find any genes or CDS features in {gff_file}")

    # Calculate overall min and max coordinates
    min_start = min(f['start'] for f in features)
    max_end = max(f['end'] for f in features)

    logger.info(f"Found {len(features)} features spanning {min_start}-{max_end}")
    return features, (min_start, max_end)


def group_cds_by_gene(features):
    """
    Group CDS features by their protein ID for handling spliced genes.

    Spliced genes (like NEP) have multiple CDS features with the same protein_id.
    This function groups them together so they can be concatenated properly.

    Args:
        features: List of feature dicts from extract_all_genes_and_cds()

    Returns:
        dict: {key: {'gene_name': str, 'protein_id': str, 'cds_list': [sorted CDS features]}}
              where key is protein_id if available, otherwise gene_name
    """
    logger = logging.getLogger(__name__)

    # Filter to only CDS features
    cds_features = [f for f in features if f['type'].lower() == 'cds']

    genes = {}

    for cds in cds_features:
        # Extract gene name and protein_id from attributes
        gene_name = extract_attribute_value(cds['attributes'], 'gene') or cds['name']
        protein_id = extract_attribute_value(cds['attributes'], 'protein_id')

        # Use protein_id as key to group spliced CDS (same protein_id = same gene)
        # Fall back to gene_name if protein_id is not available
        key = protein_id if protein_id else gene_name

        if key not in genes:
            genes[key] = {
                'gene_name': gene_name,
                'protein_id': protein_id,
                'cds_list': []
            }
        genes[key]['cds_list'].append(cds)

    # Sort CDS features by start position for each gene
    # This ensures spliced genes are concatenated in the correct order
    for gene_data in genes.values():
        gene_data['cds_list'].sort(key=lambda x: x['start'])

    logger.debug(f"Grouped {len(cds_features)} CDS features into {len(genes)} genes")

    return genes