"""
Utility functions shared across flu-usher pipeline scripts.
"""

import logging
import re


def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


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