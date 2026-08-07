"""Compatibility import for the annotated CNV VCF parser.

Generic file-format parsing is owned by :mod:`cosap_formats`. New consumers should
import from that package directly; this module preserves the historical COSAP path.
"""

from cosap_formats import AnnotatedCnvVcfParser

__all__ = ["AnnotatedCnvVcfParser"]
