"""locus2 package public API with lazy imports.

This module exposes the original top-level names but delays importing
their implementation modules until the attribute is accessed. This keeps
``import locus2`` fast and avoids requiring heavy optional dependencies
to be installed unless code paths that need them are used.

Public names preserved: Locus, Exon, CDS, UTR, Features, Tags, Transcript,
Gene, Genes, Loci, Architecture, make_genome, scan_motifs
"""

from importlib import import_module
from types import ModuleType
from typing import Dict

# Mapping of public name -> (module_path, attribute_name)
# attribute_name can be same as public name or a different symbol in module
_EXPORTS: Dict[str, tuple[str, str]] = {
	# locus
	'Locus': ('.locus', 'Locus'),
	'Exon': ('.locus', 'Exon'),
	'CDS': ('.locus', 'CDS'),
	'UTR': ('.locus', 'UTR'),
	# features
	'Tags': ('.tags', 'Tags'),
	# genes
	'Transcript': ('.genes', 'Transcript'),
	'Gene': ('.genes', 'Gene'),
	'Genes': ('.genes', 'Genes'),
	# loci
	'Loci': ('.loci', 'Loci'),
	# architecture
	'Architecture': ('.architecture', 'Architecture'),
	# motifs
	'make_genome': ('.motifs', 'make_genome'),
	'scan_motifs': ('.motifs', 'scan_motifs'),
}





__all__ = sorted(list(_EXPORTS.keys()))


def _import_from(module_name: str, attr: str):
	"""Import `attr` from package-relative `module_name` and cache in globals."""
	# module_name may be relative (starts with '.')
	mod = import_module(f"{__package__}{module_name}")
	value = getattr(mod, attr)
	globals()[attr] = value
	return value


def __getattr__(name: str):
	"""Lazy attribute loader for public names."""
	if name in _EXPORTS:
		module_name, attr = _EXPORTS[name]
		return _import_from(module_name, attr)
	raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
	return sorted(list(globals().keys()) + __all__)

