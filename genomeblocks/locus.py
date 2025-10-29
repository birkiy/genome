"""Core locus dataclasses: Locus, Exon, CDS, UTR.

Minimal, copied from original core; no function renames.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict

@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    strand: str = "."

    @property
    def uid(s) -> str: return f"{s.chrom}:{s.start}-{s.end}({s.strand})"

    @property
    def length(l) -> int: return l.end - l.start

    @property
    def center(l) -> int: return (l.start + l.end) // 2

    def __eq__(s, o: object) -> bool: return isinstance(o, Locus) and s.uid == o.uid
    def __ne__(s, o: object) -> bool: return isinstance(o, Locus) and s.uid != o.uid
    def __lt__(s, o: "Locus") -> bool: return NotImplemented if not isinstance(o, Locus) else \
        s.start < o.start if s.chrom == o.chrom else s.chrom < o.chrom
    def __le__(s, o: "Locus") -> bool: return NotImplemented if not isinstance(o, Locus) else \
        s.start <= o.start if s.chrom == o.chrom else s.chrom < o.chrom
    def __gt__(s, o: "Locus") -> bool: return NotImplemented if not isinstance(o, Locus) else \
        s.start > o.start if s.chrom == o.chrom else s.chrom < o.chrom
    def __ge__(s, o: "Locus") -> bool: return NotImplemented if not isinstance(o, Locus) else \
        s.start >= o.start if s.chrom == o.chrom else s.chrom < o.chrom

    def distance_to(s, o: "Locus") -> int: return abs(s.center - o.center) if s.chrom == o.chrom else NotImplemented
    def overlaps(s, o: "Locus") -> bool: return (s.chrom == o.chrom) and not (s.end <= o.start or s.start >= o.end)

    def sequence(s, genome: Dict, r=None) -> str:
        if r is not None:
            return genome[s.chrom][s.center-r:s.center+r]
        else:
            return genome[s.chrom][s.start:s.end]

    def copy(s) -> "Locus": return Locus(s.chrom, s.start, s.end, s.strand)


@dataclass
class Exon(Locus):
    exon_number: int = 0


@dataclass
class CDS(Exon):
    pass


@dataclass
class UTR(Exon):
    type: str = None  # "5'" or "3'"
