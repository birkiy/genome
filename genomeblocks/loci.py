"""Loci container and interval operations."""
from typing import Iterable, Optional, Union, List, Dict
from collections import defaultdict

# Local imports are delayed to avoid circular imports
from .locus import Locus
import cgranges as cr

class Loci(list):

    def __init__(self, iterable: Iterable[Locus] = (), *, filename: Optional[str] = None):
        super().__init__(iterable)
        self.filename = filename
        self._uids: Optional[dict] = None
        self._cgr: Optional[cr.cgranges] = None

    def _build_cgr(self) -> None:
        idx = cr.cgranges()
        for i, l in enumerate(self):
            idx.add(l.chrom, int(l.start), int(l.end), i)
        idx.index()
        self._cgr = idx

    @property
    def cgr(self) -> cr.cgranges:
        if self._cgr is None: self._build_cgr()
        return self._cgr


    @property
    def uids(self) -> dict:
        if self._uids is None: self._uids = {l.uid: i for i, l in enumerate(self)}
        return self._uids

    def to_pyranges(s, names=None):
        import pyranges as pr
        df = s.to_frame().rename(columns={'Chr': 'Chromosome'})
        if names is None: names = [l.uid for l in s]
        df['Name'] = names
        return pr.PyRanges(df)

    def to_frame(s, names=None):
        import pandas as pd
        if names is None: names = [l.uid for l in s]
        return pd.DataFrame([
            {"Chr": l.chrom, "Start": l.start, "End": l.end, 'Strand': l.strand, 'Name': names[i]}
            for i, l in enumerate(s)
        ])

    def to_bed(s, path=None):
        lines = []
        for l in s:
            lines.append(f"{l.chrom}\t{l.start}\t{l.end}\t{l.uid}\t0\t{l.strand}")
        content = "\n".join(lines)
        if path is not None:
            with open(path, 'w') as f:
                f.write(content)


    def __and__(s, o: "Loci") -> "Loci": return s.intersect(o)
    def __add__(s, o: "Loci") -> "Loci": return Loci([*s, *o])
    def __sub__(s, o: "Loci") -> "Loci": return s.difference(o)
    def __truediv__(s, o: "Loci") -> "Loci": return s.difference(o)
    def __or__(s, o: "Loci") -> "Loci": return Loci([*s, *o])
    def __xor__(s, o: "Loci") -> "Loci": return (s - o) + (o - s)

    def copy(self) -> "Loci":
        return Loci(list(self), filename=self.filename)

    def __str__(self):
        return f"Loci(n={len(self)})"
    __repr__ = __str__

    def __getstate__(self):
        s = self.__dict__.copy()
        s.pop("_cgr", None)
        return s

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._cgr = None

    def __getitem__(self, key: Union[int, str, slice]):
        if isinstance(key, int) or isinstance(key, slice):
            return super().__getitem__(key)
        elif isinstance(key, str):
            try:
                index = self.uids[key]
                return super().__getitem__(index)
            except KeyError:
                raise KeyError(f"No Locus with UID '{key}' found.")
        else:
            raise TypeError(f"Invalid key type: {type(key)}. Expected int, str, or slice.")


# helper functions attached to Loci (kept as in original)

def _auto_format(p: str) -> str:
    if p.endswith(".gtf"): return "Features"
    if p.endswith((".gff", ".gff3")): return "Features"
    if p.endswith(".narrowPeak"): return "Loci"
    return "Loci"

@classmethod
def make(cls, filename: str, filetype: Optional[str] = None) -> Loci:
    if filetype is None: filetype = _auto_format(filename)
    if filetype == "Features": raise NotImplementedError("GTF/GFF parsing not implemented. Please use Features.make() instead.")
    elif filetype != 'Loci': raise ValueError(f"Unsupported file type for Loci!")
    L = Loci(filename=filename)
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            strand = fields[5] if len(fields) > 5 else "."
            L.append(Locus(chrom, start, end, strand))
    return L


def subloci(s, uids: List[str]):
    sub = Loci(s[s.uids[u]] for u in uids)
    return sub

Loci.make = make
Loci.subloci = subloci


def overlaps(s, key:Union[Locus, str], start:Optional[int]=None, end:Optional[int]=None):
    if isinstance(key, Locus):
        start = key.start
        end = key.end
        key = key.chrom
    return Loci(s[i] for *_ , i in s.cgr.overlap(key, start, end))


def intersect(s, o: Loci) -> Loci:
    return Loci([l for l in s if any(True for *_ , _ in o.cgr.overlap(l.chrom, l.start, l.end))])


def difference(s, o: Loci) -> Loci:
    return Loci([l for l in s if not any(True for *_ , _ in o.cgr.overlap(l.chrom, l.start, l.end))])


def slop(s, n:int) -> Loci:
    return Loci([Locus(l.chrom, max(0, l.start - n), l.end + n, l.strand) for l in s])


def sort(s) -> Loci:
    out = Loci(sorted(s))
    out.uids
    return out


def merge(s) -> Loci:
    if not s: return Loci()
    sorted_loci = sort(s)
    merged = [sorted_loci[0].copy()]
    for l in sorted_loci[1:]:
        last = merged[-1]
        if last.chrom != l.chrom: merged.append(l.copy()); continue
        if last.overlaps(l) or last.end == l.start:
            last.end = max(last.end, l.end)
        else:
            merged.append(l.copy())
    return Loci(merged)


def nearest(s: Loci, o: Loci, s_names=None, o_names=None):
    pr_s = s.to_pyranges(s_names)
    pr_o = o.to_pyranges(o_names)
    return pr_s.nearest(pr_o).df.rename(columns={'Chromosome': 'Chr'})


from .tags import Tags
from .genes import Genes


def map(
    s: Loci,
    mapping: Union[Dict[str, "Loci"], Genes],
) -> Tags:
    """
    Convenience wrapper that materialises tag masks for the provided mapping.

    Args:
        mapping: Dict of label -> Loci overlaps to evaluate, or a `Genes` object.

    Returns:
        Tags database populated with the requested annotations.
    """
    tags = Tags.make(s)
    tags.add(mapping)
    return tags


def tag(s: Loci, o: Loci, tag: str):
    return Tags.make(s).add({tag: o})


Loci.overlaps = overlaps
Loci.intersect = intersect
Loci.difference = difference
Loci.slop = slop
Loci.sort = sort
Loci.merge = merge
Loci.nearest = nearest
Loci.tag = tag
Loci.map = map
