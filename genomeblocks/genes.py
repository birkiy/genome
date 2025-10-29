"""Transcript/Gene/Genes definitions and GTF parsing helpers."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union

# local imports delayed where necessary to avoid circular refs
from .locus import Exon, CDS, UTR, Locus

@dataclass
class Transcript(Locus):
    transcript_id: str = ""
    exons: Optional[List[Exon]] = field(default_factory=list)
    cds: Optional[List[CDS]] = field(default_factory=list)
    utr: Optional[List[UTR]] = field(default_factory=list)

    def __post_init__(s):
        # Loci import delayed to avoid circular import
        from .loci import Loci
        if s.exons is not None: s.exons = Loci(s.exons)
        if s.cds is not None: s.cds = Loci(s.cds)
        if s.utr is not None: s.utr = Loci(s.utr)

    def add_exon(s, e: Exon) -> None:
        s.exons.append(e)
    def add_cds(s, c: CDS) -> None:
        s.cds.append(c)
    def add_utr(s, u: UTR) -> None:
        s.utr.append(u)


@dataclass
class Gene(Locus):
    gene_id: str = ""
    gene_name: Optional[str] = ""
    gene_type: Optional[str] = ""
    transcripts: Optional[Dict[str, Transcript]] = field(default_factory=dict)

    def __post_init__(s):
        if s.transcripts is not None: s.transcripts = {k: v for k,v in s.transcripts.items()}
        s.tss = Locus(s.chrom, s.start, s.start+1) if s.strand == '+' else Locus(s.chrom, s.end, s.end-1)

    def add_transript(s, t_id, t: Transcript) -> None:
        s.transcripts[t_id] = t


class Genes(dict):
    """a dictionary of genes (lightweight container)

    Many original helper methods are attached as functions so we preserve names.
    """
    def __init__(self, *args, filename:Optional[str]=None, _promoter_r:int=1000, **kwargs):
        super().__init__(*args, **kwargs)
        self.filename = filename
        self._annot = None
        self._promoter_r = _promoter_r

    def _build_annot(s):
        # delay Loci import to avoid circular import
        from .loci import Loci
        s._annot = {
            'body' : Loci(s.values()),
            'prom' : Loci(s.get_tss().values()).slop(s._promoter_r).sort().merge(),
            'exon' : Loci(e for g in s.values() for t in g.transcripts.values() for e in t.exons).sort().merge(),
            'utr5' : Loci(u for g in s.values() for t in g.transcripts.values() for u in t.utr if u.type ==  "5'").sort().merge(),
            'utr3' : Loci(u for g in s.values() for t in g.transcripts.values() for u in t.utr if u.type ==  "3'").sort().merge()
        }
        return

    @property
    def annot(s):
        if s._annot is None: s._build_annot()
        return s._annot

    def table(self):
        lines = []
        if len(self) > 0:
            lines.append(f"{'Name':<20} {'Count':<10}")
            lines.append("-" * 30)
            lines.append(f"{'Transcripts':<20} {sum(len(g.transcripts) for g in self.values()):<10}")
            lines.append(f"{'Exons':<20} {sum(len(t.exons) for g in self.values() for t in g.transcripts.values()):<10}")
            lines.append(f"{'CDS':<20} {sum(len(t.cds) for g in self.values() for t in g.transcripts.values()):<10}")
            lines.append(f"{'UTR':<20} {sum(len(t.utr) for g in self.values() for t in g.transcripts.values()):<10}")
        return "\n".join(lines)

    def get_tss(s, gene_type: Optional[Union[str, List[str]]] = None):
        tss = {}
        for g in s.values():
            if gene_type is None or (isinstance(gene_type, str) and g.gene_type == gene_type) or (isinstance(gene_type, list) and g.gene_type in gene_type):
                tss[g.gene_name] = g.tss
        return tss


# Helper: parse attributes column

def _parse_attributes(attr_str: str) -> Dict[str, str]:
    attrs = {}
    for attr in attr_str.split(";"):
        if not attr.strip(): continue
        key_value = attr.strip().split(" ", 1)
        if len(key_value) == 2:
            key, value = key_value
            attrs[key] = value.strip('"')
    return attrs


# GTF/GFF make function (as classmethod)
def make(cls, filename, gene_name_key='gene_name', gene_type_key='gene_type', chr_map=None, promoter_r=1000):
    genes = Genes(filename=filename, _promoter_r=promoter_r)
    unmapped = []

    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            fields = line.strip().split("\t")
            if len(fields) == 9:
                chrom, source, feature_type, start, end, score, strand, phase, attributes = fields
            elif len(fields) == 8:
                chrom, source, feature_type, start, end, score, strand, attributes = fields

            if chr_map is not None and chrom in chr_map: chrom = chr_map[chrom]
            start, end = int(start), int(end)
            attrs = _parse_attributes(attributes)
            gene_id = attrs['gene_id']

            if feature_type == "gene":
                gene_name = attrs.get(gene_name_key, gene_id)
                gene_type = attrs.get(gene_type_key, None)
                genes[gene_id] = Gene(chrom, start, end, strand, gene_id=gene_id, gene_name=gene_name, gene_type=gene_type)

            elif feature_type == "transcript":
                t_id = attrs['transcript_id']
                t = Transcript(chrom, start, end, strand, t_id)
                if gene_id not in genes:
                    genes[gene_id] = Gene(chrom, start, end, strand, gene_id=gene_id)
                genes[gene_id].add_transript(t_id, t)

            elif feature_type == "exon":
                t_id = attrs['transcript_id']
                e_number = int(attrs.get('exon_number', 0))
                e = Exon(chrom, start, end, strand, exon_number=e_number)
                genes[gene_id].transcripts[t_id].add_exon(e)

            elif feature_type == "CDS":
                t_id = attrs['transcript_id']
                e_number = int(attrs.get('exon_number', 0))
                c = CDS(chrom, start, end, strand, exon_number=e_number)
                genes[gene_id].transcripts[t_id].add_cds(c)

            elif feature_type in ('five_prime_UTR','three_prime_UTR', 'UTR'):
                t_id = attrs['transcript_id']
                e_number = int(attrs.get('exon_number', 0))
                if feature_type == 'five_prime_UTR': utr_type = "5'"
                elif feature_type == 'three_prime_UTR': utr_type = "3'"
                else:
                    if strand == '+': utr_type = "5'" if end <= genes[gene_id].transcripts[t_id].cds[0].start else "3'"
                    if strand == '-': utr_type = "3'" if end <= genes[gene_id].transcripts[t_id].cds[0].start else "5'"
                u = UTR(chrom, start, end, strand, exon_number=e_number, type=utr_type)
                genes[gene_id].transcripts[t_id].add_utr(u)
            else:
                unmapped.append(feature_type)
    unmapped = list(set(unmapped))
    if len(unmapped) > 0: print('[INFO] Unmapped feature types:', ', '.join(unmapped))
    return genes


# Attach helper functions to Genes to preserve original API assignment
Genes.make = classmethod(make)


def annotations(s: Genes, L):
    import pandas as pd
    annot = s.annot
    labels = []
    for l in L:
        lab = "Intergenic"
        if any(True for _ in annot['prom'].overlaps(l)):
            lab = "Promoter-TSS"
        elif any(True for _ in annot['utr5'].overlaps(l)):
            lab = "5UTR"
        elif any(True for _ in annot['utr3'].overlaps(l)):
            lab = "3UTR"
        elif any(True for _ in annot['exon'].overlaps(l)):
            lab = "Exonic"
        elif any(True for _ in annot['body'].overlaps(l)):
            lab = "Intronic"
        labels.append({"uid": l.uid, "annotation": lab})
    return pd.DataFrame(labels)


def nearest_genes(s, loci):
    from .loci import Loci
    l_tss = Loci(s.get_tss().values()).slop(s._promoter_r)
    g_names = list(s.get_tss().keys())
    return loci.nearest(l_tss, o_names=g_names)


Genes.annotations = staticmethod(annotations)
Genes.nearest_genes = staticmethod(nearest_genes)
