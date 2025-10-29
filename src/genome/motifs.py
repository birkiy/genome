"""Motif scanning utilities."""
from typing import Dict


def _parse_fasta(path):
    import gzip
    name = ""
    seq = list()
    if path.endswith(".gz"):
        f = gzip.open(path, "rt")
    else:
        f = open(path, "r")
    for _, l in enumerate(f):
        l = l.rstrip("\n")
        if l[0] == ">":
            if len(seq) == 0:
                name = l[1:].split(" ")[0]
                continue
            yield (name, "".join(seq))
            seq = list()
            name = l[1:].split(" ")[0]
        else:
            seq.append(l)
    yield (name, "".join(seq))


def make_genome(path):
    return {k:v for k,v in _parse_fasta(path)}


def scan_motifs(s, genome, motif_path, motif_format='jaspar', r=250, threshold=13.0, norm=True, verbose=True):
    import lightmotif
    from tqdm import tqdm
    # ensure genome is a dict
    if not isinstance(genome, dict): genome = make_genome(genome)
    n_motif = len([_ for _ in lightmotif.load(motif_path, format=motif_format)])
    M = {}
    for m in (tqdm(lightmotif.load(motif_path, format=motif_format), total=n_motif) if verbose else lightmotif.load(motif_path, format=motif_format)):
        M[m.name] = 0
        pssm = m.counts.normalize(0.1).log_odds()
        for l in s:
            seq = lightmotif.stripe(l.sequence(genome, r=r).upper())
            for _ in lightmotif.scan(pssm, seq, threshold=threshold):
                M[m.name] += 1
        if norm: M[m.name] = M[m.name] / len(m.counts)
    return M
