"""Architecture graph utilities (graph-tool wrapper)."""
from typing import Optional

import graph_tool as gt

class Architecture(gt.Graph):

    def __init__(s, name: Optional[str]=None):
        super().__init__(directed=False); s._name = name or "Architecture"
        s.vp.uid = s.new_vertex_property("string")
        s.index = {}
        s.ep.w = s.new_edge_property("float")
        s.ep.n = s.new_edge_property("float")
        s.ep.d = s.new_edge_property("float")
        print("[INFO] Initialized an empty Architecture graph. 🏗️")

    def _add_vertex(s, uid):
        if uid not in s.index:
            v = s.add_vertex()
            s.vp.uid[v] = uid
            s.index[uid] = v
        return s.index[uid]

    @property
    def n_loci(s) -> int: return s.num_vertices()
    @property
    def n_links(s) -> int: return s.num_edges()

    def __str__(s) -> str:
        eprops = ", ".join(s.ep.keys())
        vprops = ", ".join(s.vp.keys())
        return f"Architecture(name='{s._name}', loci={s.n_loci}, links={s.n_links}, vertex_props=[{vprops}], edge_props=[{eprops}])"
    __repr__ = __str__

    def __len__(s) -> int: return s.num_vertices()
    def __contains__(s, uid: str) -> bool: return uid in s.index

    def __getitem__(s, key):
        if isinstance(key, str):
            if key not in s.index: raise KeyError(f"Locus UID '{key}' not found in the graph.")
            v = s.index[key]
            neighbors = {}
            for neighbor_v in v.all_neighbors():
                neighbor_uid = s.vp.uid[neighbor_v]
                edge = s.edge(v, neighbor_v)
                neighbors[neighbor_uid] = s.ep.w[edge]
            return neighbors
        if isinstance(key, tuple) and len(key) == 2:
            uid1, uid2 = key
            if uid1 not in s.index or uid2 not in s.index: raise KeyError(f"One or both UIDs ('{uid1}', '{uid2}') not found.")
            v1, v2 = s.index[uid1], s.index[uid2]
            edge = s.edge(v1, v2)
            if edge is None: return None
            return {name: prop[edge] for name, prop in s.ep.items()}
        raise TypeError("Key must be a Locus UID (string) or a tuple of two UIDs.")

    def copy(s) -> "Architecture":
        if s.is_directed(): raise NotImplementedError("Copy is not implemented for directed Architecture graphs.")
        new_arch = s.__class__(name=s._name)
        for v_old in s.vertices():
            uid = s.vp.uid[v_old]
            new_arch._add_vertex(uid)
        for e_old in s.edges():
            uid1 = s.vp.uid[e_old.source()]
            uid2 = s.vp.uid[e_old.target()]
            v1_new = new_arch.index[uid1]
            v2_new = new_arch.index[uid2]
            e_new = new_arch.add_edge(v1_new, v2_new)
            for key, prop_map in s.ep.items():
                value = prop_map[e_old]
                new_arch.ep[key][e_new] = value
        print(f"[INFO] Created a copy of '{s._name}'.  copier  copier 🐑")
        return new_arch

# attach methods defined elsewhere (keeps same top-level API in core shim)

def make(cls, loci, bedpe: str, *, name: str="Skeleton", r: int=2500, dmax=1e7, verbose: bool=True):
    from tqdm import tqdm
    G = cls(name=name)
    mapped_loops = 0
    total_loops = 0
    with open(bedpe) as f:
        for line in tqdm(f, desc='[INFO] Building Architecture from loops'):
            if line.startswith('#') or not line.strip(): continue
            total_loops += 1
            fields = line.strip().split()
            if len(fields) < 6: continue
            try:
                chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
                chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
            except ValueError:
                continue
            mid1 = (start1 + end1) // 2; mid2 = (start2 + end2) // 2
            if abs(mid1 - mid2) > dmax: continue
            a1 = [loci[j] for *_, j in loci.cgr.overlap(chrom1, mid1 - r, mid1 + r)]
            a2 = [loci[j] for *_, j in loci.cgr.overlap(chrom2, mid2 - r, mid2 + r)]
            if not a1 or not a2: continue
            mapped_loops += 1
            for locus1 in a1:
                for locus2 in a2:
                    if locus1.uid == locus2.uid: continue
                    v1 = G._add_vertex(locus1.uid)
                    v2 = G._add_vertex(locus2.uid)
                    G.add_edge(v1, v2)
    if verbose:
        pct_mapped = 100 * mapped_loops / max(total_loops, 1)
        print(f"[INFO] {total_loops} loops | {mapped_loops} mapped ({pct_mapped:.1f}%) | "
              f"loci={G.n_loci}, links={G.n_links}")
    return G


def add_mcool(s: Architecture, loci, mcool: str, *, resolution: Optional[int]=None, name: str="w", verbose: bool=True) -> Architecture:
    import cooler
    from tqdm import tqdm
    import numpy as np
    from collections import defaultdict
    print(f"[INFO] Adding '{name}' weights to the graph. 🏗️")
    uri = f"{mcool}::resolutions/{resolution}" if resolution else mcool
    clr = cooler.Cooler(uri)
    bins = clr.bins()[:][['chrom', 'start', 'end']].reset_index()
    from .loci import Loci
    from .locus import Locus
    bins_l = Loci(Locus(row[1], row[2], row[3]) for row in bins.itertuples(index=False))
    near = loci.nearest(bins_l)
    uid_to_bin = dict(zip(near['Name'], near['Name_b'].map(bins_l.uids)))
    bin_pair_edge_counts = defaultdict(int)
    for edge in s.edges():
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        bin1 = uid_to_bin.get(uid1)
        bin2 = uid_to_bin.get(uid2)
        if bin1 is None or bin2 is None: continue
        if bin1 > bin2: bin1, bin2 = bin2, bin1
        bin_pair_edge_counts[(bin1, bin2)] += 1
    pixels = clr.pixels()[:].set_index(['bin1_id', 'bin2_id'])['count']
    edges_set = 0
    s.ep[name] = s.new_edge_property("float")
    for edge in tqdm(s.edges(), desc='[INFO] Assigning distributed weights to edges', total=s.n_links):
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        bin1 = uid_to_bin.get(uid1)
        bin2 = uid_to_bin.get(uid2)
        if bin1 is None or bin2 is None: continue
        bin1, bin2 = (min(bin1, bin2), max(bin1, bin2))
        try:
            total_count = float(pixels[bin1][bin2].sum())
            if total_count > 0:
                num_edges_in_pair = bin_pair_edge_counts.get((bin1, bin2), 1)
                distributed_weight = total_count / num_edges_in_pair
                s.ep[name][edge] += distributed_weight
                edges_set += 1
        except KeyError:
            s.ep[name][edge] += 0
        except TypeError:
            print(pixels[bin1][bin2])
            break
    if verbose:
        print(f"[INFO] Set distributed weights for {edges_set}/{s.n_links} edges from cooler. [{name}]")
    return s


def _pl_model(x, C, alpha):
    return C * (x**(-alpha))


def _pl_expect(x, y):
    import numpy as np
    from scipy.optimize import curve_fit
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    if not m.any():
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}
    try:
        b, a = np.polyfit(np.log10(x[m]), np.log10(y[m]), 1)
        guess = [10.0**a, -b]
        popt, _ = curve_fit(_pl_model, x[m], y[m], p0=guess)
        C_fit, alpha_fit = popt
        y_expected = _pl_model(x, C_fit, alpha_fit)
        return y_expected, {"alpha": float(alpha_fit), "C": float(C_fit)}
    except (RuntimeError, np.linalg.LinAlgError, ValueError):
        return np.full_like(y, np.nan), {"alpha": np.nan, "C": np.nan}


def normalize(s: Architecture, loci, *, source: str = "w", name: str = "n", verbose: bool = True) -> Architecture:
    import numpy as np
    from tqdm import tqdm
    if verbose:
        print(f"[INFO] Normalizing '{source}' by power-law expectation. Storing in '{name}'. 📏")
    if source not in s.ep:
        raise ValueError(f"Source edge property '{source}' not found in the graph.")
    if name not in s.ep:
        s.ep[name] = s.new_edge_property("float")
    edge_list = list(s.edges())
    distances = np.zeros(len(edge_list), dtype=float)
    raw_weights = np.zeros(len(edge_list), dtype=float)
    for i, edge in tqdm(enumerate(edge_list),desc='[INFO] Calculating distances for edges', total=len(edge_list)):
        v1, v2 = edge.source(), edge.target()
        uid1, uid2 = s.vp.uid[v1], s.vp.uid[v2]
        dist = loci[uid1].distance_to(loci[uid2])
        s.ep.d[edge] = dist
        raw_weights[i] = s.ep[source][edge]
        distances[i] = dist
    expected, fit_params = _pl_expect(distances, raw_weights)
    if verbose:
        print(f"[INFO] Power-law fit complete: alpha={fit_params['alpha']:.3f}, C={fit_params['C']:.3e}")
    norm_weights = raw_weights / np.maximum(expected, 1e-12)
    for i, e in enumerate(s.edges()):
        s.ep[name][e] = norm_weights[i]
    if verbose:
        print(f"[INFO] Set O/E weights for {s.n_links} intra-chromosomal edges to `ep.{name}`.")
    return s

# attach these utilities to the Architecture class to preserve the old API
Architecture.make = make
Architecture.add_mcool = add_mcool
Architecture.normalize = normalize
