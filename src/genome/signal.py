from __future__ import annotations
import time
import shutil
from multiprocessing import Process, Value, cpu_count
from multiprocessing.shared_memory import SharedMemory
from typing import List, Tuple, Sequence, Dict
from collections.abc import Sequence as SequenceABC

import numpy as np
import pyBigWig as bw
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import gridspec

from .loci import Loci


def plan_workers(n_tracks: int, n_loci: int, *, cores: int | None = None,
                max_bw_parallel: int = 6) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """
    Plan the distribution of work across multiple processes.
    
    Args:
        n_tracks: Number of bigwig tracks
        n_loci: Number of genomic loci
        cores: Number of CPU cores to use (default: all available)
        max_bw_parallel: Maximum number of bigwig files to process in parallel
    
    Returns:
        List of work chunks as [((track_start,track_end), (loci_start,loci_end)), ...]
    """
    c = cores or cpu_count()
    t_chunks = min(n_tracks, max_bw_parallel, c)
    t_step = (n_tracks + t_chunks - 1) // t_chunks
    l_per_t = max(1, c // t_chunks)
    l_step = (n_loci + l_per_t - 1) // l_per_t

    plan = []
    for ti in range(t_chunks):
        t_lo = ti * t_step
        t_hi = min(n_tracks, (ti + 1) * t_step)
        for li in range(l_per_t):
            l_lo = li * l_step
            l_hi = min(n_loci, (li + 1) * l_step)
            plan.append(((t_lo, t_hi), (l_lo, l_hi)))
    return plan


def _worker(
    t_b: Tuple[int, int],
    l_b: Tuple[int, int],
    bwp: Sequence[str],
    loci: Loci,
    n_bins: int,
    flank: int,
    agg: str,
    shm_name: str,
    shape: Tuple[int, int, int],
    dt: np.dtype,
    ctr: Value, # type: ignore
    span: bool = False,
):
    """Worker process for parallel signal extraction from bigwig files."""
    t_lo, t_hi = t_b
    l_lo, l_hi = l_b
    shm = SharedMemory(name=shm_name)
    cube = np.ndarray(shape, dtype=dt, buffer=shm.buf)

    hs = [bw.open(p) for p in bwp[t_lo:t_hi]]
    ch = hs[0].chroms()

    try:
        for r in range(l_lo, l_hi):
            loc = loci[r]
            chrom = loc.chrom
            if chrom not in ch:  # skip unknown chrom
                with ctr.get_lock():
                    ctr.value += (t_hi - t_lo)
                continue

            if span:
                L = max(0, loc.start)
                R = min(ch[chrom], loc.end)
                pre, post = 0, 0
                core = max(1, n_bins)
            else:
                c = (loc.start + loc.end) // 2
                L, R = c - flank, c + flank
                pre = max(0, -L)
                post = max(0, R - ch.get(chrom, 0))
                L = max(0, L)
                R = min(ch.get(chrom, 0), R)
                core = n_bins - (pre + post)

            if core > 0:
                for off, h in enumerate(hs):
                    t = t_lo + off
                    xs = h.stats(chrom, L, R, nBins=core, type=agg)
                    xs = np.nan_to_num(xs).astype(dt, copy=False)
                    if span:
                        arr = xs
                    else:
                        arr = np.concatenate((
                            np.zeros(pre, dtype=dt),
                            xs,
                            np.zeros(post, dtype=dt)
                        ))
                    cube[r, t, :len(arr)] = arr

            with ctr.get_lock():
                ctr.value += (t_hi - t_lo)
    finally:
        for h in hs:
            h.close()
        shm.close()


def signal(
    loci: Loci,
    bigwigs: Sequence[str],
    *,
    n_bins: int = 200,
    flank: int = 3_000,
    agg: str = "mean",
    dtype: str | np.dtype = np.float32,
    progress: bool = True,
    max_bw_parallel: int = 6,
    max_workers: int | None = None,
    span: bool = False,
) -> np.ndarray:
    """
    Extract signal from bigwig files for given genomic loci.
    
    Args:
        loci: Genomic loci to extract signal from
        bigwigs: List of bigwig file paths
        n_bins: Number of bins to divide each region into
        flank: Number of base pairs to include on each side of locus center
        agg: Aggregation method ('mean', 'max', 'min', etc.)
        dtype: Numpy dtype for the output array
        progress: Show progress bar
        max_bw_parallel: Maximum number of bigwig files to process in parallel
        max_workers: Maximum number of worker processes
        span: Use full region span instead of center±flank
    
    Returns:
        numpy array of shape (n_loci, n_tracks, n_bins)
    """
    n_loci, n_tracks = len(loci), len(bigwigs)
    if n_loci == 0:
        raise ValueError("No loci provided.")

    bytes_need = n_loci * n_tracks * n_bins * np.dtype(dtype).itemsize
    if bytes_need > 0.5 * shutil.disk_usage("/").free:
        raise MemoryError("Cube may exceed safe RAM; try disk-chunk mode.")

    shm = SharedMemory(create=True, size=bytes_need)
    cube = np.ndarray((n_loci, n_tracks, n_bins), dtype=dtype, buffer=shm.buf)
    cube.fill(0)

    plan = plan_workers(
        n_tracks,
        n_loci,
        cores=(max_workers or cpu_count()),
        max_bw_parallel=max_bw_parallel,
    )

    ctr = Value('i', 0)
    procs = []
    for t_b, l_b in plan:
        p = Process(target=_worker, args=(
            t_b, l_b, bigwigs, loci, n_bins, flank, agg,
            shm.name, cube.shape, dtype, ctr, span
        ))
        p.start()
        procs.append(p)

    if progress:
        tot = n_loci * n_tracks
        with tqdm(total=tot, dynamic_ncols=True) as bar:
            last = 0
            while any(p.is_alive() for p in procs):
                v = ctr.value
                bar.update(v - last)
                last = v
                time.sleep(0.2)
            bar.update(ctr.value - last)

    for p in procs:
        p.join()
        if p.exitcode != 0:
            shm.close()
            shm.unlink()
            raise RuntimeError(f"worker {p.pid} exit {p.exitcode}")

    out = cube.copy()
    shm.close()
    shm.unlink()
    return out

def _bcast(x, n, name):
    """Convert scalar/str to list [x]*n; sequence of length n passes through"""
    if isinstance(x, str) or not isinstance(x, SequenceABC):
        return [x] * n
    if len(x) != n:
        raise ValueError(f"{name} must have length {n}")
    return list(x)


def plot_heatmap(
    loci: Loci,
    S: np.ndarray,
    *,
    sets: List[str] | None = None,
    samples: List[str] | None = None,
    colors: Dict[str, tuple] | None = None,
    ymax: float = 10,
    ymin: float = 0,
    height: int = 3000,
    cmap: str = "Blues",
    vmax: float = 10,
    profile: bool = True,
    no_sort: bool = False,
    dpi: int = 100,
):
    """
    Plot heatmap of genomic signals.
    
    Args:
        loci: Genomic loci with .attr['Set'] membership
        S: Signal array of shape (regions × tracks × bins)
        sets: List of set names
        samples: List of sample names
        colors: Dictionary mapping set names to colors
        ymax: Maximum y-axis value for profile plots
        ymin: Minimum y-axis value for profile plots
        height: Region height in base pairs
        cmap: Colormap for heatmaps
        vmax: Maximum value for heatmap color scaling
        profile: Include average profile above heatmaps
        no_sort: Don't sort regions by mean signal
        dpi: Figure DPI
    
    Returns:
        matplotlib Figure
    """
    if sets is None:
        sets = sorted({loc.attr.get("Set", "regions") for loc in loci})
    if samples is None:
        samples = [f"track_{i}" for i in range(S.shape[1])]

    n = len(samples)
    cmaps = _bcast(cmap, n, "cmap")
    vms = _bcast(vmax, n, "vmax")
    ys = _bcast(ymax, n, "ymax")
    yl = _bcast(ymin, n, "ymin")

    if colors is None:
        colors = {k: plt.get_cmap('tab10')(i) for i, k in enumerate(sets)}

    gidx = [np.array([loc.attr.get("Set") == s for loc in loci]) for s in sets]
    order = (np.arange(len(loci)) if no_sort
             else np.argsort(S.mean(axis=(1, 2)))[::-1])
    S_ = S[order]
    g_ = [idx[order] for idx in gidx]
    nb = S.shape[-1]

    rows = ([sum(idx.sum() for idx in g_) // 4] if profile else []) + [idx.sum() for idx in g_]
    fig = plt.figure(figsize=(3 * n, 10), dpi=dpi)
    gs = gridspec.GridSpec(len(rows), n, height_ratios=rows)
    plt.subplots_adjust(hspace=0.05, wspace=0.3)

    for i, s in enumerate(samples):
        if profile:
            ax = fig.add_subplot(gs[0, i])
            for j, idx in enumerate(g_):
                ax.plot(S_[idx, i, :].mean(0), color=colors[sets[j]], lw=2)
            ax.set_ylim(yl[i], ys[i])
            ax.set_xticks([])
            if i == 0:
                ax.set_ylabel("signal")
            ax.set_title(s)

        for j, idx in enumerate(g_):
            ax = fig.add_subplot(gs[j + (1 if profile else 0), i])
            ax.imshow(S_[idx, i, :], aspect='auto', cmap=cmaps[i], vmin=0, vmax=vms[i])
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_ylabel(sets[j], rotation=0, ha='right', va='center')
            if j == len(g_) - 1:
                ax.set_xticks([0, nb // 2, nb])
                kb = round(height / 1000, 1)
                ax.set_xticklabels([f"-{kb}kb", "center", f"+{kb}kb"])
    return fig


def plot_profiles(
    loci: Loci,
    S: np.ndarray,
    *,
    sets: List[str] | None = None,
    colors: Dict[str, tuple] | None = None,
    ylim: float | None = None,
    dpi: int = 100,
    height: int = 3000,
):
    """
    Plot average signal profiles by set.
    
    Args:
        S: Signal array of shape (regions × tracks × bins)
        loci: Genomic loci with .attr['Set'] membership
        sets: List of set names
        colors: Dictionary mapping set names to colors
        ylim: Y-axis limit
        dpi: Figure DPI
        height: Region height in base pairs
    
    Returns:
        matplotlib Figure
    """
    if sets is None:
        sets = sorted({loc.attr.get("Set", "regions") for loc in loci})
    if colors is None:
        colors = {k: plt.get_cmap('tab10')(i) for i, k in enumerate(sets)}

    gidx = [np.array([loc.attr.get("Set") == s for loc in loci]) for s in sets]
    fig, axs = plt.subplots(1, len(sets), figsize=(len(sets) * 3, 3), dpi=dpi, squeeze=False)
    axs = axs[0]

    for ax, idx, lab in zip(axs, gidx, sets):
        m = S[idx, 0, :].mean(0)
        ax.plot(m, color=colors[lab], lw=2)
        ax.set_title(lab)
        ax.set_xticks([0, S.shape[-1] // 2, S.shape[-1]])
        kb = height // 1000
        ax.set_xticklabels([f"-{kb}kb", "center", f"+{kb}kb"])
        ax.set_ylim(0, (ylim if ylim is not None else np.percentile(m, 99)))
        ax.set_yticks([])

    axs[0].set_ylabel("signal")
    return fig


Loci.signal = signal

Loci.plot_heatmap = plot_heatmap
Loci.plot_profiles = plot_profiles