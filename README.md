genome
======

Fluent building blocks for regulatory genomics.

- `Loci`: craft and manipulate candidate regulatory element (CRE) sets.
- `Architecture`: overlay chromatin contacts and mcool matrices on your loci.
- `Genes`: annotate CREs with promoter-anchored gene models.
- `genome.signal`: extract bigWig signal, make heatmaps, draw profiles.

Quick Start
-----------

```python
from genome import Architecture, Genes, Loci

cre = (Loci.make(data['ATAC-bed'])
           .slop(100)
           .sort()
           .merge())
se = cre.intersect(Loci.make(data['H3K27ac-SE']))

arch = (Architecture.make(cre, data['RNAP-loops'], r=2500)
                       .add_mcool(cre, data['RNAP-mcool'], resolution=5000)
                       .normalize(cre))

genes = Genes.make(data['Gencode-GTF'], promoter_r=1000)
counts = genes.annotations(se & cre).groupby('annotation').size()
```

Signal Utilities
----------------

```python
from genome import Loci
from genome.signal import signal, plot_heatmap, plot_profiles

loci = Loci.make("path/to/regions.bed")
signals = loci.signal(bigwigs=["sample1.bw", "sample2.bw"], n_bins=200, flank=3000)

loci.plot_heatmap(S=signals, samples=["Sample1", "Sample2"], profile=True)
lcoi.plot_profiles(S=signals)
```

Setup
-----

```bash
# managed environment
conda env create -f conda/environment.yml
conda activate genome-dev
python -m pip install .

# manual environment
conda create -n genome python=3.10 -c conda-forge numpy pandas scipy cgranges pyranges tqdm -y
conda activate genome
python -m pip install .
```

Packaging
---------

```bash
conda install -c conda-forge conda-build
conda-build conda/recipe
```

Notes
- Update metadata in `conda/recipe/meta.yaml` before releasing.
- Optional extras (`cooler`, `lightmotif`, `graph-tool`) are commented in `conda/environment.yml`.
