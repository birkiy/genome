Genomeblocks
======

Fluent building blocks for regulatory genomics.

- `Loci`: craft and manipulate candidate regulatory element (CRE) sets.
- `Architecture`: overlay chromatin contacts and mcool matrices on your loci.
- `Genes`: annotate CREs with promoter-anchored gene models.

Quick Start
-----------

```python
from genomeblocks import Architecture, Genes, Loci

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

Setup
-----

```bash
conda env create -f environment.yml
conda activate genomeblocks
```