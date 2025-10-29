
genome
======

This repository contains the `genome` Python package (refactored from the previous `locus` monolith) and a conda recipe scaffold for building and distributing the package via conda.

Installation (recommended)

1) Using conda (recommended for heavy binary deps):

```bash
# create an environment and install runtime deps from conda-forge
conda create -n genome python=3.10 -c conda-forge numpy pandas -y
conda activate genome

# install this package from source (local)
python -m pip install .
```

2) Installing from PyPI/TestPyPI is intentionally not the default here; this repo is being prepared for conda-forge distribution under the name `genome`.

Building the conda package locally

```bash
# install conda-build if needed
conda install -c conda-forge conda-build
# run conda-build from repo root
conda-build conda/recipe
```

Notes
- The conda recipe `conda/recipe/meta.yaml` is a scaffold; update `about.home`, `extra.recipe-maintainers` and `source` to point to your GitHub repository or tarball before opening a conda-forge feedstock PR.
- If you want me to open a feedstock PR on conda-forge, I can prepare the feedstock files; you'll still need to create the feedstock repository or grant rights where necessary.

