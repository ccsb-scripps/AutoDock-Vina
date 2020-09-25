[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![PyPI version fury.io](https://img.shields.io/badge/version-1.2.0-green.svg)](https://pypi.python.org/pypi/ansicolortags/) 

# AutoDock Vina
Molecular docking and virtual screening program

## Python API

### Prerequisites

You need, at a minimum (requirements):
* Python (=3.7)
* OpenBabel
* SWIG
* Boost-cpp
* Sphinx (documentation)
* Sphinx_rtd_theme (documentation)

### Installation

I highly recommand you to install the Anaconda distribution (https://www.continuum.io/downloads) if you want a clean python environnment with nearly all the prerequisites already installed. To install everything properly, you just have to do this:
```bash
$ conda create -n vina python=3.7
$ conda activate vina
$ conda install -c conda-forge openbabel swig boost-cpp sphinx sphinx_rtd_theme
```

Finally, we can install the `Vina` package
```bash
$ git clone https://github.com/ccsb-scripps/AutoDock-Vina
$ git checkout boost-python
$ cd AutoDock-Vina/build/python
$ python setup.py build install
```

### Documentation

Build python documentation with Sphinx
```bash
$ cd build/python/docs
$ make html
```

### Quick tutorial
```python
#!/usr/bin/env python

from vina import Vina


v = Vina(exhaustiveness=32)

v.set_receptor("protein.pdbqt")
v.set_ligand('ligand.pdbqt')

v.compute_vina_maps([0., 0., 0.], [30, 30, 30], 0.375)
print(v.score())
print(v.optimize())
v.dock()
v.write_docking_results("docking_results.pdbqt")
```
