[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![PyPI version fury.io](https://img.shields.io/badge/version-1.2.0-green.svg)](https://pypi.python.org/pypi/ansicolortags/) 

# AutoDock Vina
Molecular docking and virtual screening program

## Python API

### Prerequisites

You need, at a minimum (requirements):
* Python (=3.7)
* Numpy
* SWIG
* Boost-cpp
* Sphinx (documentation)
* Sphinx_rtd_theme (documentation)

### Installation

I highly recommand you to install the Anaconda distribution (https://www.continuum.io/downloads) if you want a clean python environnment with nearly all the prerequisites already installed. To install everything properly, you just have to do this:
```bash
$ conda create -n vina python=3.7
$ conda activate vina
$ conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme
```

Finally, we can install the `Vina` package
```bash
$ git clone https://github.com/ccsb-scripps/AutoDock-Vina
$ cd AutoDock-Vina
$ git checkout boost-python
$ cd build/python
$ python setup.py build install
```

### Documentation

Build python documentation with Sphinx
```bash
$ cd docs
$ make html
# Open index.html located in docs/build/html
```

### Quick tutorial
```python
#!/usr/bin/env python

from vina import Vina


v = Vina()

v.set_receptor(rigid_pdbqt_filename="protein.pdbqt")
v.set_ligand_from_file('ligand.pdbqt')

v.compute_vina_maps(center=[0., 0., 0.], box_size=[30, 30, 30])
print(v.score())
print(v.optimize())
v.dock(exhaustiveness=32)
v.write_poses(pdbqt_filename="docking_results.pdbqt")
```
