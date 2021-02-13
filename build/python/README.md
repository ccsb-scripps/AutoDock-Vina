# AutoDock Vina - Python API

### Requirements

You need, at a minimum (requirements):
* Python (>=3.5)
* Numpy
* SWIG
* Boost-cpp
* Sphinx (documentation)
* Sphinx_rtd_theme (documentation)

### Installation (from source)

I highly recommand you to install the Anaconda distribution (https://www.continuum.io/downloads) if you want a clean python environnment with nearly all the prerequisites already installed. To install everything properly, you just have to do this:
```bash
$ conda create -n vina python=3
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

### Full documentation

The installation instructions, documentation and tutorials can be found on [readthedocs.org](https://autodock-vina.readthedocs.io/en/latest/).

### Citations
* Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-461.
