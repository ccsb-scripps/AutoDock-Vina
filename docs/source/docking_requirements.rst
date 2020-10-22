Software requirements
=====================

For those tutorials, the ADFR software suite, providing a number of software tools for automated docking and peripheral tasks, and the Python package Raccoon-lite, for preparing macrocycle for example, are necessary.

ADFR software suite
-------------------

.. warning::

	MAC users please note that ADFR software suite is NOT working under the Catalina OS. We strongly advise to refrain from upgrading your OS to Catalina. We will send email to the mailing list if and when ADFR software suite or an alternative will be available for this version of Mac OS. If you already are using Catalina, we recommend install VirtualBox and running ADFR software suite inside the virtual box.

Thew ADFR software suite was developed in the Sanner lab at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_ formerly known as the Molecular Graphics Laboratory (MGL) of The Scripps Research Institute for visualization and analysis of molecular structures. You can find more information about the ADFR software suite installation process here: `https://ccsb.scripps.edu/adfr/downloads <https://ccsb.scripps.edu/adfr/downloads/>`_. The current version contains the following tools for docking:
    
    - ADFR v1.2 and associate scripts
    - AGFR v1.2
    - AutoSite v1.0 and v1.1
    - ADCP v1.0
    - AutoGrid4.2
    - prepare_ligand
    - prepare_receptor

Moreover, the ADFR software suite provides a number of software tools for automated docking and peripheral tasks. These tools are implemented using the Python, C++ and C programming languages and a re-usable component philosophy. To avoid Python packages mismatches we opted to shift ADFR suite with a self-contain Python interpreter that is isolated from the default Python interpreter installed on your computer (except for Windows installations). Details about the implementation and packages provided by the ADFR software suite can be found here: `https://ccsb.scripps.edu/adfr/implementation <https://ccsb.scripps.edu/adfr/implementation/>`_

**Citations**:
	
	- Zhang, Y., Forli, S., Omelchenko, A., & Sanner, M. F. (2019). AutoGridFR: Improvements on AutoDock Affinity Maps and Associated Software Tools. Journal of Computational Chemistry, 40(32), 2882-2886.
	- Zhang, Y., & Sanner, M. F. (2019). AutoDock CrankPep: combining folding and docking to predict protein–peptide complexes. Bioinformatics, 35(24), 5121-5127.
	- Ravindranath, P. A., & Sanner, M. F. (2016). AutoSite: an automated approach for pseudo-ligands prediction—from ligand-binding sites identification to predicting key ligand atoms. Bioinformatics, 32(20), 3142-3149.
	- Ravindranath, P. A., Forli, S., Goodsell, D. S., Olson, A. J., & Sanner, M. F. (2015). AutoDockFR: advances in protein-ligand docking with explicitly specified binding site flexibility. PLoS computational biology, 11(12), e1004586.
	- Zhao, Y., Stoffler, D., & Sanner, M. (2006). Hierarchical and multi-resolution representation of protein flexibility. Bioinformatics, 22(22), 2768-2774.

Raccoon-lite
------------

.. note::

	Raccoon-lite is not related to Raccoon, contained in MGLTools, and does not serve the same purpose. More information about Raccoon, a Virtual Screenings graphical interface for AutoDock and AutoDock Vina, can be found here: `https://ccsb.scripps.edu/autodock/raccoon <https://ccsb.scripps.edu/autodock/raccoon/>`_

The Python package Raccoon-lite is a new type of package developped in the Forli lab also at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_.  It provides tools covering other docking aspects not handled by the ADFR software suite. This package provides addtionnal tools for the following docking protocols:

	- Hydrated docking
	- Macrocycles
	- Docking with Zinc metal

This new tool can be easily installed with all its depencies using `pip` or `conda` package managers. 

**Raccoon-lite installation using Conda**:

The `Conda package manager <https://docs.conda.io/en/latest/>`_ is included as part of the Anaconda Python distribution, which can be download from `https://docs.continuum.io/anaconda/install <https://docs.continuum.io/anaconda/install/>`_. This is a Python distribution specially designed for scientific applications, with many of the most popular scientific packages preinstalled. Alternatively, you can use `Miniconda <https://conda.pydata.org/miniconda.html>`_, which includes only Python itself, plus the Conda package manager.

1. Begin by installing the most recent 64 bit, Python 3.x version of either Anaconda or Miniconda
2. (Optional, but highly suggested) If you want, you can create a dedicated environment for `raccoon-lite` and `AutoDock Vina` packages:

.. code-block:: bash

	conda create -n vina python=3
	conda activate vina
	conda install -c conda-forge numpy openbabel swig boost-cpp sphinx sphinx_rtd_theme
	conda install -c ccsb-scripps autodock-vina

3. And type the following command

.. code-block:: bash

	conda install -c ccsb-scripps raccoon-lite

**Raccoon-lite installation using pip**:

.. code-block:: bash

	pip install raccoon-lite

If the installation was successful, you should now be able to access to the following command/tools from your terminal by typing:

	- rc_prepare_macrocycle.py
	- rc_wet.py
	- rc_mapwater.py
	- rc_dry.py
	- rc_prepare_gfp4zn.py
	- rc_zinc_pseudo.py
