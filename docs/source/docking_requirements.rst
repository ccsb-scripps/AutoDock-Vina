.. _docking_requirements:

Software requirements
=====================

For those tutorials, the ADFR software suite, providing a number of software tools for automated docking and peripheral tasks, and the Python package Raccoon-lite, for preparing macrocycle for example, are necessary.

ADFR software suite
-------------------

.. warning::

    macOS users please note that ADFR software suite is NOT working under Catalina (10.15). We will send email to the mailing list if and when ADFR software suite or an alternative will be available for this version of macOS. If you already are using Catalina, we recommend install VirtualBox and running ADFR software suite inside the virtual box.

The ADFR software suite was developed in the Sanner lab at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_ formerly known as the Molecular Graphics Laboratory (MGL) of The Scripps Research Institute for visualization and analysis of molecular structures. You can find more information about the ADFR software suite installation process here: `https://ccsb.scripps.edu/adfr/downloads <https://ccsb.scripps.edu/adfr/downloads/>`_. The current version contains the following tools for docking:
    
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

Meeko
-----

The Python package ``meeko`` is a new type of package developped in the Forli lab also at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_.  It provides tools covering other docking aspects not handled by the ADFR software suite. This package provides addtionnal tools for the following docking protocols:

    - Hydrated docking
    - Macrocycles

**Meeko installation using pip**:

.. note::

    When using ``pip``, it's good pratice to use a virtual environment and also the easiest solution. An example with the `Conda package manager <https://docs.conda.io/en/latest/>`_ is available further down.

.. warning::
    
    OpenBabel executable and library must be installed first, see `installation instruction <https://open-babel.readthedocs.io/en/latest/Installation/install.html#install-binaries>`_.

.. code-block:: bash
    
    $ pip install -U numpy openbabel meeko

If the installation was successful, you should now be able to access to the following command from your terminal by typing:

    - mk_prepare_ligand.py

**Meeko installation in a Conda environment**:

.. note::

    See instructions in :ref:`installation` on how to setup and create an dedicated ``Conda`` environment.

Type the following command to install ``NumPy``, ``OpenBabel`` and ``meeko``:

.. code-block:: bash
    
    $ conda activate vina
    $ conda install -c conda-forge numpy openbabel
    $ pip install meeko
