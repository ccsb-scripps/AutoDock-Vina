.. _docking_requirements:

Software requirements
=====================

In addition to AutoDockTools, the long standing GUI for docking calculation setup and management, the following helper tools are developed at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_ and can be used for docking input preparation. 

Python package Meeko
-----

The Python package ``meeko`` is package recently developped in the Forli lab also at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_. As showcased in the Colab examples, `Meeko <https://github.com/forlilab/Meeko>`_ provides commandline scripts for ligand preparation, receptor preparation and other essential tools for the following docking protocols:

    - Docking with flexible macrocycles
    - Flexible docking
    - Covalent docking
    - Reactive docking
    - Hydrated docking

**Meeko installation using pip**:

.. note::

    When using ``pip``, it's good practice to do so in a virtual (conda, mamba, etc.) environment instead of the base environment. An example with the `Conda package manager <https://docs.conda.io/en/latest/>`_ is available further down.

.. code-block:: bash
    
    $ pip install -U numpy scipy rdkit vina meeko prody

If the installation was successful, you should now be able to access to the following command from your terminal by typing:

    - mk_prepare_ligand.py
    - mk_prepare_receptor.py
    - mk_export.py

**Meeko installation in a Conda environment**:

.. note::

    See instructions in :ref:`installation` on how to setup and create an dedicated ``Conda`` environment.

Type the following command to install ``Meeko`` and (optionally) ``ProDy``:

.. code-block:: bash
    
    $ conda activate vina
    $ conda install python=3.10    # for ProDy interoperability
    $ conda install -c conda-forge numpy scipy rdkit vina meeko
    $ pip install prody


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

