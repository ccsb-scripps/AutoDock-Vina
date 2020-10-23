.. _zinc_docking:

Docking with zinc metalloproteins
=================================

Introduction
------------

Zinc is present in a wide variety of proteins and is important in the metabolism of most organisms. Zinc metalloenzymes are therapeutically relevant targets in diseases such as cancer, heart disease, bacterial infection, and Alzheimer’s disease. In most cases a drug molecule targeting such enzymes establishes an interaction that coordinates with the zinc ion. Thus, accurate prediction of the interaction of ligands with zinc is an important aspect of computational docking and virtual screening against zinc containing proteins. The AutoDock force field was extended to include a specialized potential describing the interactions of zinc-coordinating ligands. This potential describes both the energetic and geometric components of the interaction. The new force field, named AutoDock4Zn, was calibrated on a data set of 292 crystal complexes containing zinc. Redocking experiments show that the force field provides significant improvement in performance in both free energy of binding estimation as well as in root-mean-square deviation from the crystal structure pose.

**Please cite this paper, if you are using this protocol in your work**:
    - Santos-Martins, D., Forli, S., Ramos, M. J., & Olson, A. J. (2014). AutoDock4Zn: an improved AutoDock force field for small-molecule docking to zinc metalloproteins. Journal of chemical information and modeling, 54(8), 2371-2379.

**Software requirements**

The following software must be already installed:
AutoDock version 4.2.6 or later:
http://autodock.scripps.edu/downloads
MGLTools version 1.5.7 or later 
http://mgltools.scripps.edu/downloads
Instructions refers to terminal  commands, therefore basic knowledge of UNIX  or Windows command lineusage is advisable. 

The variable $MGLROOT refers to your MGLTools installation directory, e.g.:
/usr/local/mgltools-1.5.7/  (Linux/Mac)
"C:\Program Files (x86)\MGLTools-1.5.7\”  (Windows)

WINDOWS USERS:  Use the “python.exe” command instead of “pythonsh” for all commands. The path separator in Windows is the backslash (“\”), although Unix paths below can be used by using quotes aroundthem (i.e., cd “AutoDock4Zn_toolbox/”).

**Required files**

For convenience, have the following files in the working directory:

- Data files
    - protein.pdb and ligand.mol2 (AutoDock4Zn_toolbox/Tutorial)
    - AD4Zn.dat (AutoDock4Zn_toolbox)
- Scripts
    - prepare_gpf4zn.py (AutoDock4Zn_toolbox/)
    - prepare_dpf42.py ($MGLROOT/MGLToolsPckgs/AutoDockTools/Utilities24)
    - zinc_pseudo.py (AutoDock4Zn_toolbox/)
    - prepare_receptor4.py ($MGLROOT/MGLToolsPckgs/AutoDockTools/Utilities24)
    - prepare_ligand4.py ($MGLROOT/MGLToolsPckgs/AutoDockTools/Utilities24)
- Binaries
    - autodock4 (version: 4.2.6 or newer)
    - autogrid4(version: 4.2.6 or newer)
    - pythonsh ($MGLROOT/bin)  (Linux/Mac)
    - python.exe(“C:\Program Files (x86)\MGLTools-1.5.7\”)(Windows)

Summary checklist
The following steps differ from the typical AutoDock protocol:   
Run zinc_pseudo.py   
Generate receptor PDBQT filewith prepare_gpf4zn.py   
Use of AD4Zn.dat in the parameter files   
Use of autogrid4.2.5.x.20131125

1. Preparing the receptor
-------------------------

Add polar hydrogens, gasteiger charges and set atom types:   

.. code-block:: bash

    $ pythonsh prepare_receptor4.py -r protein.pdb -o protein.pdbqt   

It is common to fine tune the receptor by modeling missing atoms, or testing different combinations of conformations and protonation states. Such tasks are not discussed in this tutorial.

Add Tetrahedral Zinc Pseudo Atoms (TZ) to the receptor

TZ atoms represent the preferred position for tetrahedral coordination by the ligand.

.. code-block:: bash

    $ pythonsh zinc_pseudo.py -r protein.pdbqt -o protein_tz.pdbqt

Step 3 - Use the modified forcefield (AD4Zn.dat)The AutoDock4Zn forcefield is mostly defined by non bonded pairwise potentials which are written to the grid parameter file (\*.gpf) in the form of npb_r_eps keywords. The file AD4Zn.dat includes the definition of the TZ atom type for the AutoDock forcefield. The keyword parameter_file in the gpf specifies AD4Zn.dat as the forcefield to be used, so AutoGrid requires a local copy of it in the working directory. Alternatively, the keyword parameter_file in the .gpf can point to the full or relative path where AD4Zn.dat is located.

2. Preparing the ligand
-----------------------

.. code-block:: bash

    $ pythonsh prepare_ligand4.py -l ligand.mol2 -o ligand.pdbqt -A  hydrogens

3. Generating affinity maps
---------------------------

The preparation script will be executed to generate the GPF to configure the grid calculation:

.. code-block:: bash

    $ pythonsh prepare_gpf4zn.py -l ligand.pdbqt -r protein_tz.pdbqt \
    -o protein_tz.gpf  -p npts=40,30,50 -p gridcenter=18,134,-1 \
    –p parameter_file=AD4Zn.dat   

The -p flag is used to set the box center (gridcenter) and size (npts) along with the parameter_file specific for this case.

The code that supports user defined pairwise potentials (npb_r_eps keyword) from the .gpf file was restored from an old version of the software and added to the modified AutoGrid4 binary provided in this tutorial. This changes will be included in next release of the standard AutoDock binaries. 

.. code-block:: bash

    $ autogrid4 -p protein_tz.gpf -o protein_tz.glg

At this stage, all forcefield information has been encoded in the maps (\*.map), and the remaining steps are the same as in the standard AutoDock protocol.

4. Running AutoDock Vina
------------------------

.. code-block:: bash

    $ vina -p ligand_protein_tz.dpf -o ligand_protein_tz.dlg

5. Results
----------
