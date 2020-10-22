Basic docking
=============

Let's start with our first example of docking, where the typical usage pattern would be to dock a single molecule into a rigid receptor. In this example we will dock the approved anticancer drug imatinib (Gleevec; PDB entry 1iep) in the structure of c-Abl using AutoDock Vina. The target for this protocol is the kinase domain of the proto-oncogene tyrosine protein kinase c-Abl. The protein is an important target for cancer chemotherapy—in particular, the treatment of chronic myelogenous leukemia.

.. note::
	This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the raccoon-lite Python package.

1. Preparing the receptor
-------------------------

During this step, we will create a PDBQT file of our receptor that will contain hydrogen atoms, as well as partial charges. As a prerequisite, a receptor coordinate file must contain all hydrogen atoms. If hydrogen atoms are absent in the protein structure file, you can use the `-A "hydrogens"` flag to the following command for preparing the receptor PDBQT file. Many tools exist to add missing hydrogen atoms to a protein, one popular choice would be to use `REDUCE <http://kinemage.biochem.duke.edu/software/reduce.php>`_. If you are using experimental structures (for instance, from the PDB), use a text editor to remove water, ligands, cofactors, ions deemed unnecessary. The file `1iep_receptorH.pdb` is provided for this protocol, and it includes receptor coordinates that are taken from PDB entry 1iep.

.. code-block:: bash

	prepare_receptor.py -r 1iep_receptorH.pdb -o 1iep_receptor.pdbqt


2. Preparing the ligand
-----------------------

This step is very similar to the previous step. We will also create a PDBQT file from a ligand molecule file (in MOL2 or PDB format). As well as for the receptor, the coordinate set must includes all hydrogen atoms according to the choosen protonation state. This may be obtained in a variety of ways, including with experimental coordinates from the `PDB <https://www.rcsb.org>`_ or `Cambridge Crystallographic Database <http://www.ccdc.cam.ac.uk>`_. The file `1iep_ligandH.pdb` is provided for use as a tutorial for this protocol. This file includes ligand coordinates taken from PDB entry 1iep, to which all hydrogen atoms have been added and manually adjusted to the known protonation state.

.. code-block:: bash

	prepare_ligand.py -r 1iep_ligandH.pdb -o 1iep_ligand.pdbqt


3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

Now, we have to define the grid space for the docking, typically, a 3D box around a the potential binding site of a receptor. During this step, we will create the input file for `AutoGrid4`, which will create an affinity map file for each atom types. The grid parameter file specifies an AutoGrid calculation, including the size and location of the grid, the atom types that will be used, the coordinate file for the rigid receptor, and other parameters for calculation of the grids.

.. code-block:: console
	:caption: Content of the grid parameter file (**1iep.gpf**) for the receptor c-Abl (**1iep_receptor.pdbqt**)

	npts 40 46 40                        # num.grid points in xyz
	gridfld 1iep.maps.fld                # grid_data_file
	spacing 0.375                        # spacing(A)
	receptor_types A C HD N OA SA        # receptor atom types
	ligand_types A C NA OA N HD          # ligand atom types
	receptor 1iep_receptor.pdbqt         # macromolecule
	gridcenter 15.19 53.903 16.917       # xyz-coordinates or auto
	smooth 0.5                           # store minimum energy w/in rad(A)
	map 1iep.A.map                       # atom-specific affinity map
	map 1iep.C.map                       # atom-specific affinity map
	map 1iep.NA.map                      # atom-specific affinity map
	map 1iep.OA.map                      # atom-specific affinity map
	map 1iep.N.map                       # atom-specific affinity map
	map 1iep.HD.map                      # atom-specific affinity map
	elecmap 1iep.e.map                   # electrostatic potential map
	dsolvmap 1iep.d.map                  # desolvation potential map
	dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

We created the grid parameter file in the previous step, and now we can use `AutoGrid` to generate the different map files that will be used for the molecular docking.

.. code-block:: bash

	autogrid4 -p 1iep.gpf -l 1iep.glg

4. Running AutoDock Vina
------------------------

The imatinib ligand used in this protocol is challenging, and Vina will occasionally not find the correct pose with the default parameters. Vina provides a parameter called `Exhaustiveness` to change the amount of computational effort used during a docking experiment. The default exhaustiveness value is 8; increasing this to about 32 will give a more consistent docking result. At this point of the tutorial, you have the choice to decide to run the molecular docking using either the AutoDock forcefield (requires affinity maps) or using the Vina/Vinardo forcefield.

4.a. Using AutoDock forcefield
______________________________

When using the AutoDock forcefield, you only need to provide the affinity maps and the ligand, while specifying that the forcefield used will be AutoDock4 using the option `--scoring ad4`.

.. code-block:: bash

	vina  --ligand 1iep_ligand.pdbqt --maps 1iep --scoring ad4 \
	      --exhaustiveness 32 --out 1iep_ligand_ad4_out.pdbqt

4.b. Using Vina forcefield
__________________________

Contrary to AutoDock4, you don't need to precalculate the affinity grid maps with `autogrid` when using the Vina forcefield. AutoDock Vina computes those maps internally before the docking. However, you still need to specify the center and dimensions (in Angstrom) of the grid space, as well as the receptor.

.. code-block:: console
	:caption: Content of the config file (**box.txt**) for AutoDock Vina

	center_x = 15.19
	center_y = 53.903
	center_z = 16.917
	size_x = 15.0
	size_y = 17.25
	size_z = 15.0

.. code-block:: bash

	vina --receptor 1iep_receptor.pdbqt --ligand 1iep_ligand.pdbqt --config box.txt \
	     --exhaustiveness=32 --out 1iep_ligand_vina_out.pdbqt

Running AutoDock Vina will write a docked coordinate file `1iep_ligand_out.pdbqt` and also present docking information to the terminal window.

5. Results
----------

The predicted free energy of binding should be about −13 kcal mol –1 for poses that are similar to the crystallographic pose. With exhaustiveness set to 32, Vina will most often give a single docked pose with this energy. With the lower default exhaustiveness, several poses flipped end to end, with less favorable energy, may be reported.
