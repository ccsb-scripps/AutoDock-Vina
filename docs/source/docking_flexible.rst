.. _flexible_docking:

Flexible docking
================

The lack of receptor flexibility is arguably the greatest limitation in these types of docking methods. However, AutoDock Vina allows some limited flexibility of selected receptor side chains. In this tutorial, we will describe the cross-docking of imatinib to c-Abl in PDB entry 1fpu, treating Thr315 as flexible. 

.. tip::

	Some systems might experience larger motions of loops or domains, the relaxed complex method has shown success by sampling a variety of receptor conformations using molecular dynamics and then performing docking simulations on these snapshots.

	- Lin, J. H., Perryman, A. L., Schames, J. R., & McCammon, J. A. (2003). The relaxed complex method: Accommodating receptor flexibility for drug design with an improved scoring scheme. Biopolymers: Original Research on Biomolecules, 68(1), 47-62.

.. note::
	This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the raccoon-lite Python package.

1. Preparing the flexible receptor
----------------------------------

During this step, we are going to split the receptor coordinates into two PDBQT files: one for the rigid portion and one for the flexible side chains. As with the rigid docking tutorial, the method requires a receptor coordinate file that includes all hydrogen atoms. 

.. code-block:: bash
	
	$ prepare_receptor -r 1fpu_receptorH.pdb -o 1fpu_receptor.pdbqt
	$ prepare_flexreceptor.py -r 1fpu_receptor.pdbqt -s A:THR315

This will create two different files, one containing only the rigid part of the protein, and the other one containing Thr315 as flexible residue:

.. code-block:: console

	1fpu_rigid.pdbqt           # rigid part
	1fpu_flex.pdbqt            # flexible sidechain of Thr315


2. Prepare ligand
-----------------

For the molecular docking with flexible sidechains, we will use the ligand file `1iep_ligand.pdbqt` from the previous tutorial :ref:`basic_docking`.

3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

As well as for the docking with a fully rigid receptor, we need to generate a grid parameter file (.gpf) to precalculate the affinity map for each atom types. However, instead of using the full receptor, affinity maps will be calculated only for the rigid part of the receptor (`1fpu_rigid.pdbqt`).

To prepare the gpf file for the rigid part of the receptor:

.. code-block:: bash

	$ prepare_gpf.py -l xxxx.pdbqt -r 1fpu_rigid.pdbqt

.. code-block:: console
	:caption: Content of the grid parameter file (**1fpu_flex.gpf**) for the receptor c-Abl (**1fpu_rigid.pdbqt**)

	npts 40 46 52                        # num.grid points in xyz
	gridfld 1fpu.maps.fld                # grid_data_file
	spacing 0.375                        # spacing(A)
	receptor_types A C HD N OA SA        # receptor atom types
	ligand_types A C NA OA N HD          # ligand atom types
	receptor 1fpu_rigid.pdbqt            # macromolecule
	gridcenter 15.19 53.903 14.661       # xyz-coordinates or auto
	smooth 0.5                           # store minimum energy w/in rad(A)
	map 1fpu.A.map                       # atom-specific affinity map
	map 1fpu.C.map                       # atom-specific affinity map
	map 1fpu.NA.map                      # atom-specific affinity map
	map 1fpu.OA.map                      # atom-specific affinity map
	map 1fpu.N.map                       # atom-specific affinity map
	map 1fpu.HD.map                      # atom-specific affinity map
	elecmap 1fpu.e.map                   # electrostatic potential map
	dsolvmap 1fpu.d.map                  # desolvation potential map
	dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

To execute `autogrid4` using `1fpu_flex.gpf`, run the folllowing command line:

.. code-block:: bash

	$ autogrid4 -p 1fpu_flex.gpf -l 1fpu_flex.glg

You should obtain as well the following files:

.. code-block:: console

	1fpu.maps.fld       # grid data file
	1fpu.*.map          # affinity maps for different atom types
	1fpu.d.map          # desolvation map
	1fpu.e.map          # electrostatic map


4. Running AutoDock Vina
------------------------

4.a. Using AutoDock forcefield
______________________________

While using the AutoDock4 forcefield, only the flex part of the receptor is necessary, as well as the affinity maps. Once the receptor (flex part), ligand and maps were prepared, you can perform the flexible side-chain docking by simply running the following command line:

.. code-block:: bash

	$ vina --flex 1fpu_flex.pdbqt --ligand 1iep_ligand.pdbqt --maps 1fpu --scoring ad4 \
	       --exhaustiveness 24 --out 1fpu_ligand_flex.pdbqt

4.b. Using Vina forcefield
__________________________

As well as for the fully rigid molecular docking, you only need to specify the center and dimensions (in Angstrom) of the grid.

.. code-block:: console
	:caption: Content of the config file (**box.txt**) for AutoDock Vina

	center_x = 15.19
	center_y = 53.903
	center_z = 14.661
	size_x = 15.0
	size_y = 17.25
	size_z = 19.5

However, when using the Vina forcefield, you will need to specify both the rigid (needed to compute internally the affinity maps) and flex part of receptor. To perform the same docking experiment but using Vina forcefield run the following command line:

.. code-block:: bash

	$ vina --receptor 1fpu_rigid.pdbqt --flex 1fpu_flex.pdbqt --ligand 1iep_ligand.pdbqt \
	       --config box.txt --exhaustiveness 24 --out 1fpu_ligand_flex.pdbqt

.. tip::

	Alternatively, you can use the Vinardo forcefield by adding the `--scoring vinardo` option.

5. Results
----------

Analyze the flexible docking results using ADT (Fig. 3b), as described in Step 5A(iv).

