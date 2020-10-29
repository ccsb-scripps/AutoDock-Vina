.. _basic_docking:

Basic docking
=============

Let's start with our first example of docking, where the typical usage pattern would be to dock a single molecule into a rigid receptor. In this example we will dock the approved anticancer drug `imatinib <https://en.wikipedia.org/wiki/Imatinib>`_ (Gleevec; PDB entry `1iep <https://www.rcsb.org/structure/1IEP>`_) in the structure of c-Abl using AutoDock Vina. The target for this protocol is the kinase domain of the proto-oncogene tyrosine protein kinase c-Abl. The protein is an important target for cancer chemotherapy—in particular, the treatment of chronic myelogenous leukemia.

.. note::
	This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the raccoon Python package.

.. note::
	The materials present is this tutorial can be also found here: `https://www.nature.com/articles/nprot.2016.051 <https://www.nature.com/articles/nprot.2016.051>`_. If you are using this tutorial for your works, you can cite the following paper:

	- Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

1. Preparing the receptor
-------------------------

During this step, we will create a PDBQT file of our receptor containing only the polar hydrogen atoms as well as partial charges. For this step, we will use the ``prepare_receptor`` command tool from the ADFR Suite. As a prerequisite, a receptor coordinate file must contain all hydrogen atoms. If hydrogen atoms are absent in the protein structure file, you can add the ``-A "hydrogens"`` flag. Many tools exist to add missing hydrogen atoms to a protein, one popular choice would be to use `REDUCE <http://kinemage.biochem.duke.edu/software/reduce.php>`_. If you are using experimental structures (for instance, from the PDB), use a text editor to remove waters, ligands, cofactors, ions deemed unnecessary for the docking. The file ``1iep_receptorH.pdb`` is provided (see ``<autodock-vina_directory>/example/basic_docking/data`` directory). This file contains the receptor coordinates taken from PDB entry ``1iep``.

.. code-block:: bash

	$ prepare_receptor -r 1iep_receptorH.pdb -o 1iep_receptor.pdbqt

Other options are available for ``prepare_receptor`` by typing ``prepare_receptor -h``. If you are not sure about this step, the output PDBQT file ``1iep_receptor.pdbqt`` is available in ``<autodock-vina_directory>/example/basic_docking/solution`` directory.


2. Preparing the ligand
-----------------------

This step is very similar to the previous step. We will also create a PDBQT file from a ligand molecule file (in MOL2 or PDB format). As well as for the receptor, the coordinate set must includes all hydrogen atoms according to the choosen protonation state. This may be obtained in a variety of ways, including with experimental coordinates from the `PDB <https://www.rcsb.org>`_ or `Cambridge Crystallographic Database <http://www.ccdc.cam.ac.uk>`_. The file ``1iep_ligandH.pdb`` is also provided (see ``<autodock-vina_directory>/example/basic_docking/data`` directory). This file includes ligand coordinates taken from PDB entry ``1iep``, to which all hydrogen atoms have been added and manually adjusted to the known protonation state.

.. code-block:: bash

	$ prepare_ligand -r 1iep_ligandH.pdb -o 1iep_ligand.pdbqt

As well, different options are available for ``prepare_ligand``, type  ``prepare_ligand -h`` for more details. If you are not sure about this step, the output PDBQT file ``1iep_ligand.pdbqt`` is available in ``<autodock-vina_directory>/example/basic_docking/solution`` directory.


3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

Now, we have to define the grid space for the docking, typically, a 3D box around a the potential binding site of a receptor. During this step, we will create the input file for AutoGrid4, which will create an affinity map file for each atom types. The grid parameter file specifies an AutoGrid4 calculation, including the size and location of the grid, the atom types that will be used, the coordinate file for the rigid receptor, and other parameters for calculation of the grids.

To prepare the gpf file for AutoGrid4, your can use the ``prepare_gpf.py`` command line tool. This Python script is available here: ``<autodock-vina_directory>/example/basic_docking/scripts``.

.. code-block:: bash

	$ pythonsh <script_directory>/prepare_gpf.py -l 1iep_ligand.pdbqt -r 1iep_receptor.pdbqt -y

The option ``-y`` specifies that we want to center automatically the grid around the ligand. For more information about ``prepare_gpf.py``, type ``pythonsh prepare_gpf.py -h``. At the end you should obtain the following GPF file ``1iep_receptor.gpf`` containing those lines:


.. code-block:: console
	:caption: Content of the grid parameter file (**1iep_receptor.gpf**) for the receptor c-Abl (**1iep_receptor.pdbqt**)

	npts 40 45 40                        # num.grid points in xyz
	gridfld 1iep_receptor.maps.fld       # grid_data_file
	spacing 0.375                        # spacing(A)
	receptor_types A C OA N SA HD        # receptor atom types
	ligand_types A C NA OA N HD          # ligand atom types
	receptor 1iep_receptor.pdbqt         # macromolecule
	gridcenter 15.662 53.211 15.546      # xyz-coordinates or auto
	smooth 0.5                           # store minimum energy w/in rad(A)
	map 1iep_receptor.A.map              # atom-specific affinity map
	map 1iep_receptor.C.map              # atom-specific affinity map
	map 1iep_receptor.NA.map             # atom-specific affinity map
	map 1iep_receptor.OA.map             # atom-specific affinity map
	map 1iep_receptor.N.map              # atom-specific affinity map
	map 1iep_receptor.HD.map             # atom-specific affinity map
	elecmap 1iep_receptor.e.map          # electrostatic potential map
	dsolvmap 1iep_receptor.d.map              # desolvation potential map
	dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

After creating the GPF file, and now we can use the ``autogrid4`` command to generate the different map files that will be used for the molecular docking:

.. code-block:: bash

	$ autogrid4 -p 1iep.gpf -l 1iep.glg

From this command you should have generated the following files:

.. code-block:: console

	1iep_receptor.maps.fld       # grid data file
	1iep_receptor.*.map          # affinity maps for A, C, HD, H, NA, N, OA atom types
	1iep_receptor.d.map          # desolvation map
	1iep_receptor.e.map          # electrostatic map

4. Running AutoDock Vina
------------------------

The imatinib ligand used in this protocol is challenging, and Vina will occasionally not find the correct pose with the default parameters. Vina provides a parameter called ``exhaustiveness`` to change the amount of computational effort used during a docking experiment. The default exhaustiveness value is ``8``; increasing this to ``32`` will give a more consistent docking result. At this point of the tutorial, you have the choice to decide to run the molecular docking using either the ``AutoDock`` forcefield (requires affinity maps, see previous step) or using the ``Vina`` forcefield (no need for affinity maps).

4.a. Using AutoDock forcefield
______________________________

When using the AutoDock4 forcefield, you only need to provide the affinity maps and the ligand, while specifying that the forcefield used will be AutoDock4 using the option ``--scoring ad4``.

.. code-block:: bash

	$ vina  --ligand 1iep_ligand.pdbqt --maps 1iep_receptor --scoring ad4 \
	        --exhaustiveness 32 --out 1iep_ligand_ad4_out.pdbqt

Running AutoDock Vina will write a PDBQT file called ``1iep_ligand_ad4_out.pdbqt`` contaning all the poses found during the molecular docking and also present docking information to the terminal window.

4.b. Using Vina forcefield
__________________________

Contrary to AutoDock4, you don't need to precalculate the affinity grid maps with ``autogrid4`` when using the Vina forcefield. AutoDock Vina computes those maps internally before the docking. However, you still need to specify the center and dimensions (in Angstrom) of the grid space, as well as the receptor. Here, instead of specifying each parameters for the grid box using the arguments ``--center_x, --center_y, --center_z`` and ``--size_x, --size_y, --size_z``, we will store all those informations in a text file ``1iep_receptor_vina_box.txt``.

.. code-block:: console
	:caption: Content of the config file (**1iep_receptor_vina_box.txt**) for AutoDock Vina

	center_x = 15.662
	center_y = 53.211
	center_z = 15.546
	size_x = 15.0
	size_y = 16.875
	size_z = 15.0

.. code-block:: bash

	$ vina --receptor 1iep_receptor.pdbqt --ligand 1iep_ligand.pdbqt \
	       --config 1iep_receptor_vina_box.txt \
	       --exhaustiveness=32 --out 1iep_ligand_vina_out.pdbqt

.. tip::

	Alternatively, you can use the Vinardo forcefield by adding the ``--scoring vinardo`` option.

Running AutoDock Vina will write a PDBQT file called ``1iep_ligand_vina_out.pdbqt``.

5. Results
----------

With ``exhaustiveness`` set to ``32``, Vina will most often give a single docked pose with this energy. With the lower default exhaustiveness, several poses flipped end to end, with less favorable energy, may be reported.

.. warning::
	
	Please don't forget that energy scores giving by the AutoDock and Vina forcefield are not comparable between each other.

5.a. Using AutoDock forcefield
______________________________

The predicted free energy of binding should be about ``-14 kcal/mol`` for poses that are similar to the crystallographic pose.

.. code-block:: console

	Scoring function : ad4
	Ligand: 1iep_ligand.pdbqt
	Exhaustiveness: 32
	CPU: 0
	Verbosity: 1

	Reading AD4.2 maps ... done.

	0%   10   20   30   40   50   60   70   80   90   100%
	|----|----|----|----|----|----|----|----|----|----|
	***************************************************

	mode |   affinity | dist from best mode
	     | (kcal/mol) | rmsd l.b.| rmsd u.b.
	-----+------------+----------+----------
	   1       -14.27          0          0
	   2       -13.37      1.446      1.994
	   3       -12.64      1.433      2.053
	   4       -11.94       3.84      11.35
	   5       -11.14       3.22      11.07
	   6       -11.09      2.271      4.243
	   7       -10.85      1.928      12.36
	   8       -10.85      2.082      12.45
	   9       -10.03      2.151      11.26

5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-10 kcal/mol``.

.. code-block:: console

	Scoring function : vina
	Rigid receptor: 1iep_receptor.pdbqt
	Ligand: 1iep_ligand.pdbqt
	Center: X 15.662 Y 53.211 Z 15.546
	Size: X 15 Y 16.875 Z 15
	Grid space: 0.375
	Exhaustiveness: 32
	CPU: 0
	Verbosity: 1

	Computing Vina grid ... done.

	0%   10   20   30   40   50   60   70   80   90   100%
	|----|----|----|----|----|----|----|----|----|----|
	***************************************************

	mode |   affinity | dist from best mode
	     | (kcal/mol) | rmsd l.b.| rmsd u.b.
	-----+------------+----------+----------
	   1       -10.28          0          0
	   2       -9.616      2.281      12.51
	   3       -8.552      1.818      12.78
	   4       -8.218      1.501      2.332
	   5       -6.968      1.683      2.594
	   6       -5.675      1.825      3.873
	   7       -5.385      1.943      12.86
	   8       -5.342      2.116      11.57
	   9       -5.335        2.8      11.35
