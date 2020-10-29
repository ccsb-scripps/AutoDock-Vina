.. _hydrated_docking:

Hydrated docking
================

Introduction
------------

In physiological environments, proteins and other biological structures are surrounded by water molecules. When a small-molecule binds to a protein, it must displace most of the waters occupying the binding cavity. However, rarely are all water molecules displaced. Some waters can be so strongly bound and conserved among similar proteins that from a ligand-docking perspective they are considered a part of the target structure, altering the binding site topography. 

Thus a method was developed that uses the existing version of AutoDock but modifies the force field to model explicit bridging water molecules. The ligand is decorated with an ensemble of water molecules (by adding dummy atoms), which may or may not then contribute to the interactions. A modified AutoGrid map is then used during docking, giving a favorable score when the water is well placed and omitting the water if it overlaps with the receptor. A final script analyzes the docked results, retaining only those waters in appropriate positions. In tests, this method has shown improvement in the prediction of bound conformations of small fragment molecules, such as those used in fragment-based drug discovery.

In this tutorial, we are going to dock a fragment-size ligand (nicotine) with explicit water molecules in the acetylcholine binding protein (AChBP) structure (PDB entry `1uw6 <https://www.rcsb.org/structure/1uw6>`_). It was shown that the absence of water molecules could have a dramatic influence in docking performance leading to either inaccurate scoring and/or incorrect pose. With hydrated docking, fragment-sized ligands show a overall RMSD improvement.

.. note::

	This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the raccoon Python package.

.. note::
	
	The materials present is this tutorial can be also found here: `https://www.nature.com/articles/nprot.2016.051 <https://www.nature.com/articles/nprot.2016.051>`_. If you are using this tutorial or this docking method for your work, you can cite the following papers:

	- Forli, S., & Olson, A. J. (2012). A force field with discrete displaceable waters and desolvation entropy for hydrated ligand docking. Journal of medicinal chemistry, 55(2), 623-638.
	- Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

1. Preparing the receptor
-------------------------

The receptor can be prepared using the method described earlier in the following tutorials: :ref:`basic_docking` or :ref:`flexible_docking` if one wants to incorporate some sidechain flexibility. The file ``1uw6_receptorH.pdb`` is provided (see ``<autodock-vina_directory>/example/hydrated_docking/data`` directory). This file contains the receptor coordinates of chain A and B taken from the PDB entry ``1uw6``.

.. code-block:: bash

	$ prepare_receptor -r 1uw6_receptorH.pdb -o 1uw6_receptor.pdbqt

The output PDBQT file ``1uw6_receptor.pdbqt`` is available in ``<autodock-vina_directory>/example/flexible_docking/solution`` directory if necessary.

2. Preparing the ligand
-----------------------

For the hydrated docking, explicit water molecules (W atoms) must be added to a PDBQT file. And for that, we will use the ``wet.py`` command tool. For this step, the file ``1uw6_ligandH.pdb`` is also provided (see ``<autodock-vina_directory>/example/hydrated_docking/data`` directory). This file includes ligand coordinates taken from PDB entry ``1uw6`` with all hydrogen atoms already present.

We first need to generate the PDBQT file, and after add water positions to the ligand. For that, type and execute the following command lines:

.. code-block:: bash
	
	$ prepare_ligand -l 1uw6_ligandH.pdb -o 1uw6_ligand.pdbqt
	$ wet.py -i 1uw6_ligand.pdbqt -o 1uw6_ligand_hydro.pdbqt

From the ``wet.py`` command, you should obtain the following output:

.. code-block:: console

	PDBQT 1uw6_ligand.pdbqt
	NAME is  1uw6_ligand

	====================================
	  ________         __   
	 |  |  |  |.-----.|  |_ 
	 |  |  |  ||  -__||   _|
	 |________||_____||____|

	   processing 1uw6_ligand.pdbqt
	 - saving the file => 1uw6_ligand_hydro.pdbqt
	 - ignoring phosphate/sulphate groups
	 - hydratable atoms : 2 / 13 
	 - 2 waters added

In total, 2 water molecules were added to the fragment. For more information about the ``wet.py`` command and all the available options, just type ``wet.py``. If you were not able to generate the ``1uw6_ligand.pdbqt`` or ``1uw6_ligand_hydro.pdbqt`` files, you can look at the ``solution`` directory.

3. Generating affinity maps
---------------------------

As well as for the :ref:`basic_docking` or :ref:`flexible_docking` tutorials, we will also need to calculate the affinity maps for each atom types present in the ligand. However, this time we will also need to craft a special ``W`` affinity map for the water molecules attached to the ligand. This maps will be a combinaison of both the ``OA`` and ``HD`` affinity maps.

The first step will be to create the GPF file and then calculate the grid maps following the standard AutoDock protocol, checking that OA and HD types are present in the ligand atom set. If not, the GPF file must be modified to include them; i.e. :
	
	- add OA and HD to the line “ligand_types …. ”
	- add lines “map protein.HD.map” and “map protein.OA.map”

.. warning::
	
	For this step, don't use ``1uw6_ligand_hydro.pdbqt`` for creating the GPF file.

.. code-block:: bash

	$ pythonsh <script_directory>/prepare_gpf.py -l 1uw6_ligand.pdbqt -r 1uw6_receptor.pdbqt -y

The option ``-y`` specifies that we want to center automatically the grid around the ligand. After manually adding the ``OA`` atom type, you should have a GPF file called ``1uw6_receptor.gpf`` that looks like this:

.. code-block:: console
	:caption: Content of the grid parameter file (**1uw6_receptor.gpf**) for the receptor (**1uw6_receptor.pdbqt**)

	npts 40 40 40                        # num.grid points in xyz
	gridfld 1uw6_receptor.maps.fld       # grid_data_file
	spacing 0.375                        # spacing(A)
	receptor_types A C NA OA N SA HD     # receptor atom types
	ligand_types A NA HD N OA            # ligand atom types * ADD OA/HD TYPES IF NOT PRESENT *
	receptor 1uw6_receptor.pdbqt         # macromolecule
	gridcenter 83.640 69.684 -10.124     # xyz-coordinates or auto
	smooth 0.5                           # store minimum energy w/in rad(A)
	map 1uw6_receptor.A.map              # atom-specific affinity map
	map 1uw6_receptor.NA.map             # atom-specific affinity map
	map 1uw6_receptor.HD.map             # atom-specific affinity map * ADD HD IF NOT PRESENT *
	map 1uw6_receptor.N.map              # atom-specific affinity map
	map 1uw6_receptor.OA.map             # atom-specific affinity map * ADD OA IF NOT PRESENT *
	elecmap 1uw6_receptor.e.map          # electrostatic potential map
	dsolvmap 1uw6_receptor.d.map              # desolvation potential map
	dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

You can now execute ``autogrid4`` using the GPF file called ``1uw6_receptor.gpf`` and generate the additional water map ``W`` by combining ``OA`` and ``HD`` affinity maps using ``mapwater.py``:

.. code-block:: bash

	$ autogrid4 -p 1uw6_receptor.gpf -l 1uw6_receptor.glg
	$ mapwater.py -r 1uw6_receptor.pdbqt -s 1uw6_receptor.W.map

For more informations about the ``mapwater.py`` command tool and all the available options, just type ``mapwater.py``. After executing this command, you should obtain a new affinity map called ``1uw6_receptor.W.map`` and the following the output:

.. code-block:: console

	ADD PWD AND FILE SUMMARY
	  receptor :  1uw6_receptor.pdbqt
	      OA map -> 1uw6_receptor.OA.map
	      HD map -> 1uw6_receptor.HD.map
	 => Water map weight : DEFAULT [ 0.60 ]

	  MapWater generator
	 =====================
	  mode      :  BEST
	  weight    :   0.6
	  HD_weight :   1.0
	  OA_weight :   1.0
	  entropy   :   -0.2

	     Output info  
	  --------------------
	  filename  : 1uw6_receptor.W.map
	  OA points : 91.73%
	  HD points : 8.27%

	  lowest  map value : -0.99
	  highest map value : -0.01

4. Running AutoDock Vina
------------------------

4.a. Using AutoDock forcefield
______________________________

Now that you generated the ligand with explicit water molecules attached (``1uw6_ligand_hydro.pdbqt``) and the extra affinity map for the  ``W`` atom type (``1uw6_receptor.W.map``), you can do the molecular docking with AutoDock Vina using the AutoDock forcefield:

.. code-block:: bash

	$ vina  --ligand 1uw6_ligand_hydro.pdbqt --maps 1uw6_receptor --scoring ad4 \
	        --exhaustiveness 32 --out 1uw6_ligand_hydro_ad4_out.pdbqt

4.b. Using Vina forcefield
__________________________

.. warning::
	
	While this method was calibrated and validated with the AutoDock4 forcefield, we strongly advice you against using this protocol with the Vina and Vinardo forcefield.

5. Results and post-processing
------------------------------

.. warning::

	Be aware that with this implementation of the method, it is difficult to compare results obtained with very diverse ligands without doing extra of post-processing on the results, because the energy estimation needs to be normalized. For this reason, the method is not suitable for virtual screenings. This doesn’t affect the structural accuracy, so comparisons within docking poses are fine. An improved scoring function to overcome this issue is in the works.

The predicted free energy of binding should be about ``-8.2 kcal/mol`` for poses that are similar to the crystallographic pose.

.. code-block:: console

	Scoring function : ad4
	Ligand: 1uw6_ligand_hydro.pdbqt
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
	   1       -8.256          0          0
	   2       -8.202      1.418      1.752
	   3       -7.821      1.929      2.596
	   4       -7.589      2.135      2.494
	   5       -7.391      2.079       2.57
	   6       -7.349      2.362      4.763
	   7       -7.346      2.065      5.306
	   8       -7.084      3.295      5.307
	   9       -6.994      3.102      5.487

Docking results are filtered by using the receptor to remove displaced waters and the W map file to rank the conserved ones as strong or weak water molecules.

.. code-block:: bash

	$ dry.py -r 1uw6_receptor.pdbqt -m 1uw6_receptor.W.map -i 1uw6_ligand_hydro_ad4_out.pdbqt

For more informations about the ``dry.py`` command tool and all the available options, just type ``dry.py``. Running the previous command should give you this output:

.. code-block:: console

	                  ____                      
	                 /\  _`\                    
	                 \ \ \/\ \  _ __  __  __    
	                  \ \ \ \ \/\`'__\\ \/\ \   
	                   \ \ \_\ \ \ \/\ \ \_\ \  
	                    \ \____/\ \_\ \/`____ \ 
	                     \/___/  \/_/  `/___/> \
	                                      /\___/
	                                      \/__/ 

	    
	========================== INPUT DATA =========================
	 importing ATOMS from  1uw6_ligand_hydro_ad4_out.pdbqt

	 [ using map file 1uw6_receptor.W.map ]
	===============================================================


	 receptor structure loaded	 		 [ 4069 atoms ]
	 receptor 5A shell extracted  			 [ 480 atoms in 5 A shell ] 
	 removing ligand/ligand overlapping waters	  [ 5 water(s) removed ]
	 removing ligand/receptor overlapping waters	  [ 8 water(s) removed ]

	 scanning grid map for conserved waters...	  [ filtered pose contains 5 waters ]

	 water grid score results [ map: 1uw6_receptor.W.map ] 
		 [ Water STRONG ( -0.92 ) +++ ]
		 [ Water STRONG ( -0.66 ) +++ ]
		 [ Water  WEAK  ( -0.50 )  +  ]
		 [ Water STRONG ( -0.83 ) +++ ]
		 [ Water STRONG ( -0.99 ) +++ ]

Waters are ranked (STRONG, WEAK) and scored inside the output file ``1uw6_ligand_hydro_ad4_out_DRY_SCORED.pdbqt`` with the calculated energy.