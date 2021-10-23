.. _basic_docking:

Basic docking
=============

Let's start with our first example of docking, where the typical usage pattern would be to dock a single molecule into a rigid receptor. In this example we will dock the approved anticancer drug `imatinib <https://en.wikipedia.org/wiki/Imatinib>`_ (Gleevec; PDB entry `1iep <https://www.rcsb.org/structure/1IEP>`_) in the structure of c-Abl using AutoDock Vina. The target for this protocol is the kinase domain of the proto-oncogene tyrosine protein kinase c-Abl. The protein is an important target for cancer chemotherapy—in particular, the treatment of chronic myelogenous leukemia.

.. note::
    This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the meeko Python package.

.. note::
    The materials present is this tutorial can be also found here: `https://www.nature.com/articles/nprot.2016.051 <https://www.nature.com/articles/nprot.2016.051>`_. If you are using this tutorial for your works, you can cite the following paper:

    - Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/basic_docking/data`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/basic_docking>`_). If you ever feel lost, you can always take a look at the solution here: ``AutoDock-Vina/example/basic_docking/solution``. All the Python scripts used here (except for ``prepare_receptor`` and ``mk_prepare_ligand.py``) are located in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.


1. Preparing the receptor
-------------------------

During this step, we will create a PDBQT file of our receptor containing only the polar hydrogen atoms as well as partial charges. For this step, we will use the ``prepare_receptor`` command tool from the ADFR Suite. As a prerequisite, a receptor coordinate file must contain all hydrogen atoms. If hydrogen atoms are absent in the protein structure file, you can add the ``-A "hydrogens"`` flag. Many tools exist to add missing hydrogen atoms to a protein, one popular choice would be to use `REDUCE <http://kinemage.biochem.duke.edu/software/reduce.php>`_. If you are using experimental structures (for instance, from the `Protein Data Bank <https://www.rcsb.org>`_), use a text editor to remove waters, ligands, cofactors, ions deemed unnecessary for the docking. This file contains the receptor coordinates taken from PDB entry `1iep <https://www.rcsb.org/structure/1IEP>`_.

.. code-block:: bash

    $ prepare_receptor -r 1iep_receptorH.pdb -o 1iep_receptor.pdbqt

Other options are available for ``prepare_receptor`` by typing ``prepare_receptor -h``. If you are not sure about this step, the output PDBQT file ``1iep_receptor.pdbqt`` is available in ``solution`` directory.


2. Preparing the ligand
-----------------------

This step is very similar to the previous step. We will also create a PDBQT file from a ligand molecule file (in MOL/MOL2 or SDF format) using the ``Meeko`` python package (see installation instruction here: :ref:`docking_requirements`). For convenience, the file ``1iep_ligand.sdf`` is provided (see ``data`` directory). But you can obtain it directly from the `PDB <https://www.rcsb.org>`_ here: `1iep <https://www.rcsb.org/structure/1IEP>`_ (see ``Download instance Coordinates`` link for the STI molecule). Since the ligand file does not include the hydrogen atoms, we are going to automatically add them and correct the protonation for a pH of 7.4.

.. warning::
  
  We strongly advice you against using PDB format for preparing small molecules, since it does not contain information about bond connections. Please don't forget to always check the protonation state of your molecules before docking. Your success can sometimes hang by just an hydrogen atom. ;-)

.. code-block:: bash

    $ mk_prepare_ligand.py -i 1iep_ligand.sdf -o 1iep_ligand.pdbqt --pH 7.4

Other options are available for ``mk_prepare_ligand.py`` by typing ``mk_prepare_ligand.py --help``. If you are not sure about this step, the output PDBQT file ``1iep_ligand.pdbqt`` is available in ``solution`` directory.


3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

Now, we have to define the grid space for the docking, typically, a 3D box around a the potential binding site of a receptor. During this step, we will create the input file for AutoGrid4, which will create an affinity map file for each atom types. The grid parameter file specifies an AutoGrid4 calculation, including the size and location of the grid, the atom types that will be used, the coordinate file for the rigid receptor, and other parameters for calculation of the grids.

To prepare the gpf file for AutoGrid4, you can use the ``prepare_gpf.py`` command line tool.

.. code-block:: bash

    $ pythonsh <script_directory>/prepare_gpf.py -l 1iep_ligand.pdbqt -r 1iep_receptor.pdbqt -y

The option ``-y`` specifies that we want to center automatically the grid around the ligand. For more information about ``prepare_gpf.py``, type ``pythonsh prepare_gpf.py -h``. At the end you should obtain the following GPF file ``1iep_receptor.gpf`` containing those lines:


.. code-block:: console
    :caption: Content of the grid parameter file (**1iep_receptor.gpf**) for the receptor c-Abl (**1iep_receptor.pdbqt**)

    npts 54 54 54                        # num.grid points in xyz
    gridfld 1iep_receptor.maps.fld       # grid_data_file
    spacing 0.375                        # spacing(A)
    receptor_types A C OA N SA HD        # receptor atom types
    ligand_types A C NA OA N HD          # ligand atom types
    receptor 1iep_receptor.pdbqt         # macromolecule
    gridcenter 15.190 53.903 16.917      # xyz-coordinates or auto
    smooth 0.5                           # store minimum energy w/in rad(A)
    map 1iep_receptor.A.map              # atom-specific affinity map
    map 1iep_receptor.C.map              # atom-specific affinity map
    map 1iep_receptor.NA.map             # atom-specific affinity map
    map 1iep_receptor.OA.map             # atom-specific affinity map
    map 1iep_receptor.N.map              # atom-specific affinity map
    map 1iep_receptor.HD.map             # atom-specific affinity map
    elecmap 1iep_receptor.e.map          # electrostatic potential map
    dsolvmap 1iep_receptor.d.map         # desolvation potential map
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

4.a. Using AutoDock4 forcefield
_______________________________

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

    center_x = 15.190
    center_y = 53.903
    center_z = 16.917
    size_x = 20.0
    size_y = 20.0
    size_z = 20.0

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
  Performing docking (random seed: -556654859) ... 
  0%   10   20   30   40   50   60   70   80   90   100%
  |----|----|----|----|----|----|----|----|----|----|
  ***************************************************

  mode |   affinity | dist from best mode
       | (kcal/mol) | rmsd l.b.| rmsd u.b.
  -----+------------+----------+----------
     1       -14.62          0          0
     2       -13.13      1.051      1.529
     3       -12.26      1.442      2.158
     4       -11.91      3.646       11.5
     5       -11.89      3.859      11.99
     6       -11.47      1.978      13.56
     7       -11.33      1.727      2.585
     8       -10.85      3.619      5.759
     9       -10.23      7.057       12.7

5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-13 kcal/mol``.

.. code-block:: console

  Scoring function : vina
  Rigid receptor: 1iep_receptor.pdbqt
  Ligand: 1iep_ligand.pdbqt
  Center: X 15.19 Y 53.903 Z 16.917
  Size: X 20 Y 20 Z 20
  Grid space: 0.375
  Exhaustiveness: 32
  CPU: 0
  Verbosity: 1

  Computing Vina grid ... done.
  Performing docking (random seed: -131415392) ... 
  0%   10   20   30   40   50   60   70   80   90   100%
  |----|----|----|----|----|----|----|----|----|----|
  ***************************************************

  mode |   affinity | dist from best mode
       | (kcal/mol) | rmsd l.b.| rmsd u.b.
  -----+------------+----------+----------
     1       -12.92          0          0
     2       -10.97      3.012      12.42
     3       -10.79      3.713      12.19
     4       -10.69      3.913      12.36
     5       -10.32      2.538      12.64
     6       -9.464      2.916      12.53
     7       -9.204       1.35      2.025
     8       -9.137      1.596      2.674
     9       -8.637      3.969      12.69
