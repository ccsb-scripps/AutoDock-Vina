.. _macrocycle_docking:

Docking with macrocycles
========================

AutoDock Vina is not able to manage directly the flexibility associated with bonds in cyclic molecules. Different approaches can be used to dock macrocyclic molecules, like identifying one or more low energy conformations and docking them as different ligands, but generating them and docking them separately can be a time-consuming task. But AutoDock Vina (and AD4) has a specialized protocol to dock macrocycles while modeling their flexibility on-the-fly. 

The current implementation is described in our paper on the D3R Grand Challenge 4 (see below) and is summarized herein:

    1. One of the bonds in the ring structure is broken, resulting in an open form of the macrocycle that removes the need for correlated torsional variations, and enabling torsional degrees of freedom to be explored independently. 
    2. To each of the atoms of the broken bond, a dummy-atom is added. The distance between a dummy-atom and its parent atom is the same as the length of the broken bond, and the 1-3 angle also matches the original geometry. 
    3. During the docking, a linear attractive potential is applied on each dummy-atom to restore the bond resulting in the closed ring form. Thus, macrocycle conformations are sampled while adapting to the binding pocket, at the cost of increased search complexity with the added extra rotatable bonds. 

.. note::
    This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the meeko Python package.

.. note::
    If you are using this tutorial for your works, you can cite the following paper:

    - Forli, S., & Botta, M. (2007). Lennard-Jones potential and dummy atom settings to overcome the AUTODOCK limitation in treating flexible ring systems. Journal of chemical information and modeling, 47(4), 1481-1492.
    - Santos-Martins, D., Eberhardt, J., Bianco, G., Solis-Vasquez, L., Ambrosio, F. A., Koch, A., & Forli, S. (2019). D3R Grand Challenge 4: prospective pose prediction of BACE1 ligands with AutoDock-GPU. Journal of Computer-Aided Molecular Design, 33(12), 1071-1081.

For this tutorial, we are going to use the dataset from the Drug Design Data Resource (D3R) Grand Challenge 4 (GC4) organized in 2018. For more information about the D3R organization and it role in the computer-aided drug discovery community: `https://drugdesigndata.org <https://drugdesigndata.org/>`_. In second stage (1b) of this challenge, participants were asked to predict the bound poses of 20 BACE1 ligands knowing the crystallographic structure of the protein. A macrocycle is present in 19 out of 20 ligands, making pose predictions particularly challenging for the reasons described earlier. The docking procedure being the same for all 20 ligands, we are going to redock only one of the macrocycle in this tutorial.

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/docking_with_macrocycles/data`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/docking_with_macrocycles>`_). If you ever feel lost, you can always take a look at the solution here: ``AutoDock-Vina/example/docking_with_macrocycles/solution``. All the Python scripts used here (except for ``prepare_receptor`` and ``mk_prepare_ligand.py``) are located in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.

1. Preparing the receptor
-------------------------

The receptor can be prepared using the method described earlier in the following tutorials: :ref:`basic_docking`. The file ``BACE_1_receptorH.pdb`` is provided (see ``data`` directory located at ``<autodock-vina_directory>/example/macrocycle_docking/``). This file contains the receptor coordinates of chain A.

.. code-block:: bash

    $ prepare_receptor -r BACE_1_receptorH.pdb -o BACE_1_receptor.pdbqt

The output PDBQT file ``BACE_1_receptor.pdbqt`` is available in ``solution`` directory if necessary.

2. Preparing the ligand
-----------------------

For this particular tutorial we are going to use the python package ``Meeko`` for preparing the macrocycle ligand. If you didn't install yet ``Meeko`` you can check the instruction here: :ref:`docking_requirements`. As well, the second step consists to prepare the ligand, by converting the MOL2 file ``BACE_1_ligand.mol2`` to a PDBQT file readable by AutoDock Vina. The MOL2 file is located in the ``data`` directory.

.. code-block:: bash

    $ mk_prepare_ligand.py -i BACE_1_ligand.mol2 -m -o BACE_1_ligand.pdbqt

Different options are available for ``mk_prepare_ligand.py``, type  ``mk_prepare_ligand.py --help`` for more details. If you are not sure about this step, the output PDBQT file ``BACE_1_ligand.pdbqt`` is available in ``solution`` directory.

3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

Now, we have to define the grid space for the docking, typically, a 3D box around a the potential binding site of a receptor. During this step, we will create the input file for AutoGrid4, which will create an affinity map file for each atom types. The grid parameter file specifies an AutoGrid4 calculation, including the size and location of the grid, the atom types that will be used, the coordinate file for the rigid receptor, and other parameters for calculation of the grids.

To prepare the gpf file for AutoGrid4, your can use the ``prepare_gpf.py`` command line tool. This Python script is available here: ``<autodock-vina_directory>/example/autodock_scripts``. However, no extra maps are needed for the ``G0`` and ``CG0`` atoms because, for sake of evaluation of ligand-protein interaction, they are considered as normal carbon atoms. Therefore, C maps are used in their place and so we are going to ignore those atoms by manually specifying the ligand atom types ``-p ligand_types='A,C,OA,N,HD'``.

.. code-block:: bash

    $ pythonsh ../../autodock_scripts/prepare_gpf.py -l BACE_1_ligand.pdbqt \
               -r BACE_1_receptor.pdbqt -y -p ligand_types='A,C,OA,N,HD' \
               -p npts='54,54,54'

The option ``-y`` specifies that we want to center automatically the grid around the ligand. For more information about ``prepare_gpf.py``, type ``pythonsh prepare_gpf.py -h``. At the end you should obtain the following GPF file ``1iep_receptor.gpf`` containing those lines:

.. code-block:: console
    :caption: Content of the grid parameter file (**BACE_1_receptor.gpf**) for the receptor BACE (**BACE_1_receptor.pdbqt**)

    npts 54 54 54                        # num.grid points in xyz
    gridfld BACE_1_receptor.maps.fld     # grid_data_file
    spacing 0.375                        # spacing(A)
    receptor_types A C NA OA N SA HD     # receptor atom types
    ligand_types A C OA N HD             # ligand atom types
    receptor BACE_1_receptor.pdbqt       # macromolecule
    gridcenter 30.103 6.152 15.584       # xyz-coordinates or auto
    smooth 0.5                           # store minimum energy w/in rad(A)
    map BACE_1_receptor.A.map            # atom-specific affinity map
    map BACE_1_receptor.C.map            # atom-specific affinity map
    map BACE_1_receptor.OA.map           # atom-specific affinity map
    map BACE_1_receptor.N.map            # atom-specific affinity map
    map BACE_1_receptor.HD.map           # atom-specific affinity map
    elecmap BACE_1_receptor.e.map        # electrostatic potential map
    dsolvmap BACE_1_receptor.d.map              # desolvation potential map
    dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

After creating the GPF file, and now we can use the ``autogrid4`` command to generate the different map files that will be used for the molecular docking:

.. code-block:: bash

    $ autogrid4 -p BACE_1_receptor.gpf -l BACE_1_receptor.glg

From this command you should have generated the following files:

.. code-block:: console

    1iep_receptor.maps.fld       # grid data file
    1iep_receptor.*.map          # affinity maps for A, C, HD, N, OA atom types
    1iep_receptor.d.map          # desolvation map
    1iep_receptor.e.map          # electrostatic map

4. Running AutoDock Vina
------------------------

4.a. Using AutoDock4 forcefield
_______________________________

When using the AutoDock4 forcefield, you only need to provide the affinity maps and the ligand, while specifying that the forcefield used will be AutoDock4 using the option ``--scoring ad4``.

.. code-block:: bash

    $ vina --ligand BACE_1_ligand.pdbqt --maps BACE_1_receptor --scoring ad4 \
           --exhaustiveness 32 --out BACE_1_ligand_ad4_out.pdbqt

Running AutoDock Vina will write a PDBQT file called ``BACE_1_ligand_ad4_out.pdbqt``.

4.b. Using Vina forcefield
__________________________

As well as for the fully rigid molecular docking, you only need to specify the center and dimensions (in Angstrom) of the grid. Here, instead of specifying each parameters for the grid box using the arguments ``--center_x, --center_y, --center_z`` and ``--size_x, --size_y, --size_z``, we will also store all those informations in a text file ``BACE_1_receptor_vina_box.txt``.

.. code-block:: console
    :caption: Content of the config file (**BACE_1_receptor_vina_box.txt**) for AutoDock Vina

    center_x = 30.103
    center_y = 6.152
    center_z = 15.584
    size_x = 20
    size_y = 20
    size_z = 20

However, when using the Vina forcefield, you will need to specify the receptor ``BACE_1_receptor.pdbqt`` (needed to compute internally the affinity maps). To perform the same docking experiment but using Vina forcefield run the following command line:

.. code-block:: bash

    $ vina --receptor BACE_1_receptor.pdbqt --ligand BACE_1_ligand.pdbqt \
           --config BACE_1_receptor_vina_box.txt \
           --exhaustiveness 32 --out BACE_1_ligand_vina_out.pdbqt

.. tip::

    Alternatively, you can use the Vinardo forcefield by adding the ``--scoring vinardo`` option.

Running AutoDock Vina will write a PDBQT file called ``BACE_1_ligand_vina_out.pdbqt``.

5. Results
----------

.. warning::
    
    Please don't forget that energy scores giving by the AutoDock and Vina forcefield are not comparable between each other.

5.a. Using AutoDock forcefield
______________________________

The predicted free energy of binding should be about ``-13 kcal/mol`` for the best pose and should corresponds to the crystallographic pose ``BACE_1_ligand.mol2``. 

.. code-block:: console

    Scoring function : ad4
    Ligand: BACE_1_ligand.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: -226896966) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -12.98          0          0
       2       -12.18      1.046      1.292
       3       -12.16      1.719      3.126
       4       -10.92      2.008      3.019
       5       -10.35      3.209      5.331
       6       -9.841      3.061      4.709
       7       -9.476       3.08      8.727
       8       -9.094      3.794       5.73
       9       -8.438       3.43      8.918

5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-11 kcal/mol``. Using the Vina scoring function, the best pose shows also an excellent overlap with the crystallographic coordinates.

.. code-block:: console

    Scoring function : vina
    Rigid receptor: BACE_1_receptor.pdbqt
    Ligand: BACE_1_ligand.pdbqt
    Center: X 30.103 Y 6.152 Z 15.584
    Size: X 20 Y 20 Z 20
    Grid space: 0.375
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Computing Vina grid ... done.
    Performing docking (random seed: -1673486704) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -11.17          0          0
       2        -9.66       3.78      7.156
       3       -9.638      3.147      5.394
       4       -9.563      1.416      2.844
       5       -9.442      5.028      8.249
       6       -9.374      2.683      8.933
       7       -9.342      2.713      9.145
       8       -9.226      3.415      5.818
       9       -9.107      4.728      8.382
