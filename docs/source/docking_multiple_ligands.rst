.. _multiple_ligands_docking:

Multiple ligands docking
========================

Vina is now able to dock simultaneously multiple ligands. This functionality may find application in fragment based drug design, where small molecules that bind the same target can be grown or combined into larger compounds with potentially better affinity.

The protein PDE in complex with two inhibitors (pdb id: `5x72 <https://www.rcsb.org/structure/5X72>`_) was used as an example to demonstrate the ability of the AutoDock Vina to dock successfully multiple ligands. The two inhibitors in this structure are stereoisomers, and only the R-isomer is able to bind in a specific region of the pocket, while both the R- and S-isomers can bind to the second location. 

.. note::
    This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the meeko Python package.

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/mulitple_ligands_docking/data`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/mulitple_ligands_docking>`_). If you ever feel lost, you can always take a look at the solution here: ``AutoDock-Vina/example/mulitple_ligands_docking/solution``. All the Python scripts used here (except for ``prepare_receptor`` and ``mk_prepare_ligand.py``) are located in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.

1. Preparing the flexible receptor
----------------------------------

Exactly like the :ref:`basic_docking` tutorial, the method requires a receptor coordinate file that includes all hydrogen atoms. The file ``5x72_receptorH.pdb`` is provided (see ``data`` directory located at ``<autodock-vina_directory>/example/multiple_ligands_docking/``). This file contains the receptor coordinates taken from the PDB entry ``5x72``. It was manually obtained by extracting the receptor coordinates (using an text editor) from the original PDB file ``5x72.pdb`` in the ``data`` directory, and the hydrogen atoms added using `reduce <http://kinemage.biochem.duke.edu/software/reduce.php>`_.

.. code-block:: bash
    
    $ prepare_receptor -r 5x72_receptorH.pdb -o 5x72_receptor.pdbqt

If you are not sure about this step, the output PDBQT file ``5x72_receptor.pdbqt`` is available in the ``solution`` directory.

2. Prepare ligands
------------------

Here, we will prepare two ligands instead of only one. We will start from the SDF files ``5x72_ligand_p59.sdf`` and ``5x72_ligand_p69.sdf`` located in the ``data`` directory. They were also obtained directly from the `PDB <https://www.rcsb.org>`_ here: `5x72 <https://www.rcsb.org/structure/5X72>`_ (see ``Download instance Coordinates`` link for the P59 and P69 molecules). Since the ligand files do not include the hydrogen atoms, we are going to automatically add them.

.. warning::
  
  We strongly advice you against using PDB format for preparing small molecules, since it does not contain information about bond connections. Please don't forget to always check the protonation state of your molecules before docking. Your success can sometimes hang by just an hydrogen atom. ;-)

.. code-block:: bash

    $ mk_prepare_ligand.py -i 5x72_ligand_p59.sdf -o 5x72_ligand_p59.pdbqt --add_hydrogen
    $ mk_prepare_ligand.py -i 5x72_ligand_p69.sdf -o 5x72_ligand_p69.pdbqt --add_hydrogen

The output PDBQT ``5x72_ligand_p59.pdbqt`` and ``5x72_ligand_p69.pdbqt`` can be found in the ``solution`` directory.

3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

As well as for the docking with a fully rigid receptor, we need to generate a GPF file to precalculate the affinity maps for each atom type. However, instead of using the full receptor, affinity maps will be calculated only for the rigid part of the receptor (``5x72_receptor.pdbqt``).

To prepare the GPF file for the rigid part of the receptor:

.. code-block:: bash

    $ pythonsh <script_directory>/prepare_gpf.py -l 5x72_ligand_p59.pdbqt -r 5x72_receptor.pdbqt \ 
               -p npts='80,64,64' -p gridcenter='-15 15 129' -o 5x72_receptor.gpf

This time we manually specified the center of the grid ``-p gridcenter='-15 15 129'`` as well as its size ``-p npts='80,64,64'``.

.. code-block:: console
    :caption: Content of the grid parameter file (**5x72_receptor.gpf**) for the receptor (**5x72_receptor.pdbqt**)

    npts 80 64 64                        # num.grid points in xyz
    gridfld 5x72_receptor.maps.fld       # grid_data_file
    spacing 0.375                        # spacing(A)
    receptor_types A C NA OA N SA HD     # receptor atom types
    ligand_types A C F OA N HD           # ligand atom types
    receptor 5x72_receptor.pdbqt         # macromolecule
    gridcenter -15.000 15.000 129.000    # xyz-coordinates or auto
    smooth 0.5                           # store minimum energy w/in rad(A)
    map 5x72_receptor.A.map              # atom-specific affinity map
    map 5x72_receptor.C.map              # atom-specific affinity map
    map 5x72_receptor.F.map              # atom-specific affinity map
    map 5x72_receptor.OA.map             # atom-specific affinity map
    map 5x72_receptor.N.map              # atom-specific affinity map
    map 5x72_receptor.HD.map             # atom-specific affinity map
    elecmap 5x72_receptor.e.map          # electrostatic potential map
    dsolvmap 5x72_receptor.d.map              # desolvation potential map
    dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

.. warning::

    You might have to manually edit the GPF file and add addtional atom types if the second ligand contains different atom types not present in the ligand used for creating the GPF file.

To execute ``autogrid4`` using ``5x72_receptor.gpf``, run the folllowing command line:

.. code-block:: bash

    $ autogrid4 -p 5x72_receptor.gpf -l 5x72_receptor_rigid.glg

You should obtain as well the following files:

.. code-block:: console

    1fpu_receptor.maps.fld       # grid data file
    1fpu_receptor.*.map          # affinity maps for A, C, HD, NA, N, OA atom types
    1fpu_receptor.d.map          # desolvation map
    1fpu_receptor.e.map          # electrostatic map

4. Running AutoDock Vina
------------------------

4.a. Using AutoDock4 forcefield
_______________________________

When using the AutoDock4 forcefield, you only need to provide the affinity maps and the ligand, while specifying that the forcefield used will be AutoDock4 using the option ``--scoring ad4``.

.. code-block:: bash

    $ vina --ligand 5x72_ligand_p59.pdbqt 5x72_ligand_p69.pdbqt --maps 5x72_receptor \ 
           --scoring ad4 --exhaustiveness 32 --out 5x72_ligand_ad4_out.pdbqt

4.b. Using Vina forcefield
__________________________

As well as for the fully rigid molecular docking, you only need to specify the center and dimensions (in Angstrom) of the grid. Here, instead of specifying each parameters for the grid box using the arguments ``--center_x, --center_y, --center_z`` and ``--size_x, --size_y, --size_z``, we will also store all those informations in a text file ``5x72_receptor_vina_box.txt``.

.. code-block:: console
    :caption: Content of the config file (**5x72_receptor_vina_box.txt**) for AutoDock Vina

    center_x = -15
    center_y = 15
    center_z = 129
    size_x = 30
    size_y = 24
    size_z = 24

However, when using the Vina forcefield, you will need to specify the receptor ``5x72_receptor.pdbqt`` (needed to compute internally the affinity maps). To perform the same docking experiment but using Vina forcefield run the following command line:

.. code-block:: bash

    $ vina --receptor 5x72_receptor.pdbqt --ligand 5x72_ligand_p59.pdbqt 5x72_ligand_p69.pdbqt \
           --config 5x72_receptor_vina_box.txt \
           --exhaustiveness 32 --out 5x72_ligand_vina_out.pdbqt

.. tip::

    Alternatively, you can use the Vinardo forcefield by adding the ``--scoring vinardo`` option.

Running AutoDock Vina will write a PDBQT file called ``5x72_ligand_flex_vina_out.pdbqt``.

5. Results
----------

.. warning::
    
    Please don't forget that energy scores giving by the AutoDock and Vina forcefield are not comparable between each other.

5.a. Using AutoDock forcefield
______________________________

The predicted free energy of binding should be about ``-18 kcal/mol`` for poses that are similar to the crystallographic pose. Using the AutoDock4 scoring function, the first two sets of poses (top 2) need to be considered to show also a good overlap with the crystallographic poses
 
.. code-block:: console

    Scoring function : ad4
    Ligands:
      - 5x72_ligand_p59.pdbqt
      - 5x72_ligand_p69.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: 1295744643) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -18.94          0          0
       2       -18.62      1.634      3.349
       3        -18.4      1.413      3.312
       4       -18.24      1.341      3.921
       5       -18.03      1.599      9.262
       6       -17.93      1.631      9.166
       7       -17.84      1.928      4.933
       8       -17.74       1.74      8.879
       9       -17.74          2      9.433

5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-21 kcal/mol``. Using the Vina scoring function, the best set of poses (top 1) shows an excellent overlap with the crystallographic coordinates for one of the isomers.

.. code-block:: console

    Scoring function : vina
    Rigid receptor: 5x72_receptor.pdbqt
    Ligands:
      - 5x72_ligand_p59.pdbqt
      - 5x72_ligand_p69.pdbqt
    Center: X -15 Y 15 Z 129
    Size: X 30 Y 24 Z 24
    Grid space: 0.375
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Computing Vina grid ... done.
    Performing docking (random seed: -2141167371) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -21.32          0          0
       2       -20.94      1.061      3.648
       3       -20.73      1.392      3.181
       4       -19.93      1.744      4.841
       5       -19.34      1.384      3.352
       6       -19.05      1.185      9.184
       7        -18.9      1.198      3.586
       8       -18.76      1.862      8.986
       9       -18.63      1.749      9.194
