.. _flexible_docking:

Flexible docking
================

The lack of receptor flexibility is arguably the greatest limitation in these types of docking methods. However, AutoDock Vina allows some limited flexibility of selected receptor side chains. In this tutorial, we will describe the cross-docking of the `imatinib molecule <https://en.wikipedia.org/wiki/Imatinib>`_ to c-Abl in PDB entry `1fpu <https://www.rcsb.org/structure/1FPU>`_, treating Thr315 as flexible. 

.. note::
    This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the meeko Python package.

.. note::
    The materials present is this tutorial can be also found here: `https://www.nature.com/articles/nprot.2016.051 <https://www.nature.com/articles/nprot.2016.051>`_. If you are using this tutorial for your works, you can cite the following paper:

    - Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/flexible_docking/data`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/flexible_docking>`_). If you ever feel lost, you can always take a look at the solution here: ``AutoDock-Vina/example/flexible_docking/solution``. All the Python scripts used here (except for ``prepare_receptor`` and ``mk_prepare_ligand.py``) are located in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.

1. Preparing the flexible receptor
----------------------------------

During this step, we are going to split the receptor coordinates into two PDBQT files: one for the rigid portion and one for the flexible side chains. As with the :ref:`basic_docking` tutorial, the method requires a receptor coordinate file that includes all hydrogen atoms. The file ``1fpu_receptorH.pdb`` is provided (see ``<autodock-vina_directory>/example/flexible_docking/data`` directory). It contains the receptor coordinates taken from PDB entry ``1fpu``.

.. code-block:: bash
    
    $ prepare_receptor -r 1fpu_receptorH.pdb -o 1fpu_receptor.pdbqt
    $ pythonsh <script_directory>/prepare_flexreceptor.py -r 1fpu_receptor.pdbqt -s THR315

Other options are available for ``prepare_flexreceptor.py`` with the ``-h`` option. This will create two different files, one containing only the rigid part of the protein, and the other one containing Thr315 as flexible residue:

.. code-block:: console

    1fpu_receptor_rigid.pdbqt           # rigid part
    1fpu_receptor_flex.pdbqt            # flexible sidechain of Thr315

If you are not sure about this step, the output PDBQT files ``1fpu_receptor_rigid.pdbqt`` and ``1fpu_receptor_flex.pdbqt`` are available in ``solution`` directory.


2. Prepare ligand
-----------------

For the molecular docking with flexible sidechains, we will use the ligand file ``1iep_ligand.pdbqt`` from the previous tutorial :ref:`basic_docking`.

3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

As well as for the docking with a fully rigid receptor, we need to generate a GPF file to precalculate the affinity map for each atom types. However, instead of using the full receptor, affinity maps will be calculated only for the rigid part of the receptor (``1fpu_receptor_rigid.pdbqt``).

To prepare the GPF file for the rigid part of the receptor:

.. code-block:: bash

    $ pythonsh <script_directory>/prepare_gpf.py -l 1iep_ligand.pdbqt -r 1fpu_receptor_rigid.pdbqt -y

Luckily for us, the structure ``1iep`` and ``1fpu`` are almost perfectly superposed already, so we can also center the grid around ``1iep_ligand.pdbqt``. Otherwise, the center of the grid can be specified using the option ``-p gridcenter='X,X,X'``.

.. code-block:: console
    :caption: Content of the grid parameter file (**1fpu_receptor_rigid.gpf**) for the receptor c-Abl (**1fpu_receptor_rigid.pdbqt**)

    npts 54 54 54                        # num.grid points in xyz
    gridfld 1fpu_receptor_rigid.maps.fld # grid_data_file
    spacing 0.375                        # spacing(A)
    receptor_types A C NA OA N SA HD     # receptor atom types
    ligand_types A C NA OA N HD          # ligand atom types
    receptor 1fpu_receptor_rigid.pdbqt   # macromolecule
    gridcenter 15.190 53.903 16.917      # xyz-coordinates or auto
    smooth 0.5                           # store minimum energy w/in rad(A)
    map 1fpu_receptor_rigid.A.map        # atom-specific affinity map
    map 1fpu_receptor_rigid.C.map        # atom-specific affinity map
    map 1fpu_receptor_rigid.NA.map       # atom-specific affinity map
    map 1fpu_receptor_rigid.OA.map       # atom-specific affinity map
    map 1fpu_receptor_rigid.N.map        # atom-specific affinity map
    map 1fpu_receptor_rigid.HD.map       # atom-specific affinity map
    elecmap 1fpu_receptor_rigid.e.map    # electrostatic potential map
    dsolvmap 1fpu_receptor_rigid.d.map   # desolvation potential map
    dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

To execute ``autogrid4`` using ``1fpu_receptor_rigid.gpf``, run the folllowing command line:

.. code-block:: bash

    $ autogrid4 -p 1fpu_receptor_rigid.gpf -l 1fpu_receptor_rigid.glg

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

While using the AutoDock4 forcefield, only the flex part of the receptor is necessary, as well as the affinity maps. Once the receptor (flex part ``1fpu_receptor_flex.pdbqt``), ligand ``1iep_ligand.pdbqt`` and maps ``1fpu_receptor_rigid`` were prepared, you can perform the flexible side-chain docking by simply running the following command line:

.. code-block:: bash

    $ vina --flex 1fpu_receptor_flex.pdbqt --ligand 1iep_ligand.pdbqt \
           --maps 1fpu_receptor_rigid --scoring ad4 \
           --exhaustiveness 32 --out 1fpu_ligand_flex_ad4_out.pdbqt

Running AutoDock Vina will write a PDBQT file called ``1fpu_ligand_flex_ad4_out.pdbqt`` contaning all the poses found during the molecular docking as well as the Thr315 sidechain conformations, and also present docking information to the terminal window.

4.b. Using Vina forcefield
__________________________

As well as for the fully rigid molecular docking, you only need to specify the center and dimensions (in Angstrom) of the grid. Here, instead of specifying each parameters for the grid box using the arguments ``--center_x, --center_y, --center_z`` and ``--size_x, --size_y, --size_z``, we will also store all those informations in a text file ``1fpu_receptor_rigid_vina_box.txt``.

.. code-block:: console
    :caption: Content of the config file (**1fpu_receptor_rigid_vina_box.txt**) for AutoDock Vina

    center_x = 15.190
    center_y = 53.903
    center_z = 16.917
    size_x = 20.0
    size_y = 20.0
    size_z = 20.0

However, when using the Vina forcefield, you will need to specify both the rigid ``1fpu_receptor_rigid.pdbqt`` (needed to compute internally the affinity maps) and flex part ``1fpu_receptor_flex.pdbqt`` of the receptor. To perform the same docking experiment but using Vina forcefield run the following command line:

.. code-block:: bash

    $ vina --receptor 1fpu_receptor_rigid.pdbqt --flex 1fpu_receptor_flex.pdbqt \
           --ligand 1iep_ligand.pdbqt --config 1fpu_receptor_rigid_vina_box.txt \
           --exhaustiveness 32 --out 1fpu_ligand_flex_vina_out.pdbqt

.. tip::

    Alternatively, you can use the Vinardo forcefield by adding the ``--scoring vinardo`` option.

Running AutoDock Vina will write a PDBQT file called ``1fpu_ligand_flex_vina_out.pdbqt``.

5. Results
----------

.. warning::
    
    Please don't forget that energy scores giving by the AutoDock and Vina forcefield are not comparable between each other.

5.a. Using AutoDock forcefield
______________________________

The predicted free energy of binding should be about ``-15 kcal/mol`` for poses that are similar to the crystallographic pose.

.. code-block:: console

    Scoring function : ad4
    Flex receptor: 1fpu_receptor_flex.pdbqt
    Ligand: ../data/1iep_ligand.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: -1132104431) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -15.41          0          0
       2       -14.95      1.164      1.803
       3       -13.92      1.112      1.744
       4       -13.39      3.975      6.038
       5       -13.08       1.48      2.166
       6       -12.13      3.877      11.74
       7       -12.13      5.806      9.094
       8       -11.89      1.251      1.971
       9       -11.55      2.804      10.81

5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-12 kcal/mol``.

.. code-block:: console

    Scoring function : vina
    Rigid receptor: 1fpu_receptor_rigid.pdbqt
    Flex receptor: 1fpu_receptor_flex.pdbqt
    Ligand: ../data/1iep_ligand.pdbqt
    Center: X 15.19 Y 53.903 Z 16.917
    Size: X 20 Y 20 Z 20
    Grid space: 0.375
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Computing Vina grid ... done.
    Performing docking (random seed: 1973662971) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -12.17          0          0
       2       -11.41       3.23       12.1
       3       -11.22      1.512      2.137
       4       -11.19       4.07         12
       5       -10.64      3.833      11.99
       6        -10.2      2.537      12.12
       7       -9.547      2.493      12.26
       8       -9.367      2.476      12.41
       9       -9.051      3.809      11.72
