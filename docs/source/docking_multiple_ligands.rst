.. _multiple_ligands_docking:

Multiple ligands docking
========================

Vina is now able to dock simultaneously multiple ligands. This functionality may find application in fragment based drug design, where small molecules that bind the same target can be grown or combined into larger compounds with potentially better affinity.

The protein PDE in complex with two inhibitors (pdb id: `5x72 <https://www.rcsb.org/structure/5X72>`_) was used as an example to demonstrate the ability of the AutoDock Vina to dock successfully multiple ligands. The two inhibitors in this structure are stereoisomers, and only the R-isomer is able to bind in a specific region of the pocket, while both the R- and S-isomers can bind to the second location. 


**System and software requirements**

This is a command-line-based tutorial for a basic docking experiment with AutoDock-Vina. It can be done on macOS, Linux, and Windows Subsystem for Linux (WSL). 

This tutorial uses python package Meeko for receptor and ligand preparation. Installation guide and advanced usage can be found from the documentation: `https://meeko.readthedocs.io/en/readthedocs/ <https://meeko.readthedocs.io/en/readthedocs/>`_.

**Input and expected output files**

The input and expected output files can be found here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/_basic_docking>`_.


1. Preparing receptor
----------------------------------

Exactly like the :ref:`basic_docking` tutorial, the method requires a receptor coordinate file that includes all hydrogen atoms. The file ``5x72_receptorH.pdb`` is provided. This file contains the receptor coordinates taken from the PDB entry ``5x72``. It was manually obtained by extracting the receptor coordinates (using an text editor) from the original PDB file ``5x72.pdb`` in the ``data`` directory, and the hydrogen atoms added using `reduce <https://github.com/rlabduke/reduce>`_. To create the PDBQT file and the vina box TXT (configuration) file, we will use the ``mk_prepare_receptor.py`` command-line script: 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 5x72_receptorH.pdb -o 5x72_receptor -p -v \
    --box_center -15.000 15.000 129.000 --box_size 30 24 24 \                     
    --default_altloc A -a

Here again we're using the ``-a`` option to ignore the partially resolve residue, LYS29 in chain A. Alternatively, it can be conviniently removed by ``--delete_residues A:29``. 

2. Prepare ligands
------------------

Here, we will prepare two ligands instead of only one. We will start from the SDF files ``5x72_ligand_p59.sdf`` and ``5x72_ligand_p69.sdf`` located in the ``data`` directory. They were also obtained directly from the `PDB <https://www.rcsb.org>`_ here: `5x72 <https://www.rcsb.org/structure/5X72>`_ (see ``Download instance Coordinates`` link for the P59 and P69 molecules). Since the ligand files do not include the hydrogen atoms, we are going to add them using ``scrub.py`` from python package Scrubber.

.. code-block:: bash

    $ scrub.py 5x72_ligand_p59.sdf -o 5x72_ligand_p59H.sdf 
    $ mk_prepare_ligand.py -i 5x72_ligand_p59H.sdf -o 5x72_ligand_p59.pdbqt
    $ scrub.py 5x72_ligand_p69.sdf -o 5x72_ligand_p69H.sdf 
    $ mk_prepare_ligand.py -i 5x72_ligand_p69H.sdf -o 5x72_ligand_p69.pdbqt


3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

Similar to what's done in :ref:`basic_docking`, to use the AutoDock FF for docking calculation, we need to generate a GPF file to precalculate the affinity map for each atom types. 

To prepare the GPF file for the rigid part of the receptor, you can run (or rerun) ``mk_prepare_receptor.py`` with the additional option, ``-g`` that will enable the writing of the GPF file. 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 5x72_receptorH.pdb -o 5x72_receptor -p -v -g \
    --box_center -15.000 15.000 129.000 --box_size 30 24 24 \                     
    --default_altloc A -a

After creating the GPF file, and now we can use the ``autogrid4`` command to generate the different map files that will be used for the molecular docking: 

.. code-block:: console
    :caption: Content of the grid parameter file (**5x72_receptor.gpf**) for the receptor (**5x72_receptor.pdbqt**)
    parameter_file boron-silicon-atom_par.dat
    npts 80 64 64
    gridfld 5x72_receptor.maps.fld
    spacing 0.375
    receptor_types HD C A N NA OA F P SA S Cl Br I Mg Ca Mn Fe Zn
    ligand_types HD C A N NA OA F P SA S Cl CL Br BR I Si B
    receptor 5x72_receptor.pdbqt
    gridcenter -15.000 15.000 129.000
    smooth 0.500
    map 5x72_receptor.HD.map
    map 5x72_receptor.C.map
    map 5x72_receptor.A.map
    map 5x72_receptor.N.map
    map 5x72_receptor.NA.map
    map 5x72_receptor.OA.map
    map 5x72_receptor.F.map
    map 5x72_receptor.P.map
    map 5x72_receptor.SA.map
    map 5x72_receptor.S.map
    map 5x72_receptor.Cl.map
    map 5x72_receptor.CL.map
    map 5x72_receptor.Br.map
    map 5x72_receptor.BR.map
    map 5x72_receptor.I.map
    map 5x72_receptor.Si.map
    map 5x72_receptor.B.map
    elecmap 5x72_receptor.e.map
    dsolvmap 5x72_receptor.d.map
    dielectric -42.000

To execute ``autogrid4`` using ``5x72_receptor.gpf``, run the folllowing command line:

.. code-block:: bash

    $ autogrid4 -p 5x72_receptor.gpf -l 5x72_receptor_rigid.glg

You should obtain as well the following files:

.. code-block:: console

    5x72_receptor.maps.fld       # grid data file
    5x72_receptor.*.map          # affinity maps for A, C, HD, NA, N, OA atom types
    5x72_receptor.d.map          # desolvation map
    5x72_receptor.e.map          # electrostatic map

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

Contrary to AutoDock4, you don't need to precalculate the affinity grid maps with ``autogrid4`` when using the Vina forcefield. AutoDock Vina computes those maps internally before the docking. If you did not make the box dimension file when preparing receptor in the previous step, you could specify the center and dimensions (in Angstrom) of the grid box in a new TXT file:  

.. code-block:: console
    :caption: Content of the config file (**5x72_receptor.box.txt**) for AutoDock Vina

    center_x = 15.190
    center_y = 53.903
    center_z = 16.917
    size_x = 20.0
    size_y = 20.0
    size_z = 20.0

And then run the following command to execute the docking calculation: 

.. code-block:: bash

    $ vina --receptor 5x72_receptor.pdbqt --ligand 5x72_ligand_p59.pdbqt 5x72_ligand_p69.pdbqt \
           --config 5x72_receptor.box.txt \
           --exhaustiveness=32 --out 5x72_ligand_vina_out.pdbqt

.. tip::

    Alternatively, you can use the Vinardo forcefield by adding the ``--scoring vinardo`` option.

Running AutoDock Vina will write a PDBQT file called ``5x72_ligand_vina_out.pdbqt``.

5. Results
----------

.. warning::
    
    Please don't forget that energy scores giving by the AutoDock and Vina forcefield are not comparable between each other.

5.a. Using AutoDock forcefield
______________________________

The predicted free energy of binding should be near ``-18 kcal/mol`` for poses that are similar to the crystallographic pose. Using the AutoDock4 scoring function, the first two sets of poses (top 2) need to be considered to show also a good overlap with the crystallographic poses
 
.. code-block:: console

    Scoring function : ad4
    Ligands:
    - 5x72_ligand_p59.pdbqt
    - 5x72_ligand_p69.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: -1370364650) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -17.67          0          0
    2       -17.61      1.124      3.731
    3       -17.45      1.837      3.718
    4       -17.41      1.981      9.343
    5       -17.17      1.242      3.802
    6       -17.17      1.436      9.123
    7       -17.11      1.478       5.26
    8        -17.1       1.62      8.954
    9          -17      1.669       9.66

5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-19 kcal/mol``. Using the Vina scoring function, the best set of poses (top 1) shows an excellent overlap with the crystallographic coordinates for one of the isomers.

.. code-block:: console

    Scoring function : vina
    Rigid receptor: 5x72_receptor.pdbqt
    Ligands:
    - 5x72_ligand_p59.pdbqt
    - 5x72_ligand_p69.pdbqt
    Grid center: X -15 Y 15 Z 129
    Grid size  : X 30 Y 24 Z 24
    Grid space : 0.375
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Computing Vina grid ... done.
    Performing docking (random seed: -1632509975) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -19.04          0          0
    2       -18.33       1.22       3.81
    3       -17.27      1.247      3.007
    4       -17.22      1.432      3.286
    5       -16.45      1.099      3.717
    6       -16.35        1.7      4.839
    7       -16.24      1.335      5.195
    8          -16      2.332      9.449
    9       -15.29      7.079       13.5
