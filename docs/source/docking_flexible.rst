.. _flexible_docking:

Flexible docking
================

The lack of receptor flexibility is arguably the greatest limitation in these types of docking methods. However, AutoDock Vina allows some limited flexibility of selected receptor side chains. In this tutorial, we will describe the cross-docking of the `imatinib molecule <https://en.wikipedia.org/wiki/Imatinib>`_ to c-Abl in PDB entry `1fpu <https://www.rcsb.org/structure/1FPU>`_, treating Thr315 as flexible. 


**System and software requirements**

This is a command-line-based tutorial for a basic docking experiment with AutoDock-Vina. It can be done on macOS, Linux, and Windows Subsystem for Linux (WSL). 

This tutorial uses python package Meeko for receptor and ligand preparation. Installation guide and advanced usage can be found from the documentation: `https://meeko.readthedocs.io/en/readthedocs/ <https://meeko.readthedocs.io/en/readthedocs/>`_.

**Input and expected output files**

The input and expected output files can be found here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/_basic_docking>`_.


1. Preparing the flexible receptor
----------------------------------

During this step, we are going to split the receptor coordinates into two PDBQT files: one for the rigid portion and one for the flexible side chains. As with the :ref:`basic_docking` tutorial, the method requires a receptor coordinate file that includes all hydrogen atoms. The file ``1fpu_receptorH.pdb`` is provided, and we will use the ``mk_prepare_receptor.py`` command-line script from Meeko: 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 1fpu_receptorH.pdb -o 1fpu_receptor -p -v \
    --box_size 20 20 20 --box_center 15.190 53.903 16.917 \
    -f A:315 -a

Similar to the usage in :ref:`basic_docking`, the command specifies that the provided ``1fpu_receptorH.pdb`` is the input file, and ``1fpu_receptor`` will be the basename of the output files. As requested by the ``-p`` option, a receptor PDBQT will be generated. And as requested by the ``-v`` option along with the box specification arguments ``--box_size`` and ``--box_center``, a TXT file and a box PDB file containing the box dimension will be generated. Additionally, ``-f A:315`` specifies the flexible residue in chain A and with residue ID 315. Lastly, ``-a`` is the optional to ignore residues that are only partially resolved and those fail to match with the templates. 

The two generated PDBQT files will be named as: 

.. code-block:: console

    1fpu_receptor_rigid.pdbqt           # rigid part
    1fpu_receptor_flex.pdbqt            # flexible sidechain of Thr315


2. Prepare ligand
-----------------

For the molecular docking with flexible sidechains, we will reuse the ligand file ``1iep_ligand.pdbqt`` from the previous tutorial :ref:`basic_docking`.

3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

As well as for the docking with a fully rigid receptor, we need to generate a GPF file to precalculate the affinity map for each atom types. However, instead of using the full receptor, affinity maps will be calculated only for the rigid part of the receptor (``1fpu_receptor_rigid.pdbqt``).

To prepare the GPF file for the rigid part of the receptor, you can run (or rerun) ``mk_prepare_receptor.py`` with the additional option, ``-g`` that will enable the writing of the GPF file. 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 1fpu_receptorH.pdb -o 1fpu_receptor -p -v -g \
    --box_size 20 20 20 --box_center 15.190 53.903 16.917 \
    -f A:315 -a

After creating the GPF file, and now we can use the ``autogrid4`` command to generate the different map files that will be used for the molecular docking:

.. code-block:: console
    :caption: Content of the grid parameter file (**1fpu_receptor_rigid.gpf**) for the receptor c-Abl parameter_file boron-silicon-atom_par.dat
    npts 52 52 52
    gridfld 1fpu_receptor_rigid.maps.fld
    spacing 0.375
    receptor_types HD C A N NA OA F P SA S Cl Br I Mg Ca Mn Fe Zn
    ligand_types HD C A N NA OA F P SA S Cl CL Br BR I Si B
    receptor 1fpu_receptor_rigid.pdbqt
    gridcenter 15.190 53.903 16.917
    smooth 0.500
    map 1fpu_receptor_rigid.HD.map
    map 1fpu_receptor_rigid.C.map
    map 1fpu_receptor_rigid.A.map
    map 1fpu_receptor_rigid.N.map
    map 1fpu_receptor_rigid.NA.map
    map 1fpu_receptor_rigid.OA.map
    map 1fpu_receptor_rigid.F.map
    map 1fpu_receptor_rigid.P.map
    map 1fpu_receptor_rigid.SA.map
    map 1fpu_receptor_rigid.S.map
    map 1fpu_receptor_rigid.Cl.map
    map 1fpu_receptor_rigid.CL.map
    map 1fpu_receptor_rigid.Br.map
    map 1fpu_receptor_rigid.BR.map
    map 1fpu_receptor_rigid.I.map
    map 1fpu_receptor_rigid.Si.map
    map 1fpu_receptor_rigid.B.map
    elecmap 1fpu_receptor_rigid.e.map
    dsolvmap 1fpu_receptor_rigid.d.map
    dielectric -42.000

Note that when the dielectric is negative, AutoGrid will use distance-dependent dielectric of Mehler and Solmajer regardless of the number. To execute ``autogrid4`` using ``1fpu_receptor_rigid.gpf``, run the folllowing command line:

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

The flexible-receptor docking calculation using the AutoDock4 forcefield will require the flex part of the receptor as well as the affinity maps. Once the receptor (flex part ``1fpu_receptor_flex.pdbqt``), ligand ``1iep_ligand.pdbqt`` and maps ``1fpu_receptor_rigid`` were prepared, you can perform the flexible side-chain docking by simply running the following command line:

.. code-block:: bash

    $ vina --flex 1fpu_receptor_flex.pdbqt --ligand 1iep_ligand.pdbqt \
           --maps 1fpu_receptor_rigid --scoring ad4 \
           --exhaustiveness 32 --out 1fpu_ligand_flex_ad4_out.pdbqt

Running AutoDock Vina will write a PDBQT file called ``1fpu_ligand_flex_ad4_out.pdbqt`` contaning all the poses found during the molecular docking as well as the Thr315 sidechain conformations, and also present docking information to the terminal window.

4.b. Using Vina forcefield
__________________________

As explained in :ref:`basic_docking`, AutoDock Vina computes those maps internally before the docking. Therfore, you may simply execute the docking calculation with: 

.. code-block:: bash

    $ vina --receptor 1fpu_receptor_rigid.pdbqt --flex 1fpu_receptor_flex.pdbqt \
           --ligand 1iep_ligand.pdbqt --config 1fpu_receptor.box.txt \ 
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

The predicted free energy of binding should be about ``-14 kcal/mol`` for poses that are similar to the crystallographic pose.

.. code-block:: console

    Scoring function : ad4
    Flex receptor: 1fpu_receptor_flex.pdbqt
    Ligand: 1iep_ligand.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: 711073774) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1        -14.2          0          0
    2          -14      1.163      1.777
    3        -13.4      1.182      1.616
    4       -11.92      1.521       2.19
    5       -11.76      1.963      3.295
    6       -11.68      2.872       10.9
    7       -11.13      3.933      10.75
    8       -10.99      3.702      11.82
    9       -10.72       2.08      11.06


5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score near ``-12 kcal/mol``.

.. code-block:: console

    Scoring function : vina
    Rigid receptor: 1fpu_receptor_rigid.pdbqt
    Flex receptor: 1fpu_receptor_flex.pdbqt
    Ligand: 1iep_ligand.pdbqt
    Grid center: X 15.19 Y 53.903 Z 16.917
    Grid size  : X 20 Y 20 Z 20
    Grid space : 0.375
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Computing Vina grid ... done.
    Performing docking (random seed: 1431646130) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -11.63          0          0
    2       -10.57      3.242      12.12
    3        -10.3      3.974       11.9
    4       -9.906      3.903      11.97
    5       -9.895      2.609      12.28
    6       -9.854      1.958      13.07
    7       -8.849      2.059      12.19
    8       -8.758      3.259      12.12
    9       -8.543      3.981      12.35
