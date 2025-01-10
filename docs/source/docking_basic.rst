.. _basic_docking:

Basic docking
=============

Let's start with our first example of docking, where the typical usage pattern would be to dock a single molecule into a rigid receptor. In this example we will dock the approved anticancer drug `imatinib <https://en.wikipedia.org/wiki/Imatinib>`_ (Gleevec; PDB entry `1iep <https://www.rcsb.org/structure/1IEP>`_) in the structure of c-Abl using AutoDock Vina. The target for this protocol is the kinase domain of the proto-oncogene tyrosine protein kinase c-Abl. The protein is an important target for cancer chemotherapy—in particular, the treatment of chronic myelogenous leukemia.

System and software requirements
--------

This is a command-line tutorial for a basic docking experiment with AutoDock-Vina. It can be done on macOS, Linux, and Windows Subsystem for Linux (WSL). 

This tutorial uses python package **Meeko for receptor and ligand preparation**. Installation guide and advanced usage can be found from the `Meeko documentation <https://meeko.readthedocs.io/en/release>`_.

The **input and expected output files** can be found here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/basic_docking>`_.

.. note:: 
    
    If you are using this tutorial for your works, please cite the following paper:

    - Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

1. Preparing the receptor
-------------------------

During this step, we will create a PDBQT file of our receptor containing only the polar hydrogen atoms as well as partial charges. For this step, we will use the ``mk_prepare_receptor.py`` command-line script: 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 1iep_receptorH.pdb -o 1iep_receptor -p -v \
    --box_size 20 20 20 --box_center 15.190 53.903 16.917

The command specifies that the provided ``1iep_receptorH.pdb`` is the input file, and ``1iep_receptor`` will be the basename of the output files. As requested by the ``-p`` option, a receptor PDBQT will be generated. And as requested by the ``-v`` option along with the box specification arguments ``--box_size`` and ``--box_center``, a TXT file and a box PDB file containing the box dimension will be generated. The TXT file can be used as the config file for the docking calculation. And the PDB file can be used to visualize the box in other programs such as PyMOL. 

The provided PDB file contains the receptor coordinates taken from PDB entry ``1iep``. If you are using experimental structures (for instance, from the `Protein Data Bank <https://www.rcsb.org>`_), you may want to remove waters, ligands, cofactors, ions deemed unnecessary for the docking. During receptor preparation with ``mk_prepare_receptor.py``, there is a ``--delete_residues`` option to conviniently ignore those unwanted components from the input structure. 

As an alternate to ``mk_prepare_receptor.py`` , you may use the ``prepare_receptor`` command tool from the ADFR Suite. As a prerequisite, a receptor coordinate file must contain all hydrogen atoms. Many tools exist to add missing hydrogen atoms to a protein, one popular choice would be to use `REDUCE <https://github.com/rlabduke/reduce>`_. 

.. code-block:: bash

    $ prepare_receptor -r 1iep_receptorH.pdb -o 1iep_receptor.pdbqt


2. Preparing the ligand
-----------------------

This step is very similar to the previous step. We will also create a PDBQT file from a ligand molecule file (preferably in the SDF format) using the ``Meeko`` python package. For convenience, the file ``1iep_ligand.sdf`` is provided. Alternatively, the 3D structure of the ligand is available from Protein Data Bank, PubChem and other online databases. 

.. warning::
  
  We strongly advice you against using PDB format for preparing small molecules, since it does not contain information about bond connections. Please don't forget to always check the protonation state of your molecules before docking. Your success can sometimes hang by just an hydrogen atom. ;-)

In case your starting ligand structure does not contain hydrogens, you may consider the command-line script ``scrub.py`` from the python package `Scrubber <https://github.com/forlilab/scrubber>`_ to protonate the ligand. If given the Smiles string of the ligand, ``scrub.py`` is also able to generate 3D conformers and enumerate tautomeric and protonation states of ligand. 

.. code-block:: bash

    $ mk_prepare_ligand.py -i 1iep_ligand.sdf -o 1iep_ligand.pdbqt

Other options are available for ``mk_prepare_ligand.py`` by typing ``mk_prepare_ligand.py --help``. 


3. (Optional) Generating affinity maps for AutoDock FF
------------------------------------------------------

To use the AutoDock FF in docking, we need an additional input file, which is the grid parameter file (GPF) for AutoGrid4, to create an affinity map file for each atom types. The GPF file specifies an AutoGrid4 calculation, including the size and location of the grid, the atom types that will be used, the coordinate file for the rigid receptor, and other parameters for calculation of the grids.

To prepare the gpf file for AutoGrid4, you can run (or rerun) ``mk_prepare_receptor.py`` with the additional option, ``-g`` that will enable the writing of the GPF file. 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 1iep_receptorH.pdb -o 1iep_receptor -p -v -g \
    --box_size 20 20 20 --box_center 15.190 53.903 16.917

After creating the GPF file, and now we can use the ``autogrid4`` command to generate the different map files that will be used for the molecular docking:

.. code-block:: bash

    $ autogrid4 -p 1iep_receptor.gpf -l 1iep_receptor.glg

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

Contrary to AutoDock4, you don't need to precalculate the affinity grid maps with ``autogrid4`` when using the Vina forcefield. AutoDock Vina computes those maps internally before the docking. If you did not make the box dimension file when preparing receptor in the previous step, you could specify the center and dimensions (in Angstrom) of the grid box in a new TXT file:  

.. code-block:: console
    :caption: Content of the config file (**1iep_receptor.box.txt**) for AutoDock Vina

    center_x = 15.190
    center_y = 53.903
    center_z = 16.917
    size_x = 20.0
    size_y = 20.0
    size_z = 20.0

And then run the following command to execute the docking calculation: 

.. code-block:: bash

    $ vina --receptor 1iep_receptor.pdbqt --ligand 1iep_ligand.pdbqt \
           --config 1iep_receptor.box.txt \
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
    Performing docking (random seed: 1045208650) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -14.72          0          0
    2       -14.63      0.862      1.051
    3       -13.12      1.152      1.877
    4        -11.7      4.989      11.38
    5       -11.44      3.619      11.51
    6       -11.39       1.36      2.222
    7       -11.21      3.773      12.06
    8       -10.71      2.043      13.49
    9       -10.41      1.748      2.955


5.b. Using Vina forcefield
__________________________

Using the vina forcefield, you should obtain a similar output from Vina with the best score around ``-13 kcal/mol``.

.. code-block:: console

    Scoring function : vina
    Rigid receptor: 1iep_receptor.pdbqt
    Ligand: 1iep_ligand.pdbqt
    Grid center: X 15.19 Y 53.903 Z 16.917
    Grid size  : X 20 Y 20 Z 20
    Grid space : 0.375
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Computing Vina grid ... done.
    Performing docking (random seed: -1622165383) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -13.23          0          0
    2       -11.29     0.9857      1.681
    3       -11.28      3.044      12.41
    4       -11.15      3.813      12.24
    5       -9.746      3.313      12.36
    6       -9.132      1.736      13.47
    7       -9.079      2.559      12.78
    8       -8.931      3.951      12.69
    9       -8.762      3.541      12.21
