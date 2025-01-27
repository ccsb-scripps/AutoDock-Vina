.. _hydrated_docking:

Hydrated docking
================

In physiological environments, proteins and other biological structures are surrounded by water molecules. When a small-molecule binds to a protein, it must displace most of the waters occupying the binding cavity. However, rarely are all water molecules displaced. Some waters can be so strongly bound and conserved among similar proteins that from a ligand-docking perspective they are considered a part of the target structure, altering the binding site topography. 

Thus a method was developed that uses the existing version of AutoDock4 and now the new version AutoDock Vina 1.2.x but modifies the force field to model explicit bridging water molecules. In tests, this method has shown improvement in the prediction of bound conformations of small fragment molecules, such as those used in fragment-based drug discovery. The protocol can be summarized by those steps:

    1. The ligand is decorated with an ensemble of water molecules (represented by dummy atoms), which may or may not then contribute to the intermolecular interactions. 
    2. A modified AutoGrid map is then used during docking, giving a favorable score when the water is well placed and omitting the water if it overlaps with the receptor. 
    3. Finally, docked results are analyzing and poses are rescored using only water molecules that were retained.

In this tutorial, we are going to dock a fragment-size ligand (nicotine) with explicit water molecules in the acetylcholine binding protein (AChBP) structure (PDB entry `1uw6 <https://www.rcsb.org/structure/1UW6>`_). It was shown that the absence of water molecules could have a dramatic influence in docking performance leading to either inaccurate scoring and/or incorrect pose. With hydrated docking, fragment-sized ligands show a overall RMSD improvement.

System and software requirements
--------

This is a command-line tutorial for a hydrated docking experiment with AutoDock-Vina. It can be done on macOS, Linux, and Windows Subsystem for Linux (WSL). 

This tutorial uses python package **Meeko for receptor and ligand preparation**. Installation guide and advanced usage can be found from the `Meeko documentation <https://meeko.readthedocs.io/en/release-doc>`_.

The **input and expected output files** can be found here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/hydrated_docking>`_. 

This tutorial requires the following **specialized Python scripts** for preparation and result processing: ``mapwater.py`` and ``dry.py``. These scripts are provided in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.

.. note:: 
    
    If you are using this tutorial or this docking method for your work, please cite the following papers:

    - Forli, S., & Olson, A. J. (2012). A force field with discrete displaceable waters and desolvation entropy for hydrated ligand docking. Journal of medicinal chemistry, 55(2), 623-638.
    - Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

1. Preparing the receptor
-------------------------

The receptor can be prepared using the method described earlier in the following tutorials: :ref:`basic_docking` or :ref:`flexible_docking` if one wants to incorporate some sidechain flexibility. The file ``1uw6_receptorH.pdb`` is provided. This file contains the receptor coordinates from the PDB entry ``1uw6``. To create the PDBQT file and the grid parameter file (GPF) for autogrid4, we will use the ``mk_prepare_receptor.py`` command-line script:  

.. code-block:: bash

    $ mk_prepare_receptor.py -i 1uw6_receptorH.pdb -o 1uw6_receptor -p -g \
    --box_center 83.640 69.684 -10.124 --box_size 15 15 15


2. Preparing the ligand
-----------------------

For the hydrated docking, explicit water molecules (W atoms) must be added to the molecule. And for that, we will use ``Meeko`` (see installation instruction here: :ref:`docking_requirements`). For convenience, the molecule file ``1uw6_ligand.sdf`` is provided (see ``data`` directory). But you can obtain it directly from the `PDB <https://www.rcsb.org>`_ here: `1uw6 <https://www.rcsb.org/structure/1UW6>`_ (see ``Download instance Coordinates`` link for the NCT molecule (chain U [A])). Since the ligand file does not include the hydrogen atoms, we are going to add them using ``scrub.py`` from python package Molscrub. The option ``-w`` is use to add explicit water molecule to the molecule.

.. warning::
  
  We strongly advice you against using PDB format for preparing small molecules, since it does not contain information about bond connections. Please don't forget to always check the protonation state of your molecules before docking. Your success can sometimes hang by just an hydrogen atom. ;-)

.. code-block:: bash
    
    $ scrub.py 1uw6_ligand.sdf -o 1uw6_ligandH.sdf
    $ mk_prepare_ligand.py -i 1uw6_ligandH.sdf -o 1uw6_ligand.pdbqt -w

In total, 2 water molecules were added to the fragment. 

3. Generating affinity maps
---------------------------

The hydrated docking method was calibrated and validated with the AutoDock4 forcefield. Therefore, we need to generate a GPF file to precalculate the affinity map for each atom types. 

In case you haven't already made the GPF file, you could rerun ``mk_prepare_receptor.py`` with the additional option, ``-g`` that will enable the writing of the GPF file. 

.. code-block:: bash
    
    $ mk_prepare_receptor.py -i 1uw6_receptorH.pdb -o 1uw6_receptor -p -g \
    --box_center 83.640 69.684 -10.124 --box_size 15 15 15

After creating the GPF file, and now we can use the ``autogrid4`` command to generate the different map files that will be used for the molecular docking: 

.. code-block:: console
    :caption: Content of the grid parameter file (**1uw6_receptor.gpf**) for the receptor (**1uw6_receptor.pdbqt**)

    parameter_file boron-silicon-atom_par.dat
    npts 40 40 40
    gridfld 1uw6_receptor.maps.fld
    spacing 0.375
    receptor_types HD C A N NA OA F P SA S Cl Br I Mg Ca Mn Fe Zn
    ligand_types HD C A N NA OA F P SA S Cl CL Br BR I Si B
    receptor 1uw6_receptor.pdbqt
    gridcenter 83.640 69.684 -10.124
    smooth 0.500
    map 1uw6_receptor.HD.map
    map 1uw6_receptor.C.map
    map 1uw6_receptor.A.map
    map 1uw6_receptor.N.map
    map 1uw6_receptor.NA.map
    map 1uw6_receptor.OA.map
    map 1uw6_receptor.F.map
    map 1uw6_receptor.P.map
    map 1uw6_receptor.SA.map
    map 1uw6_receptor.S.map
    map 1uw6_receptor.Cl.map
    map 1uw6_receptor.CL.map
    map 1uw6_receptor.Br.map
    map 1uw6_receptor.BR.map
    map 1uw6_receptor.I.map
    map 1uw6_receptor.Si.map
    map 1uw6_receptor.B.map
    elecmap 1uw6_receptor.e.map
    dsolvmap 1uw6_receptor.d.map
    dielectric -42.000

You can now execute ``autogrid4`` using the GPF file called ``1uw6_receptor.gpf`` and generate the additional water map ``W`` by combining ``OA`` and ``HD`` affinity maps using ``mapwater.py``:

.. code-block:: bash

    $ autogrid4 -p 1uw6_receptor.gpf -l 1uw6_receptor.glg
    $ python <script_directory>/mapwater.py -r 1uw6_receptor.pdbqt -s 1uw6_receptor.W.map

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
    OA points : 91.66%
    HD points : 8.34%

    lowest  map value : -0.98
    highest map value : -0.01

4. Running AutoDock Vina
------------------------

4.a. Using AutoDock4 forcefield
_______________________________

Now that you generated the ligand with explicit water molecules attached (``1uw6_ligand.pdbqt``) and the extra affinity map for the ``W`` atom type (``1uw6_receptor.W.map``), you can do the molecular docking with Vina using the AutoDock4 forcefield:

.. code-block:: bash

    $ vina  --ligand 1uw6_ligand.pdbqt --maps 1uw6_receptor --scoring ad4 \
            --exhaustiveness 32 --out 1uw6_ligand_ad4_out.pdbqt

4.b. Using Vina forcefield
__________________________

.. warning::
    
    While this method was calibrated and validated with the AutoDock4 forcefield, we strongly advice you against using this protocol with the Vina and Vinardo forcefield.

5. Results and post-processing
------------------------------

.. warning::

    Be aware that with this implementation of the method, it is difficult to compare results obtained with very diverse ligands without doing extra of post-processing on the results, because the energy estimation needs to be normalized. For this reason, the method is not suitable for virtual screenings. This doesn’t affect the structural accuracy, so comparisons within docking poses are fine. An improved scoring function to overcome this issue is in the works.

The predicted free energy of binding should be about ``-8 kcal/mol`` for poses that are similar to the crystallographic pose.

.. code-block:: console

    Scoring function : ad4
    Ligand: 1uw6_ligand.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: 1952347903) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -8.261          0          0
    2       -7.673      1.124      1.239
    3       -7.489      2.051       2.49
    4       -7.225      2.441      3.621
    5       -7.211      1.905      2.479
    6       -7.065      2.469       5.79
    7       -6.978      3.059      5.719
    8       -6.968      2.339      3.029
    9       -6.931      3.448      5.773

Docking results are filtered by using the receptor to remove displaced waters and the W map file to rank the conserved ones as strong or weak water molecules.

.. code-block:: bash

    $ python <script_directory>/dry.py -r 1uw6_receptor.pdbqt -m 1uw6_receptor.W.map -i 1uw6_ligand_ad4_out.pdbqt

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
    importing ATOMS from  1uw6_ligand_ad4_out.pdbqt

    [ using map file 1uw6_receptor.W.map ]
    ===============================================================


    receptor structure loaded	 		 [ 4069 atoms ]
    receptor 5A shell extracted  			 [ 485 atoms in 5 A shell ] 
    removing ligand/ligand overlapping waters	  [ 0 water(s) removed ]
    removing ligand/receptor overlapping waters	  [ 0 water(s) removed ]

    scanning grid map for conserved waters...	  [ filtered pose contains 18 waters ]

    water grid score results [ map: 1uw6_receptor.W.map ] 
        [ Water STRONG ( -0.92 ) +++ ]
        [ Water DISPLC ( -0.25 )  D  ]
        [ Water STRONG ( -0.89 ) +++ ]
        [ Water DISPLC ( -0.20 )  D  ]
        [ Water DISPLC ( -0.20 )  D  ]
        [ Water DISPLC ( -0.25 )  D  ]
        [ Water STRONG ( -0.65 ) +++ ]
        [ Water DISPLC ( -0.21 )  D  ]
        [ Water STRONG ( -0.92 ) +++ ]
        [ Water  WEAK  ( -0.32 )  +  ]
        [ Water  WEAK  ( -0.49 )  +  ]
        [ Water DISPLC ( -0.20 )  D  ]
        [ Water STRONG ( -0.53 ) +++ ]
        [ Water  WEAK  ( -0.39 )  +  ]
        [ Water STRONG ( -0.89 ) +++ ]
        [ Water  WEAK  ( -0.47 )  +  ]
        [ Water STRONG ( -0.81 ) +++ ]
        [ Water DISPLC ( -0.20 )  D  ]

Waters are ranked (STRONG, WEAK) and scored inside the output file ``1uw6_ligand_ad4_out_DRY_SCORED.pdbqt`` with the calculated energy.
