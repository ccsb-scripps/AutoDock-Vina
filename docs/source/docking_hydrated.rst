.. _hydrated_docking:

Hydrated docking
================

In physiological environments, proteins and other biological structures are surrounded by water molecules. When a small-molecule binds to a protein, it must displace most of the waters occupying the binding cavity. However, rarely are all water molecules displaced. Some waters can be so strongly bound and conserved among similar proteins that from a ligand-docking perspective they are considered a part of the target structure, altering the binding site topography. 

Thus a method was developed that uses the existing version of AutoDock4 and now the new version AutoDock Vina 1.2.x but modifies the force field to model explicit bridging water molecules. In tests, this method has shown improvement in the prediction of bound conformations of small fragment molecules, such as those used in fragment-based drug discovery. The protocol can be summarized by those steps:

    1. The ligand is decorated with an ensemble of water molecules (represented by dummy atoms), which may or may not then contribute to the intermolecular interactions. 
    2. A modified AutoGrid map is then used during docking, giving a favorable score when the water is well placed and omitting the water if it overlaps with the receptor. 
    3. Finally, docked results are analyzing and poses are rescored using only water molecules that were retained.

In this tutorial, we are going to dock a fragment-size ligand (nicotine) with explicit water molecules in the acetylcholine binding protein (AChBP) structure (PDB entry `1uw6 <https://www.rcsb.org/structure/1UW6>`_). It was shown that the absence of water molecules could have a dramatic influence in docking performance leading to either inaccurate scoring and/or incorrect pose. With hydrated docking, fragment-sized ligands show a overall RMSD improvement.

.. note::

    This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the meeko Python package.

.. note::
    
    The materials present is this tutorial can be also found here: `https://www.nature.com/articles/nprot.2016.051 <https://www.nature.com/articles/nprot.2016.051>`_. If you are using this tutorial or this docking method for your work, you can cite the following papers:

    - Forli, S., & Olson, A. J. (2012). A force field with discrete displaceable waters and desolvation entropy for hydrated ligand docking. Journal of medicinal chemistry, 55(2), 623-638.
    - Forli, S., Huey, R., Pique, M. E., Sanner, M. F., Goodsell, D. S., & Olson, A. J. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. Nature protocols, 11(5), 905-919.

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/hydrated_docking/data`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/hydrated_docking>`_). If you ever feel lost, you can always take a look at the solution here: ``AutoDock-Vina/example/hydrated_docking/solution``. All the Python scripts used here (except for ``prepare_receptor`` and ``mk_prepare_ligand.py``) are located in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.

1. Preparing the receptor
-------------------------

The receptor can be prepared using the method described earlier in the following tutorials: :ref:`basic_docking` or :ref:`flexible_docking` if one wants to incorporate some sidechain flexibility. The file ``1uw6_receptorH.pdb`` is provided (see ``data`` directory located at ``<autodock-vina_directory>/example/hydrated_docking/``). This file contains the receptor coordinates of chain A and B taken from the PDB entry ``1uw6``.

.. code-block:: bash

    $ prepare_receptor -r 1uw6_receptorH.pdb -o 1uw6_receptor.pdbqt

The output PDBQT file ``1uw6_receptor.pdbqt`` is available in ``solution`` directory if necessary.

2. Preparing the ligand
-----------------------

For the hydrated docking, explicit water molecules (W atoms) must be added to the molecule. And for that, we will use ``Meeko`` (see installation instruction here: :ref:`docking_requirements`). For convenience, the molecule file ``1uw6_ligand.sdf`` is provided (see ``data`` directory). But you can obtain it directly from the `PDB <https://www.rcsb.org>`_ here: `1uw6 <https://www.rcsb.org/structure/1UW6>`_ (see ``Download instance Coordinates`` link for the NCT molecule (chain U [A])). Since the ligand file does not include the hydrogen atoms, we are going to automatically add them and correct the protonation for a pH of 7.4. The option ``-w`` is use to add explicit water molecule to the molecule.

.. warning::
  
  We strongly advice you against using PDB format for preparing small molecules, since it does not contain information about bond connections. Please don't forget to always check the protonation state of your molecules before docking. Your success can sometimes hang by just an hydrogen atom. ;-)

.. code-block:: bash
    
    $ mk_prepare_ligand.py -i 1uw6_ligand.sdf -o 1uw6_ligand.pdbqt --add_hydrogen --pH 7.4 -w

In total, 2 water molecules were added to the fragment. If you were not able to generate the ``1uw6_ligand.pdbqt`` file, you can look at the ``solution`` directory.

3. Generating affinity maps
---------------------------

As well as for the :ref:`basic_docking` or :ref:`flexible_docking` tutorials, we will also need to calculate the affinity maps for each atom types present in the ligand. However, this time we will also need to craft a special ``W`` affinity map for the water molecules attached to the ligand. This ``W`` affinity map is obtained by combining the ``OA`` and ``HD`` grid maps, therefore we are going to use the ``-p ligand_types='A,C,OA,N,HD'`` option to be sure those atom types are included in the GPF, in addition of the ligand atom types, while ignoring the ``W`` atom type:

.. code-block:: bash

    $ pythonsh <script_directory>/prepare_gpf.py -l 1uw6_ligand.pdbqt -r 1uw6_receptor.pdbqt -y \
               -p ligand_types='A,NA,C,HD,N,OA' \

The option ``-y`` specifies that we want to center automatically the grid around the ligand. After manually adding the ``OA`` atom type, you should have a GPF file called ``1uw6_receptor.gpf`` that looks like this:

.. code-block:: console
    :caption: Content of the grid parameter file (**1uw6_receptor.gpf**) for the receptor (**1uw6_receptor.pdbqt**)

    npts 40 40 40                        # num.grid points in xyz
    gridfld 1uw6_receptor.maps.fld       # grid_data_file
    spacing 0.375                        # spacing(A)
    receptor_types A C NA OA N SA HD     # receptor atom types
    ligand_types A NA C HD N OA          # ligand atom types
    receptor 1uw6_receptor.pdbqt         # macromolecule
    gridcenter 83.640 69.684 -10.124     # xyz-coordinates or auto
    smooth 0.5                           # store minimum energy w/in rad(A)
    map 1uw6_receptor.A.map              # atom-specific affinity map
    map 1uw6_receptor.NA.map             # atom-specific affinity map
    map 1uw6_receptor.C.map              # atom-specific affinity map
    map 1uw6_receptor.HD.map             # atom-specific affinity map * ADD OA IF NOT PRESENT *
    map 1uw6_receptor.N.map              # atom-specific affinity map
    map 1uw6_receptor.OA.map             # atom-specific affinity map * ADD OA IF NOT PRESENT *
    elecmap 1uw6_receptor.e.map          # electrostatic potential map
    dsolvmap 1uw6_receptor.d.map              # desolvation potential map
    dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant

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
      OA points : 91.73%
      HD points : 8.27%

      lowest  map value : -0.99
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
    Performing docking (random seed: -655217817) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1       -8.077          0          0
       2        -7.63      2.038      2.684
       3       -7.382      2.378      2.747
       4        -7.27      2.063      2.538
       5       -7.138      1.861      5.391
       6       -7.129          2      2.542
       7       -7.078      3.307      5.442
       8       -7.065       2.22      4.872
       9       -7.051      3.135      5.636

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


     receptor structure loaded           [ 4069 atoms ]
     receptor 5A shell extracted             [ 480 atoms in 5 A shell ] 
     removing ligand/ligand overlapping waters    [ 5 water(s) removed ]
     removing ligand/receptor overlapping waters      [ 8 water(s) removed ]

     scanning grid map for conserved waters...    [ filtered pose contains 5 waters ]

     water grid score results [ map: 1uw6_receptor.W.map ] 
         [ Water STRONG ( -0.92 ) +++ ]
         [ Water STRONG ( -0.66 ) +++ ]
         [ Water  WEAK  ( -0.50 )  +  ]
         [ Water STRONG ( -0.83 ) +++ ]
         [ Water STRONG ( -0.99 ) +++ ]

Waters are ranked (STRONG, WEAK) and scored inside the output file ``1uw6_ligand_ad4_out_DRY_SCORED.pdbqt`` with the calculated energy.
