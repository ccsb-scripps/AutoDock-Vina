.. _zinc_docking:

Docking with zinc metalloproteins
=================================

Zinc is present in a wide variety of proteins and is important in the metabolism of most organisms. Zinc metalloenzymes are therapeutically relevant targets in diseases such as cancer, heart disease, bacterial infection, and Alzheimer’s disease. In most cases a drug molecule targeting such enzymes establishes an interaction that coordinates with the zinc ion. Thus, accurate prediction of the interaction of ligands with zinc is an important aspect of computational docking and virtual screening against zinc containing proteins. 

The AutoDock4 force field was extended to include a specialized potential describing the interactions of zinc-coordinating ligands. This potential describes both the energetic and geometric components of the interaction. The new force field, named AutoDock4Zn, was calibrated on a data set of 292 crystal complexes containing zinc. Redocking experiments show that the force field provides significant improvement in performance in both free energy of binding estimation as well as in root-mean-square deviation from the crystal structure pose.

.. note::
    This tutorial requires a certain degree of familiarity with the command-line interface. Also, we assume that you installed the ADFR software suite as well as the meeko Python package.

.. note::

    Please cite this paper if you are using this protocol in your work:

    - Santos-Martins, D., Forli, S., Ramos, M. J., & Olson, A. J. (2014). AutoDock4Zn: an improved AutoDock force field for small-molecule docking to zinc metalloproteins. Journal of chemical information and modeling, 54(8), 2371-2379.

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/docking_with_zinc_metalloproteins/data`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/docking_with_zinc_metalloproteins>`_). If you ever feel lost, you can always take a look at the solution here: ``AutoDock-Vina/example/docking_with_zinc_metalloproteins/solution``. All the Python scripts used here (except for ``prepare_receptor`` and ``mk_prepare_ligand.py``) are located in the ``AutoDock-Vina/example/autodock_scripts`` directory, alternatively you can also find them here on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts>`_.

1. Preparing the receptor
-------------------------

During this step we will create the PDBQT file of the receptor using the PDB file called ``proteinH.pdb``, containing all the hydrogen atoms, and add the tetrahedral zinc pseudo atoms (``TZ``) around the zinc ion. TZ atoms represent the preferred position for tetrahedral coordination by the ligand. This file contains the receptor coordinates of chain A and B taken from the PDB entry `1s63 <https://www.rcsb.org/structure/1S63>`_.

To prepare the receptor, execute the following command lines:

.. code-block:: bash

    $ prepare_receptor -r protein.pdb -o protein.pdbqt
    $ pythonsh <script_directory>/zinc_pseudo.py -r protein.pdbqt -o protein_tz.pdbqt

The execution of these two commands should output these two messages. One informing us that the charge for the zinc ion was not set by ``prepare_receptor``. In this context, the message can be safely ignored since the ligand will interact preferentially with the zinc pseudo atoms (TZ). The PDBQT output files can be found in the ``solution`` directory.

.. code-block:: console

    Sorry, there are no Gasteiger parameters available for atom proteinH:B: ZN1001:ZN

The second message is telling us that only one zinc pseudo atom (TZ) was added to the receptor.

.. code-block:: console

    Wrote 1 TZ atoms on protein_tz.pdbqt.

2. Preparing the ligand
-----------------------

The second step consists to prepare the ligand, by converting the SDF file ``1s63_ligand.msdf`` to a PDBQT file readable by AutoDock Vina. As usual, we will use the ``mk_prepare_ligand.py`` Python script from ``Meeko`` (see installation instruction here: :ref:`docking_requirements`) for this task. For your convenience, the molecule file ``1s63_ligand.sdf`` is provided (see ``data`` directory). But you can obtain it directly from the `PDB <https://www.rcsb.org>`_ here: `1s63 <https://www.rcsb.org/structure/1S63>`_ (see ``Download instance Coordinates`` link for the 778 molecule. Since the ligand file does not include the hydrogen atoms, we are going to automatically add them and correct the protonation for a pH of 7.4.

.. code-block:: bash

    $ mk_prepare_ligand.py -i 1s63_ligand.sdf -o 1s63_ligand.pdbqt --add_hydrogen --pH 7.4

The output PDBQT  ``1s63_ligand.pdbqt`` can be found in the ``solution`` directory.

3. Generating affinity maps
---------------------------

The preparation script ``prepare_gpf4zn.py`` will be used to generate a special GPF file for docking with zinc pseudo atoms:

.. code-block:: bash

    $ pythonsh <script_directory>/prepare_gpf4zn.py -l 1s63_ligand.pdbqt -r protein_tz.pdbqt \
    -o protein_tz.gpf  -p npts=40,30,50 -p gridcenter=18,134,-1 \
    –p parameter_file=AD4Zn.dat   

The ``-p`` flag is used to set the box center (``gridcenter``) and size (``npts``) along with the ``parameter_file`` specific for this case. After execution, you should obtain a GPF file called ``protein_tz.gpf`` containing this:

.. code-block:: console

    npts 40 30 50                        # num.grid points in xyz
    parameter_file AD4Zn.dat             # force field default parameter file
    gridfld protein_tz.maps.fld          # grid_data_file
    spacing 0.375                        # spacing(A)
    receptor_types A C TZ NA ZN OA N P SA HD # receptor atom types
    ligand_types A C Cl NA OA N HD       # ligand atom types
    receptor protein_tz.pdbqt            # macromolecule
    gridcenter 18 134 -1                 # xyz-coordinates or auto
    smooth 0.5                           # store minimum energy w/in rad(A)
    map protein_tz.A.map                 # atom-specific affinity map
    map protein_tz.C.map                 # atom-specific affinity map
    map protein_tz.Cl.map                # atom-specific affinity map
    map protein_tz.NA.map                # atom-specific affinity map
    map protein_tz.OA.map                # atom-specific affinity map
    map protein_tz.N.map                 # atom-specific affinity map
    map protein_tz.HD.map                # atom-specific affinity map
    elecmap protein_tz.e.map             # electrostatic potential map
    dsolvmap protein_tz.d.map              # desolvation potential map
    dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant
    nbp_r_eps 0.25 23.2135 12 6 NA TZ
    nbp_r_eps 2.1   3.8453 12 6 OA Zn
    nbp_r_eps 2.25  7.5914 12 6 SA Zn
    nbp_r_eps 1.0   0.0    12 6 HD Zn
    nbp_r_eps 2.0   0.0060 12 6 NA Zn
    nbp_r_eps 2.0   0.2966 12 6  N Zn


The AutoDock4Zn forcefield is mostly defined by non bonded pairwise potentials which are written to the GPF file ``protein_tz.gpf`` in the form of ``npb_r_eps`` keywords. The file ``AD4Zn.dat`` includes the definition of the TZ atom type for the AutoDock forcefield. The keyword ``parameter_file`` in the GPF file specifies ``AD4Zn.dat`` as the forcefield to be used, so AutoGrid requires a local copy of it in the working directory. Alternatively, the keyword ``parameter_file`` in the GPF can point to the full or relative path where ``AD4Zn.dat`` is located.

.. warning::

    The behavior of the ``npb_r_eps`` keyword changed between autogrid 4.2.6 and 4.2.7. Be sure that you are using the latest version (AutoGrid 4.2.7.x.2019-07-11) of autogrid4 shipped with the ADFR Suite.

.. code-block:: bash

    $ autogrid4 -p protein_tz.gpf -o protein_tz.glg

At this stage, all forcefield information has been encoded in the affinity maps, and the remaining steps are the same as in the standard AutoDock protocol.

4. Running AutoDock Vina
------------------------

4.a. Using AutoDock4 forcefield
_______________________________

When using the AutoDock4 forcefield, you only need to provide the affinity maps and the ligand, while specifying that the forcefield used will be AutoDock4 using the option ``--scoring ad4``.

.. code-block:: bash

    $ vina --ligand 1s63_ligand.pdbqt --maps protein_tz --scoring ad4 \
           --exhaustiveness 32 --out 1s63_ligand_ad4_out.pdbqt

5. Results
----------

The predicted free energy of binding should be about ``-13 kcal/mol`` for the best pose and should corresponds to the crystallographic pose.

.. code-block:: console

    Scoring function : ad4
    Ligand: 1s63_ligand.pdbqt
    Exhaustiveness: 32
    CPU: 0
    Verbosity: 1

    Reading AD4.2 maps ... done.
    Performing docking (random seed: 1984557646) ... 
    0%   10   20   30   40   50   60   70   80   90   100%
    |----|----|----|----|----|----|----|----|----|----|
    ***************************************************

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1        -13.5          0          0
       2          -13      2.518      4.707
       3       -12.56      2.116      2.499
       4       -12.44      3.041      4.021
       5       -12.12      2.975      6.211
       6       -11.96      2.814      6.336
       7       -11.91      3.244      6.477
       8       -11.32      3.783      5.654
       9       -11.31      2.856      3.867
