.. _docking_in_batch:

Docking in batch mode
================

Docking many ligands sequentially, one after the other, is possible with AutoDock Vina. This is useful to perform virtual screening of a large number of ligands against a receptor. The following sections describe how to do this using the command line. For possible ways of doing this with the Python API, see :ref:`python_dockingn`. 

Do not confuse this with multiple ligand docking (:ref:`multiple_ligands_docking`), in which multiple ligands are docked simultaneously. 

Command line options
--------------------

The command line usage is similar to the one for a single ligand. The only difference is that instead of specifying a single ligand using the ``--ligand`` option, we use the ``--batch`` option once or multiple times to include all ligands. 

Use ``--dir`` to specify the output directory. The output filename will be ``<ligand_name>_out.pdbqt`` under the specified directory, where ``<ligand_name>`` is the name of the input ligand file without the path. When the name of the ligand is duplicated, an index will be added to the name to avoid the name collision. 

Assuming we have multiple ligand PDBQT files in the directory named ``ligands`` in the current path. For Linux and macOS, the command line would look like this: 

.. code-block:: bash

    $ vina --receptor 1iep_receptor.pdbqt --batch ligands/1iep_ligand_1.pdbqt --batch ligands/1iep_ligand_2.pdbqt --batch ligands/1iep_ligand_3.pdbqt --config config.txt --dir poses
    $ vina --receptor 1iep_receptor.pdbqt --batch ligands/1iep_ligand_*.pdbqt --config config.txt --dir poses

For Windows, more precautions should be used with the wildcard expansion and the type of the path separators: 

.. code-block:: bash

    $ vina --receptor 1iep_receptor.pdbqt --batch $(ls .\ligands\*.pdbqt | % {$_.FullName}) --config config.txt --dir poses

Looping over a list of receptors
--------------------------------

Currently, there is no option to iterate multiple receptors in a series of docking runs. However, it is possible to loop over a list of receptors and ligands using a shell script or a batch file. 

Assuming we have multiple receptor PDBQT files in the directory named ``receptors`` in the current path, the following is an example of a shell script for Linux and macOS:  

.. code-block:: bash

    #!/bin/bash
    for receptor in receptors/*.pdbqt; do
        receptor_name=$(basename "$receptor" .pdbqt)
        vina --receptor "$receptor" --batch ligands/1iep_ligand_*.pdbqt --config config.txt --dir "poses/$receptor_name"
    done

For Windows, a similar batch file or PowerShell script could be written to achieve the same result. The syntax may vary, but the logic of the routine remains the same. 
