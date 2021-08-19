.. _python_docking:

Python scripting
================

Materials for this tutorial
---------------------------

For this tutorial, all the basic material are provided and can be found in the ``AutoDock-Vina/example/python_scripting`` directory (or on `GitHub <https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/python_scripting>`_).

First example
-------------

Let's begin with a very simple example of an AutoDock Vina script. It loads affinity maps for the AutoDock forcefield with their filenames starting by `1iep`, then load a PDBQT file containing a ligand called `1iep_ligand.pdbqt`, score the current pose, minimize it, dock it and save the minimized pose to a PDBQT file called `1iep_ligand_minimized.pdbqt`. For the full documentation about the `Vina` Python package, see :ref:`Python documentation<vina>`.

.. code-block:: python

    #! /usr/bin/env python
    # -*- coding: utf-8 -*-
    #
    # My first example with AutoDock Vina in python
    #

    from vina import Vina


    v = Vina(sf_name='vina')

    v.set_receptor('1iep_receptor.pdbqt')

    v.set_ligand_from_file('1iep_ligand.pdbqt')
    v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
    
    # Score the current pose
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

    # Dock the ligand
    v.dock(exhaustiveness=32, n_poses=20)
    v.write_poses('1iep_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)

You can find this script in the `example` folder of AutoDock-Vina available on Github. To execute it from a command line, go to your terminal/console/command prompt window. Navigate to the `examples` folder by typing

.. code-block:: console

    $ cd <examples_directory>/python_scripting
    $ python first_example.py

Let's go through the script line by line and see how it works.

.. code-block:: python

    from vina import Vina

This line is just for Python to tell it to load the Python package ``vina``.

.. code-block:: python

    v = Vina(sf_name='vina')

This line creates the ``Vina`` object that will be use for the subsequent taks. More precisely, it specifies also the forcefield that will be used for the molecular docking, here the Vina forcefield will be used. In the case, you want to use another forcefield like AutoDock4 or Vinardo, just replace ``vina`` by ``ad4`` or ``vinardo``. Other options can be define, like the number of cpu that will be used during the molecule docking. By default, Vina will use the avaialble cpus on you machine (``cpu=0``), but you could use only one with ``cpu=1``.

.. code-block:: python

    v.set_receptor('1iep_receptor.pdbqt')

Here, we are loading a PDBQT file called ``1iep_receptor.pdbqt`` containing the receptor. If necessary, PDBQT file containing the flexible sidechains can be also load at the same time by doing ``v.set_receptor('1iep_rigid.pdbqt', '1iep_flex.pdbqt')``.

.. code-block:: python

    v.set_ligand_from_file('1iep_ligand.pdbqt')
    v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])

The next lines are used to first load a PDBQT file containing the ligand called ``1iep_ligand.pdbqt`` and then compute the affinity maps for each ligand atom types accroding to the Vina forcefield. You might need to read first the tutorial :ref:`basic_docking` to learn how to create a PDBQT file of a ligand. There is a small subility here, the behavior of the ``compute_vina_maps()`` function changes if the ligand was loaded before or after computing the vina maps. If no ligand was initialized, ``compute_vina_maps()`` will compute the affinity map for each atom types defined in the Vina forcefield (22 in total). This is very useful when we want to dock ligands in batch (a.k.a virtual screening) but we don't necessarily know beforehand what atom types will be necessary for thoses ligands. Alternately to ``set_ligand_from_file()``, you could also load a molecule using a molecule string in PDBQT format using the ``set_ligand_from_string()`` function.

.. code-block:: python

    # Score the current pose
    energy = v.score()
    print('Score before minimization: %.3f (kcal/mol)' % energy[0])

Next, we simply ask Vina to calculate the energy (`score`) of the current pose using the forcefield defined at the beginning, and retrieve the energy of each component in a numpy array. In this case, we are going to print to the screen only the total energy of the current pose. This task is often useful when you want to get the energy from the specific pose.

.. code-block:: python

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

This line tells AutoDock Vina to perform a local energy minimization and show the total energy. It is useful sometimes to perform a quick energy minization after manually placing a ligand in a pocket and to remove possible steric clashes with itself and the receptor.

.. code-block:: python

    # Dock the ligand
    v.dock(exhaustiveness=32, n_poses=20)
    v.write_poses('1iep_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)

Finally, we run the molecular docking. Here we will ask `Vina` to run 32 consecutive Monte-Carlo samplings using the ``exhaustiveness`` argument and store 20 poses (``n_poses``) during the search. At the end, we will write a PDBQT file called ``1iep_ligand_vina_out.pdbqt`` containing only the 5 first poses (``n_poses``), ranked by score. Of course, this can be change to 20 to include all the poses that were saved during the calculations, at the condition that the energy difference between the best pose and the 20th pose if less than 3 kcal/mol. This behavior can be changed using the ``energy_range`` argument to an higher value.

