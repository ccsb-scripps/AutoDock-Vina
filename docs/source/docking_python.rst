Python scripting
================

.. code-block:: python
	
	#! /usr/bin/env python

	from vina import Vina

	v = Vina(sf_name='ad4')
	v.load_maps('1iep')

	v.set_ligand('1iep_ligand.pdbqt')
	v.dock(exhaustiveness=32)

	v.write_poses('1iep_ligandH_out.pdbqt')

.. code-block:: console

	python basic_docking.py