#!/usr/bin/env python

import sys

from openbabel import openbabel as ob
from vina import Vina

ligand_pdbqt_filename = "1gkc/flex-xray.pdbqt"
receptor_pdbqt_filename = "1gkc/protein.pdbqt"

obconv = ob.OBConversion()
obconv.SetInFormat("pdbqt")
OBMol = ob.OBMol()
obconv.ReadFile(OBMol, ligand_pdbqt_filename)

v = Vina()

v.set_receptor(receptor_pdbqt_filename)
v.set_forcefield()
v.set_box(66.380, 30.749, 117.885, 30, 30, 30, 0.375)

v.set_ligand(OBMol)
v.score()
v.optimize()
v.score()
v.write_pose("1gkc/output_optimize.pdbqt", "global search")
print("")

v.set_ligand(OBMol)
v.compute_vina_grid()
v.global_search(30, 1.0)
v.write_results("1gkc/output.pdbqt", 20, 3)
print("")

