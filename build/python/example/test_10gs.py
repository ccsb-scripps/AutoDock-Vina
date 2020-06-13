#!/usr/bin/env python

import sys

from openbabel import openbabel as ob
from vina import Vina

ligand_pdbqt_filenames = ["10gs/flex-xray.pdbqt", "10gs/6gss_flex-xray.pdbqt"]
receptor_pdbqt_filename = "10gs/protein.pdbqt"

obconv = ob.OBConversion()
obconv.SetInFormat("pdbqt")
OBMol = ob.OBMol()
obconv.ReadFile(OBMol, ligand_pdbqt_filenames[0])

v = Vina()
print(dir(v))

v.set_obmol(OBMol)
o = v.get_obmol()
print(o)

v.set_receptor(receptor_pdbqt_filename)
v.set_forcefield()
v.set_box(10.472, 6.967, 27.934, 30, 30, 30, 0.375)

for i, ligand_pdbqt_filename in enumerate(ligand_pdbqt_filenames):
    v.set_ligand(ligand_pdbqt_filenames[i])
    v.score()
    v.optimize(0)
    v.score()
    v.write_pose("10gs/output_%s_optimize.pdbqt" % i, "random")
    print("")

for i, ligand_pdbqt_filename in enumerate(ligand_pdbqt_filenames):
    v.set_ligand(ligand_pdbqt_filename)
    v.compute_vina_grid() # We have to reocompute the grid in order to have the right atom types...
    v.global_search(30, 1.0)
    v.write_results("10gs/output_%s.pdbqt" % i, 20, 1)
    print("")
