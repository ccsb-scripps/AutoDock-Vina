#!/usr/bin/env python

from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from ringtail import RingtailCore

import argparse
import json
import logging
import multiprocessing
from os import linesep
import pathlib
import sys

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolInterchange

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.MolProps |
                                Chem.PropertyPickleOptions.PrivateProps)
RDLogger.DisableLog("rdApp.*")
rt_logger = logging.getLogger("ringtail")


class MolSupplier:
    """wraps other suppliers (e.g. Chem.SDMolSupplier) to change non-integer
        molecule names to integers, and to set rdkit mol names from properties
    """

    def __init__(self, supplier, name_from_prop=None, rename_to_int=False, nr_digits=10):
        self.supplier = supplier
        self.name_from_prop = name_from_prop
        self.rename_to_int = rename_to_int
        self.nr_digits = nr_digits
        self.names = {}
        self.counter = 0
        
    def __iter__(self):
        self.supplier.reset()
        return self

    def __next__(self):
        mol = self.supplier.__next__()
        if mol is None:
            return mol
        if self.name_from_prop:
            name = mol.GetProp(self.name_from_prop)
            mol.SetProp("_Name", name)
        if self.rename_to_int:
            name = mol.GetProp("_Name")
            newname = self._rename(name)
            mol.SetProp("_Name", newname)
        return mol
        
    def _rename(self, name):
        """rename if name is not an integer, or a sequence of alphabet chars
            followed by an integer."""

        # special case for Enamine's molecules
        if name.startswith("PV-") and name[3:].isdigit():
            return "PV" + name[3:] # remove dash from Enamine's PV-000000000000
        is_good = False
        if name.isalnum():
            # make sure all letters preceed the decimals, no mix
            is_good = True
            num_started = False
            for c in name:
                num_started |= c.isdecimal()
                if num_started and not c.isdecimal():
                    is_good = False
                    break
        if is_good:
            return name
        
        self.counter += 1
        #if name in self.names:
        #    raise RuntimeError("repeated molecule name: %s" % name)
        #self.names[name] = self.counter
        self.names[self.counter] = name
        tmp = "RN%0" + "%d" % self.nr_digits + "d"
        return tmp % self.counter


parser_essential = argparse.ArgumentParser(description="Run vina from SDF to SQLite", add_help=False)

parser_essential.add_argument("-l", "--ligands", help="input filename (.sdf)", required=True)
parser_essential.add_argument("-r", "--receptor", help="receptor filename (.pdbqt)")
parser_essential.add_argument("-m", "--maps", help="base filename of grid maps")
parser_essential.add_argument(      "--scoring_function", choices=["vina", "vinardo", "ad4"], default="vina")
parser_essential.add_argument("--out_db", help="output sqlite3 filename (.sqlite3/.db)", required=True)
parser_essential.add_argument("--size", help="size of search space (grid maps)", type=float, nargs=3)
parser_essential.add_argument("--center", help="center of search space (grid maps)", type=float, nargs=3)

basic = parser_essential.add_argument_group("options")
basic.add_argument("--out_sdf", help="output SD filename")
basic.add_argument("--name_from_prop", help="set molecule name from RDKit/SDF property")

args = parser_essential.parse_args()

v = Vina(sf_name=args.scoring_function, cpu=1)
v.set_receptor(args.receptor)
v.compute_vina_maps(center=args.center, box_size=args.size)
supplier = Chem.SDMolSupplier(args.ligands, removeHs=False)
if args.name_from_prop:
    supplier = MolSupplier(supplier, name_from_prop=args.name_from_prop)

# prepare ligand pdbqt
meeko_prep = MoleculePreparation()

rtc = RingtailCore(args.out_db)
rtc.save_receptor(args.receptor)

rt_logger.setLevel("WARNING")

def dock(mol):
    molsetups = meeko_prep.prepare(mol)
    if len(molsetups) != 1:
        return None
    lig_pdbqt = meeko_prep.write_pdbqt_string()
    v.set_ligand_from_string(lig_pdbqt)
    v.dock(max_evals=128000)
    output_pdbqt = v.poses(2)
    mol_name = mol.GetProp("_Name")
    vina_strings = {mol_name: output_pdbqt}
    return vina_strings
    
nr_cores = multiprocessing.cpu_count()
pool = multiprocessing.Pool(nr_cores - 1) # leave 1 for ringtail

for vina_strings in pool.imap_unordered(dock, supplier):
    rtc.add_results_from_vina_string(
            results_strings=vina_strings,
            save_receptor=False,
            add_interactions=True,
            )
