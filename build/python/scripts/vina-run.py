#!/usr/bin/env python

from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
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


parser = argparse.ArgumentParser(description="Run vina from SDF to SQLite")

parser.add_argument("-l", "--ligands", help="input filename (.sdf)", required=True)
parser.add_argument("-r", "--receptor", help="receptor filename (.pdbqt)")
parser.add_argument("-m", "--maps", help="base filename of grid maps")
parser.add_argument(      "--scoring_function", choices=["vina", "vinardo", "ad4"], default="vina")
parser.add_argument("--out_db", help="output sqlite3 filename (.sqlite3/.db)")
parser.add_argument("--size", help="size of search space (grid maps)", type=float, nargs=3)
parser.add_argument("--center", help="center of search space (grid maps)", type=float, nargs=3)
parser.add_argument("--out_sdf", help="output SD filename")
parser.add_argument("--name_from_prop", help="set input molecule name from RDKit/SDF property")
parser.add_argument("--log_filename", help="write log to this filename")
parser.add_argument("--cpu", type=int, default=0)
args = parser.parse_args()

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.MolProps |
                                Chem.PropertyPickleOptions.PrivateProps)
RDLogger.DisableLog("rdApp.*")
rt_logger = logging.getLogger("ringtail")
if args.log_filename is not None:
    root_logger = logging.getLogger()
    root_logger.setLevel("INFO")
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler = logging.FileHandler(args.log_filename)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

# check there is output
if args.out_db is None and args.out_sdf is None:
    print("Use --out_db and/or --out_sdf")
    sys.exit(2)

# we always need receptor to insert in DB
if args.scoring_function == "ad4":
    if args.maps is None:
        print("AD4 scoring function requires --maps to be passed in")
        sys.exit(2)
    v = Vina(sf_name="ad4", cpu=1)
    v.load_maps(args.maps)
else:
    v = Vina(sf_name=args.scoring_function, cpu=1)
    if args.maps is not None:
        v.load_maps(args.maps)
    else:
        v.set_receptor(args.receptor)
        v.compute_vina_maps(center=args.center, box_size=args.size)

supplier = Chem.SDMolSupplier(args.ligands, removeHs=False)
if args.name_from_prop:
    supplier = MolSupplier(supplier, name_from_prop=args.name_from_prop)

# prepare ligand pdbqt
meeko_prep = MoleculePreparation()

if args.out_db is not None:
    rtc = RingtailCore(args.out_db)
    rtc.save_receptor(args.receptor)
    rt_logger.setLevel("WARNING")

def dock(mol):
    try:
        mol_name = mol.GetProp("_Name")
        molsetups = meeko_prep.prepare(mol)
        if len(molsetups) != 1:
            return None
        molsetup = molsetups[0]
        lig_pdbqt, is_ok, err = PDBQTWriterLegacy.write_string(molsetup) #, add_index_map=True, remove_smiles=True)
        if not is_ok:
            raise RuntimeError(f'ligand not ok {mol.GetProp("_Name")=}')
        v.set_ligand_from_string(lig_pdbqt)
        v.dock(max_evals=128000)
        output_pdbqt = v.poses(2)
        vina_strings = {mol_name: output_pdbqt}
        return vina_strings
    except Exception as error:
        return error
    
if args.cpu == 0:
    nr_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(nr_cores - 1) # leave 1 for writing
    map_fn = pool.imap_unordered
elif args.cpu > 1:
    pool = multiprocessing.Pool(args.cpu - 1) # leave 1 for writing
    map_fn = pool.imap_unordered 
elif args.cpu == 1:
    map_fn = map
else:
    print("args.cpu can't be negative")
    sys.exit(2)

if args.out_sdf is not None:
    w = Chem.SDWriter(args.out_sdf)

for vina_strings in map_fn(dock, supplier):
    if isinstance(vina_strings, Exception):
        root_logger.error(str(vina_strings))
    try:
        if args.out_db is not None:
            rtc.add_results_from_vina_string(
                results_strings=vina_strings,
                save_receptor=False,
                add_interactions=True,
            )
        if args.out_sdf is not None:
            mol_name = list(vina_strings.keys())[0]
            root_logger.info(f"writing {mol_name=} to {args.out_sdf}")
            pdbqt_mol = PDBQTMolecule(vina_strings[mol_name])
            output_rdmol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0] # ignore sidechains
            w.write(output_rdmol)
    except Exception as error:
        root_logger.error(str(error))

if args.out_sdf is not None:
    w.close()
