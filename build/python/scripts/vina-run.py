#!/usr/bin/env python

from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from meeko import Polymer
try:
    from ringtail import RingtailCore
    _got_ringtail = True
except ImportError as err:
    _got_ringtail = False
    _ringtail_import_err = err

import argparse
import json
import logging
import multiprocessing
from os import linesep
import pathlib
import sys
import tempfile

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


def parse_vina_box(text):
    center_x = None
    center_y = None
    center_z = None
    size_x = None
    size_y = None
    size_z = None
    spacing = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("center_x"):
            center_x = float(line.split("=")[1])
        elif line.startswith("center_y"):
            center_y = float(line.split("=")[1])
        elif line.startswith("center_z"):
            center_z = float(line.split("=")[1])
        elif line.startswith("size_x"):
            size_x = float(line.split("=")[1])
        elif line.startswith("size_y"):
            size_y = float(line.split("=")[1])
        elif line.startswith("size_z"):
            size_z = float(line.split("=")[1])
        elif line.startswith("spacing"):
            spacing = float(line.split("=")[1])
    center = (center_x, center_y, center_z)
    size = (size_x, size_y, size_z)
    return center, size, spacing

DEFAULT_SPACING = 0.375

parser = argparse.ArgumentParser(description="Run vina from SDF to SQLite")

parser.add_argument("-l", "--ligands", help="input filename (.sdf)", required=True)
parser.add_argument("-r", "--receptor", help="filename of Meeko Polymer serialized to JSON")
parser.add_argument("-m", "--maps", help="base filename of grid maps")
parser.add_argument(      "--scoring_function", choices=["vina", "vinardo", "ad4"], default="vina")
parser.add_argument("--out_db", help="output sqlite3 filename (.sqlite3/.db)")
parser.add_argument("--size", help="size of search space (grid maps)", type=float, nargs=3)
parser.add_argument("--center", help="center of search space (grid maps)", type=float, nargs=3)
parser.add_argument("--spacing", help=f"distance between grid points (default: {DEFAULT_SPACING} Angstrom)", type=float)
parser.add_argument('-b', '--vina_box', help="filename of vina config with box size and center")
parser.add_argument("--out_sdf", help="output SD filename")
parser.add_argument("--name_from_prop", help="set input molecule name from RDKit/SDF property")
parser.add_argument("--log_filename", help="write log to this filename")
parser.add_argument("--cpu", type=int, default=0)
args = parser.parse_args()

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.MolProps |
                                Chem.PropertyPickleOptions.PrivateProps)
RDLogger.DisableLog("rdApp.*")
root_logger = logging.getLogger()
root_logger.setLevel("INFO")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
if args.log_filename is not None:
    file_handler = logging.FileHandler(args.log_filename)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)
else:
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

# check there is output
if args.out_db is None and args.out_sdf is None:
    print("Use --out_db and/or --out_sdf")
    sys.exit(2)

def grid_usage_error():
    print("use both --center and --size, or --vina_box, or --maps")
    sys.exit(2)
    return

if (args.center is None) != (args.size is None):
    grid_usage_error()
if args.vina_box is not None and args.center is not None:
    grid_usage_error()
if args.vina_box is None and (args.center is None or args.size is None):
    grid_usage_error()
if args.maps is not None and (args.center is not None or args.vina_box is not None):
    grid_usage_error()
spacing = DEFAULT_SPACING
if args.vina_box is not None:
    with open(args.vina_box) as f:
        txt = f.read()
    center, size, spacing_from_vina_box = parse_vina_box(txt)
    if spacing_from_vina_box is not None:
        spacing = spacing_from_vina_box
    if args.spacing is not None:
        spacing = args.spacing
elif args.center is not None:
    center = args.center
    size = args.size
elif args.maps is not None:
    center = None
    size = None
else:
    print("logic error in determining where box size/center is coming from")
    sys.exit(1)

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
        with open(args.receptor) as f:
            json_str = f.read()
        polymer = Polymer.from_json(json_str)
        with tempfile.NamedTemporaryFile(mode="wt", suffix=".pdbqt") as tmp:
            rigid, flex_dict = PDBQTWriterLegacy.write_from_polymer(polymer)
            tmp.write(rigid)
            v.set_receptor(tmp.name)
            if flex_dict:
                raise NotImplementedError("receptor has flexres, which are not passed along yet")
        v.compute_vina_maps(center=center, box_size=size)

supplier = Chem.SDMolSupplier(args.ligands, removeHs=False)
if args.name_from_prop:
    supplier = MolSupplier(supplier, name_from_prop=args.name_from_prop)

# prepare ligand pdbqt
meeko_prep = MoleculePreparation()

if args.out_db is not None:
    if not _got_ringtail:
        raise ImportError from _ringtail_import_err
    rtc = RingtailCore(args.out_db)
    rtc.save_receptor(args.receptor)  # TODO JSON?
    rt_logger = logging.getLogger("ringtail")
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
        v.dock()
        output_pdbqt = v.poses()
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
            output_rdmol.SetDoubleProp("VinaScore", pdbqt_mol[0].score)
            w.write(output_rdmol)
    except Exception as error:
        root_logger.error(str(error))

if args.out_sdf is not None:
    w.close()
