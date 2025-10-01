#!/usr/bin/env python

from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from meeko import Polymer
from meeko import gridbox
from meeko import pdbutils
try:
    from ringtail import RingtailCore
    _got_ringtail = True
except ImportError as err:
    _got_ringtail = False
    _ringtail_import_err = err

import argparse
import contextlib
import json
import logging
import multiprocessing
from os import linesep
from os import getcwd
from os import chdir
import pathlib
import subprocess
import sys
import tempfile
import shutil

import numpy as np 

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit import RDLogger
from rdkit.Chem import rdMolInterchange


logger = logging.getLogger(__name__)

@contextlib.contextmanager
def temporary_directory(suffix=None, prefix=None, dir=None, clean=True):
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = tempfile.mkdtemp(suffix, prefix, dir)
    cwd = getcwd()
    chdir(temp_dir)
    try:
        yield temp_dir
    finally:
        chdir(cwd)
        if clean:
            shutil.rmtree(temp_dir)

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

def get_parameter_text(vdw, hb, elec, dsolv):
    txt =  f"FE_coeff_vdW    {vdw:.4f}\n"
    txt += f"FE_coeff_hbond  {hb:.4f}\n"
    txt += f"FE_coeff_estat  {elec:.4f}\n"
    txt += f"FE_coeff_desolv {dsolv:.4f}\n"
    return txt

def create_gpf_dir(gpf_text, dest_folder, new_gpf_fn, vdw, hb, elec, dsolv):
    p = pathlib.Path(dest_folder)
    p.mkdir(exist_ok=True)
    weights_filename = "weights.dat"
    weights_text = get_parameter_text(vdw, hb, elec, dsolv)
    with open(p / weights_filename, "w") as f:
        f.write(weights_text)
    gpf_text = f"parameter_file {weights_filename}" + "\n" + gpf_text
    with open(p / new_gpf_fn, "w") as f:
        f.write(gpf_text)
    return

def wrap_autogrid(
    rec_path, box_center, box_size, dest_folder, grid_spacing, autogrid_path,
    rec_types, lig_types,
    vdw=0.1662, hb=0.1209, elec=0.1406, dsolv=0.1322,
):
    rec_fn = pathlib.Path(rec_path).name
    gpf_string, _npts = gridbox.get_gpf_string(
        box_center,
        box_size,
        rec_fn,
        rec_types,
        lig_types,
        dielectric=-42,
        smooth=0.5, 
        spacing=grid_spacing,
        ff_param_fname=None,
    )
    create_gpf_dir(gpf_string, dest_folder, "autogrid.gpf", vdw, hb, elec, dsolv)
    if len(pathlib.Path(rec_path).parents) > 1:
        shutil.copy(rec_path, str(pathlib.Path(dest_folder) / rec_fn))
    cmds = [autogrid_path, "-p", "autogrid.gpf", "-l", "autogrid.glg"]
    logger.info(f"subprocess run: {cmds}")
    o = subprocess.run(cmds, cwd=dest_folder, capture_output=True)
    logger.info(f"subprocess stdout: {o}")

    #fld_fn = [fn for fn in pathlib.Path(f"grids_{term.replace('ad4_', '')}").glob("*.maps.fld")]
    fld_fn = [fn for fn in pathlib.Path(dest_folder).glob("*.maps.fld")]
    if len(fld_fn) != 1:
        raise RuntimeError("expected 1 file eding with .maps.fld, got {len(fld_fn)=} {fld_fn=}")
    maps_fn = str(fld_fn[0]).replace(".maps.fld", "")
    return maps_fn

def _get_types_from_pdbqt(fname):
    atypes = set()
    with open(fname) as f:
        for line in f:
            is_atom = line.startswith("ATOM") or line.startswith("HETATM")
            if not is_atom:
                continue
            atype = line[77:].strip()
            atypes.add(atype)
    return atypes


def get_positions_from_molecule_file(filename):
    ext = filename.split('.')[-1]
    suppliers = {
        "pdb": None,  # overriden below, needed here as valid type
        "mol": Chem.MolFromMolFile,
        "mol2": Chem.MolFromMol2File,
        "sdf": Chem.SDMolSupplier,
        "pdbqt": None,
    }
    if ext not in suppliers.keys():
        print(f"File type given to --box_enveloping must be [.pdb/.mol/.mol2/.sdf/.pdbqt]")
        sys.exit(2)
    elif ext == ".pdb":
        pdbstr = pdbutils.strip_altloc_from_pdb_file(filename)
        mol = Chem.MolFromPDBBlock(pdbstr, removeHs=False, sanitize=False)
    elif ext == ".sdf":
        mol = next(suppliers[ext](filename, removeHs=False, sanitize=False))
    elif ext == ".pdbqt":
        pdbqtmol = PDBQTMolecule.from_file(filename)
        mol = RDKitMolCreate.from_pdbqt_mol(pdbqtmol)
    else:
        mol = suppliers[ext](filename, removeHs=False, sanitize=False)
    positions = mol.GetConformer.GetPositions() 
    return positions


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
DEFAULT_PADDING = 10

parser = argparse.ArgumentParser(description="Run vina from SDF to SQLite")

parser.add_argument("-l", "--ligands", help="input filename (.sdf)", required=True)
parser.add_argument("-r", "--receptor", help="filename of Meeko Polymer serialized to JSON")
parser.add_argument("-m", "--maps", help="base filename of grid maps")
parser.add_argument("-s", "--scoring", choices=["vina", "vinardo", "ad4"], default="vina")
# do not set default padding here at argparse level so we can test if the user passed a value
# deliberately with args.padding is None
parser.add_argument("--padding", help=f"Angstroms between box and atoms passed to --box_enveloping (default: {DEFAULT_PADDING})", type=float)
parser.add_argument("--box_enveloping", help="Box will envelop atoms in this file [.sdf .mol .mol2 .pdb .pdbqt]")
parser.add_argument("--flexible_amides", action="store_true")
parser.add_argument("--out_db", help="output sqlite3 filename (.sqlite3/.db)")
parser.add_argument("--size", help="size of search space (grid maps)", type=float, nargs=3)
parser.add_argument("--center", help="center of search space (grid maps)", type=float, nargs=3)
parser.add_argument("--spacing", help=f"distance between grid points (default: {DEFAULT_SPACING} Angstrom)", type=float)
parser.add_argument('-b', '--vina_box', help="filename of vina config with box size and center")
parser.add_argument("--out_sdf", help="output SD filename")
parser.add_argument("--name_from_prop", help="set input molecule name from RDKit/SDF property")
parser.add_argument("--log_filename", help="write log to this filename")
parser.add_argument("--vina_cpp_threads", type=int, default=1)
parser.add_argument("--nr_process", type=int, default=0)
parser.add_argument("--exhaustiveness", type=int)
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
    print("use one of these combinations:")
    print(f"    1) --box_enveloping                 (will pad with {DEFAULT_PADDING} Angstrom)")
    print(f"    2) --box_enveloping and --padding")
    print(f"    3) --box_enveloping and --size")
    print(f"    4) --vina_box")
    print(f"    5) --vina_box and --size            (to override size in --vina_box)")
    print(f"    6) --maps")
    print(f"    7) --center and --size")
    sys.exit(2)
    return

nr_box_options = 0
nr_box_options += int(args.box_enveloping is not None)
nr_box_options += int(args.maps is not None)
nr_box_options += int(args.vina_box is not None)
nr_box_options += int(args.center is not None)
if nr_box_options != 1:
    grid_usage_error()
if args.padding is not None and args.box_enveloping is None:
    print("--padding requires --box_enveloping")
    grid_usage_error()
if args.center is not None and args.size is None:
    print("--center requires --size")
    grid_usage_error()
if args.size is not None and (args.box_envelopping is None and args.center is None):
    print("--size requires either --center or --box_enveloping or --vina_box")
    grid_usage_error()
if args.size is not None and args.padding is not None:
    print("can't use both --size and --padding")
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
    if args.size is not None:
        size = args.size
elif args.center is not None:
    center = args.center
    size = args.size
elif args.box_enveloping is not None:
    padding = DEFAULT_PADDING if args.padding is None else args.padding
    positions = get_positions_from_molecule_file(args.box_enveloping)
    center, size = gridbox.calc_box(positions, padding)
    if args.size is not None:
        size = args.size
else:
    print("logic error in determining where box size/center is coming from, please report on github")
    sys.exit(1)

# we always need receptor to insert in DB
if args.scoring == "ad4":
    v = Vina(sf_name="ad4", cpu=args.vina_cpp_threads)
    if args.maps is None:
        ligtypes = ["HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S", "Cl", "Br", "I", "Si"]
        rec_fn = str(pathlib.Path(args.receptor).resolve())
        with open(args.receptor) as f:
            json_str = f.read()
        polymer = Polymer.from_json(json_str)
        with temporary_directory() as tmpdir:
            pdbqt_tuple = PDBQTWriterLegacy.write_from_polymer(polymer)
            rigid_pdbqt, flex_dict = pdbqt_tuple
            if flex_dict:
                raise NotImplementedError("receptor has flexres, which are not passed along yet")
            with open("receptor.pdbqt", "w") as f:
                f.write(rigid_pdbqt)
            rectypes = _get_types_from_pdbqt("receptor.pdbqt")
            maps_fn = wrap_autogrid(
                "receptor.pdbqt",
                center,
                size,
                tmpdir,
                spacing,
                "autogrid4",
                rec_types=rectypes,
                lig_types=ligtypes,
            ) 
            v.load_maps(maps_fn)
    else:
        v.load_maps(args.maps)
else:
    v = Vina(sf_name=args.scoring, cpu=args.vina_cpp_threads)
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

mol_supplier = Chem.SDMolSupplier(args.ligands, removeHs=False)
if args.name_from_prop:
    mol_supplier = MolSupplier(mol_supplier, name_from_prop=args.name_from_prop)

def wrap_mol_supplier(mol_supplier, *more_args):
    for mol in mol_supplier:
        yield (mol, *more_args)

# prepare ligand pdbqt
meeko_prep = MoleculePreparation(flexible_amides=args.flexible_amides)

if args.out_db is not None:
    if not _got_ringtail:
        raise ImportError from _ringtail_import_err
    rtc = RingtailCore(args.out_db)
    rtc.save_receptor(args.receptor)  # TODO JSON?
    rt_logger = logging.getLogger("ringtail")
    rt_logger.setLevel("WARNING")

def dock(args):
    mol, exhaustiveness = args
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
        if exhaustiveness is None:
            v.dock()
        else:
            v.dock(exhaustiveness=exhaustiveness)
        output_pdbqt = v.poses()
        vina_strings = {mol_name: output_pdbqt}
        return vina_strings
    except Exception as error:
        return error
    
if args.nr_process == 0:
    nr_cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(nr_cores - 1) # leave 1 for writing
    map_fn = pool.imap_unordered
elif args.nr_process > 1:
    pool = multiprocessing.Pool(args.nr_process - 1) # leave 1 for writing
    map_fn = pool.imap_unordered 
elif args.nr_process == 1:
    map_fn = map
else:
    print("--nr_process can't be negative")
    sys.exit(2)

if args.out_sdf is not None:
    w = Chem.SDWriter(args.out_sdf)

supplier = wrap_mol_supplier(mol_supplier, args.exhaustiveness)
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
            output_rdmol.SetProp("_Name", mol_name)
            w.write(output_rdmol)
    except Exception as error:
        root_logger.error(str(error))

if args.out_sdf is not None:
    w.close()
