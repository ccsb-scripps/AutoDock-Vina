#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Vina
#

import os
import glob
import stat
import sys

from .vina_wrapper import Vina as _Vina
from . import utils


class Vina:
    def __init__(self, sf_name='vina', exhaustiveness=8, max_evals=0, cpu=0, seed=0, verbosity=1):
        """Initialize a Vina object.

        Args:
            sf_name (str): Scoring function name to use (Vina or ad4) (default: vina)
            exhaustiveness (int): Number of MC run (default: 8)
            max_evals (int): Maximum number of evaluation (default: 0; use heuristics rules)
            cpu (int): Number of CPU to use (default: 0; use all of them)
            seed (int): Random seed (default: 0; ramdomly choosed)
            verbosity (int): verbosity 0: not output, 1: normal, 2: verbose (default: 1; some output)

        """
        sf_name = sf_name.lower()
        assert sf_name in ('vina', 'ad4'), 'Error: Scoring function %s not recognized. (only vina or ad4)' % sf_name

        self._vina = _Vina(sf_name, exhaustiveness, max_evals, cpu, seed, verbosity)
        
        self._sf_name = sf_name
        if sf_name == 'vina':
            self._weights = (-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 0.05846)
        else:
            self._weights = (0.1662, 0.1209, 0.1406, 0.1322, 50)
        self._rigid_receptor = None
        self._flex_receptor = None
        self._ligands = None
        self._center = None
        self._box_size = None
        self._spacing = None

    def __str__(self):
        """Print basic information about Vina object."""
        try:
            info = "Receptor (rigid): %s\n" % self._rigid_receptor
            info += "Receptor (flex): %s\n" % self._flex_receptor
            if isinstance(self._ligands, (list, tuple)):
                info += "Ligands: %s\n" % ", ".join(self._ligands)
            else:
                info += "Ligand: %s\n" % self._ligands
            info += "Scoring function: %s\n" % self._sf_name
            info += "Weights: %s\n" % " ".join(["%.6f" % i for i in self._weights])
            if self._center is not None:
                info += "Box center: %s\n" % " ".join(["%.3f" % i for i in self._center])
                info += "Box dimensions: %s\n" % " ".join(["%.2f" % i for i in self._box_size])
                info += "Box spacing: %.3f\n" % self._spacing
        except AttributeError:
            info = "Vina object is not defined."

        return info
    
    def cite(self):
        """Print citation message."""
        self._vina.cite()

    def set_receptor(self, rigid_pdbqt_filename, flex_pdbqt_filename=None):
        """Set receptor.

        Args:
            rigid_pdbqt_filename (str): rigid pdbqt receptor filename
            flex_pdbqt_filename (str): flexible residues pdbqt filename

        """
        # For the rigid part of the receptor
        assert os.path.exists(rigid_pdbqt_filename), 'Error: file %s does not exist.' % rigid_pdbqt_filename
        _, extension = os.path.splitext(rigid_pdbqt_filename)
        assert extension == '.pdbqt', 'Error: Vina requires a PDBQT file for the (rigid) receptor.'
        # For the flex part of the receptor
        if flex_pdbqt_filename is not None:
            assert os.path.exists(flex_pdbqt_filename), 'Error: file %s does not exist.' % flex_pdbqt_filename
            _, extension = os.path.splitext(flex_pdbqt_filename)
            assert extension == '.pdbqt', 'Error: Vina requires a PDBQT file for the (flex) receptor.'

        if flex_pdbqt_filename is None:
            self._vina.set_receptor(rigid_pdbqt_filename)
        else:
            self._vina.set_receptor(rigid_pdbqt_filename, flex_pdbqt_filename)

        self._rigid_receptor = rigid_pdbqt_filename
        self._flex_receptor = flex_pdbqt_filename

    def set_ligand(self, pdbqt_filename):
        """Set ligand(s).

        Args:
            pdbqt_filename (str or list): Name or list of PDBQT filename(s)

        """
        if not isinstance(pdbqt_filename, (list, tuple)):
            pdbqt_filename = [pdbqt_filename]

        for pf in pdbqt_filename:
            assert os.path.exists(pf), 'Error: file %s does not exist.' % pf
            _, extension = os.path.splitext(pf)
            assert extension == '.pdbqt', 'Error: Vina requires a PDBQT filename for the ligand.'

        if len(pdbqt_filename) == 1:
            self._vina.set_ligand(pdbqt_filename[0])
        else:
            self._vina.set_ligand(pdbqt_filename)

        self._ligands = pdbqt_filename

    def set_weights(self, weights):
        """Set potential weights for vina or ad4 scoring function.

        Args:
            weights (list): list or weights

        """
        assert isinstance(weights, (list, tuple)), 'Error: Cannot set weights (%s).' % weights

        if self._sf_name == 'vina':
            assert len(weights) == 7, 'Error: Number of weights does not correspond to Vina scoring function.' 
            self._vina.set_vina_weights(weights)
        else:
            assert len(weights) == 6, 'Error: Number of weights does not correspond to AD4 scoring function.'
            self._vina.set_ad4_weights(weights)

        self._weights = weights

    def compute_vina_maps(self, center, box_size, spacing=0.375):
        """Compute affinity maps using Vina scoring function.

        Args:
            center (list): center position
            box_siwe (list): size of the box in Angstrom
            spacing (float): grid spacing (default: 0.375)

        """
        assert len(center) == 3, 'Error: center of the box needs to be defined by (x, y, z) in Angstrom.'
        assert len(box_size) == 3, 'Error: box size needs to be defined by (a, b, c) in Angstrom.'
        assert all([i > 0 for i in box_size]), 'Error: box dimensions are required to be positive.'
        assert spacing > 0, 'Error: spacing should be positive.'

        x, y, z = center
        a, b, c = box_size
        self._vina.compute_vina_maps(x, y, z, a, b, c, spacing)
        self._center = center
        self._box_size = box_size
        self._spacing = spacing

    def load_maps(self, map_prefix_filename):
        """Load vina or ad4 affinity maps.

        Args:
            map_prefix_filename (str): affinity map prefix filename

        """
        existing_maps = glob.glob('%s.*.map' % map_prefix_filename)
        assert existing_maps, 'Error: Cannot find affinity maps with %s' % map_prefix_filename
        self._vina.load_maps(map_prefix_filename)
    
    def write_maps(self, map_prefix_filename='receptor', gpf_filename='NULL',
                   fld_filename='NULL', receptor_filename='NULL', overwrite=False):
        """Write affinity maps.

        Args:
            map_prefix_filename (str): affinity map pathname (path directory + prefix)
            gpf_filename (str): grid protein filename (default: NULL)
            fld_filename (str): fld filename (default: NULL)
            receptor filename (str): receptor filename (default: NULL)
            overwrite (bool): allow overwriting (default: false)

        """
        assert self._center is not None, 'Error: no affinity maps were defined.'
        if not overwrite:
            existing_maps = glob.glob('%s.*.map' % map_prefix_filename)
            assert not existing_maps, 'Error: Cannot overwrite existing affinity maps (%s)' % existing_maps
        self._vina.write_maps(map_prefix_filename, gpf_filename, fld_filename, receptor_filename)
    
    def write_pose(self, pdbqt_filename, remarks='', overwrite=False):
        """Write pose (after randomize or optimize).

        Args:
            pdbqt_filename (str): output PDBQT filename
            remarks (str): REMARKS to add in the PDBQT filename
            overwrite (bool): allow overwriting (default: false)

        """
        assert utils.check_file_writable(pdbqt_filename), 'Error: Cannot write pose at %s.' % pdbqt_filename
        if not overwrite:
            assert not os.path.exists(pdbqt_filename), 'Error: Cannot overwrite %s, already exists.' % pdbqt_filename
        self._vina.write_pose(pdbqt_filename, remarks)

    def write_docking_results(self, dlg_filename, n_poses=9, energy_range=3.0, overwrite=False):
        """Write poses from docking.

        Args:
            dlg_filename (str): docking ligand filename (PDBQT)
            n_pose (int): number of poses to write (default: 9)
            energy_range (float): maximum energy difference from best pose (default: 3.0 kcal/mol)
            overwrite (bool): allow overwriting (default: false)

        """
        assert utils.check_file_writable(dlg_filename), 'Error: Cannot write docking results at %s.' % dlg_filename
        if not overwrite:
            assert not os.path.exists(dlg_filename), 'Error: Cannot overwrite %s, already exists.' % dlg_filename
        assert n_poses > 0, 'Error: number of poses written must be positive.'
        assert energy_range > 0., 'Error: energy range must be positive.'
        self._vina.write_results(dlg_filename, n_poses, energy_range)

    def randomize(self):
        """Randomize the input ligand conformation."""
        self._vina.randomize()
    
    def score(self):
        """Score current pose.

        Returns:
            list: list of energies (total, lig_inter, flex_inter, other_inter, lig_intra, conf_independent)

        """
        return self._vina.score()

    def optimize(self):
        """Quick local BFGS optimization.

        Returns:
            list: list of energies (total, lig_inter, flex_inter, other_inter, lig_intra, conf_independent)

        """
        return self._vina.optimize()

    def dock(self, n_poses=20, min_rmsd=1.0):
        """Docking: global search optimization.

        Args:
            n_poses (int): number of pose to generate (default: 20)
            min_rmsd (float): minimum RMSD difference between poses (default: 1.0 Ansgtrom)

        """
        assert n_poses > 0, 'Error: number of poses to generate must be positive.'
        assert min_rmsd > 0., 'Error: minimal RMSD must be positive.'
        self._vina.global_search(n_poses, min_rmsd)
