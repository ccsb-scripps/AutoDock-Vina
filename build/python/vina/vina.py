#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Vina
#

import os
import glob
import stat
import sys

import numpy as np

from .vina_wrapper import Vina as _Vina
from . import utils


class Vina:
    def __init__(self, sf_name='vina', cpu=0, seed=0, no_refine=False, verbosity=1):
        """Initialize a Vina object.

        Args:
            sf_name (str): Scoring function name to use (Vina or ad4) (default: vina)
            cpu (int): Number of CPU to use (default: 0; use all of them)
            seed (int): Random seed (default: 0; ramdomly choosed)
            no_refine (boolean): when receptor is provided, do not use explicit receptor atoms
                (instead of precalculated grids) for: (1) local optimization and scoring after docking,
                (2) --local_only jobs, and (3) --score_only jobs (default: False)
            verbosity (int): verbosity 0: not output, 1: normal, 2: verbose (default: 1; some output)

        """
        sf_name = sf_name.lower()
        if not sf_name in ('vina', 'vinardo', 'ad4'):
            raise ValueError('Error: Scoring function %s not recognized. (only vina, vinardo or ad4)' % sf_name)

        self._vina = _Vina(sf_name, cpu, seed, verbosity, no_refine)
        
        self._sf_name = sf_name
        if sf_name == 'vina':
            self._weights = (-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 0.05846)
        elif sf_name == 'vinardo':
            self._weights = (-0.045, 0.8, -0.035, -0.6, 50, 0.05846)
        else:
            self._weights = (0.1662, 0.1209, 0.1406, 0.1322, 50)
        self._rigid_receptor = None
        self._flex_receptor = None
        self._ligands = None
        self._center = None
        self._box_size = None
        self._spacing = None
        self._even_nelements = True
        if sf_name in ('vina', 'vinardo'):
            self._no_refine = no_refine
        else:
            self._no_refine = False
        self._seed = self._vina.seed()

    def __str__(self):
        """Print basic information about Docking configuration (rigid receptor, flex receptor, 
            ligands, scoring function, weights, no_refine, box center, box dimensions, 
            box spacing, box even nelements, seed).
            
        """
        info = "Receptor (rigid): %s\n" % self._rigid_receptor
        info += "Receptor (flex): %s\n" % self._flex_receptor

        if isinstance(self._ligands, (list, tuple)):
            info += "Ligands: %s\n" % ", ".join(self._ligands)
        else:
            info += "Ligand: %s\n" % self._ligands

        info += "Scoring function: %s\n" % self._sf_name
        info += "Weights: %s\n" % " ".join(["%.6f" % i for i in self._weights])
        info += "No refinement with receptor atoms: %s\n" % self._no_refine

        if self._center is not None:
            info += "Box center: %s\n" % " ".join(["%.3f" % i for i in self._center])
            info += "Box dimensions: %s\n" % " ".join(["%.2f" % i for i in self._box_size])
            info += "Box spacing: %.3f\n" % self._spacing
            info += "Box even NELEMENTS: %s\n" % self._even_nelements

        info += "Seed: %s" % self._seed

        return info
    
    def cite(self):
        """Print citation message."""
        self._vina.cite()

    def info(self):
        """Return information about the docking configuration.
        
        Returns:
            dict (str): Dictionary of information about the docking configuration.
            
            Information:
                rigid_receptor (str),
                flex_receptor (str),
                ligands (list),
                scoring_function (str),
                weights (tuple),
                no_refine (bool),
                box_center (list),
                box_size (list),
                box_spacing (float),
                box_even_elements (bool),
                seed (int)
        
        """
        info = {}

        info['rigid_receptor'] = self._rigid_receptor
        info['flex_receptor'] = self._flex_receptor
        info['ligands'] = self._ligands
        info['scoring_function'] = self._sf_name
        info['weights'] = self._weights
        info['no_refine'] = self._no_refine
        info['box_center'] = self._center
        info['box_size'] = self._box_size
        info['box_spacing'] = self._spacing
        info['box_even_elements'] = self._even_nelements
        info['seed'] = self._seed
        
        return info

    def set_receptor(self, rigid_pdbqt_filename=None, flex_pdbqt_filename=None):
        """Set receptor.

        Args:
            rigid_pdbqt_filename (str): rigid pdbqt receptor filename (default: None)
            flex_pdbqt_filename (str): flexible residues pdbqt filename (default: None)

        """
        if rigid_pdbqt_filename is None and flex_pdbqt_filename is None:
            raise ValueError('Error: No (rigid) receptor or flexible residues were specified')

        # For the rigid part of the receptor
        if rigid_pdbqt_filename is not None:
            if not os.path.exists(rigid_pdbqt_filename):
                raise RuntimeError('Error: file %s does not exist.' % rigid_pdbqt_filename)
            _, extension = os.path.splitext(rigid_pdbqt_filename)
            if not extension == '.pdbqt':
                raise TypeError('Error: Vina requires a PDBQT file for the (rigid) receptor.')

        # For the flex part of the receptor
        if flex_pdbqt_filename is not None:
            if not os.path.exists(flex_pdbqt_filename):
                raise RuntimeError('Error: file %s does not exist.' % flex_pdbqt_filename)
            _, extension = os.path.splitext(flex_pdbqt_filename)
            if not extension == '.pdbqt':
                raise TypeError('Error: Vina requires a PDBQT file for the (flex) receptor.')

        if rigid_pdbqt_filename is not None:
            if flex_pdbqt_filename is not None:
                self._vina.set_receptor(rigid_pdbqt_filename, flex_pdbqt_filename)
            else:
                self._vina.set_receptor(rigid_pdbqt_filename)
        else:
            self._vina.set_receptor('', flex_pdbqt_filename)

        self._rigid_receptor = rigid_pdbqt_filename
        self._flex_receptor = flex_pdbqt_filename

    def set_ligand_from_file(self, pdbqt_filename):
        """Set ligand(s) from a file. The chemical file format must be PDBQT.

        Args:
            pdbqt_filename (str or list): Name or list of PDBQT filename(s)

        """
        if not isinstance(pdbqt_filename, (list, tuple)):
            pdbqt_filename = [pdbqt_filename]

        for pf in pdbqt_filename:
            if not os.path.exists(pf):
                raise RuntimeError('Error: file %s does not exist.' % pf)
            _, extension = os.path.splitext(pf)
            if not extension == '.pdbqt':
                raise TypeError('Error: Vina requires a PDBQT file for the ligand.')

        if len(pdbqt_filename) == 1:
            self._vina.set_ligand_from_file(pdbqt_filename[0])
        else:
            self._vina.set_ligand_from_file(pdbqt_filename)

        self._ligands = pdbqt_filename
    
    def set_ligand_from_string(self, pdbqt_string):
        """Set ligand(s) from a string. The chemical file format must be PDBQT.

        Args:
            pdbqt_string (str or list): string or list of PDBQT strings

        """
        if not isinstance(pdbqt_string, (list, tuple)):
            pdbqt_string = [pdbqt_string]
        
        for ps in pdbqt_string:
            if not isinstance(ps, str):
                raise TypeError('Error: %s is not a string.' % ps)

        if len(pdbqt_string) == 1:
            self._vina.set_ligand_from_string(pdbqt_string[0])
        else:
            self._vina.set_ligand_from_string(pdbqt_string)

        self._ligands = pdbqt_string

    def set_weights(self, weights):
        """Set potential weights for vina or ad4 scoring function.

        Args:
            weights (list): list or weights

        """
        if not isinstance(weights, (list, tuple)):
            raise TypeError('Error: Cannot set weights (%s).' % weights)
        if self._sf_name == 'vina':
            if len(weights) != 7:
                raise ValueError('Error: Number of weights does not correspond to Vina scoring function.' )
            self._vina.set_vina_weights(*weights)
        else:
            if len(weights) != 6:
                raise ValueError('Error: Number of weights does not correspond to AD4 or Vinardo scoring function.')
                if self._sf_name == 'ad4':
                    self._vina.set_ad4_weights(*weights)
                else:
                    self._vina.set_vinardo_weights(*weights)

        self._weights = weights

    def compute_vina_maps(self, center, box_size, spacing=0.375, force_even_voxels=False):
        """Compute affinity maps using Vina scoring function.

        Args:
            center (list): center position
            box_size (list): size of the box in Angstrom
            spacing (float): grid spacing (default: 0.375)
            force_even_voxels (boolean): Force the number of voxels (NPTS/NELEMENTS) to be an even number
                (and forcing the number of grid points to be odd) (default: False)

        """
        if len(center) != 3:
            raise ValueError('Error: center of the box needs to be defined by (x, y, z) in Angstrom.')
        elif len(box_size) != 3:
            raise ValueError('Error: box size needs to be defined by (a, b, c) in Angstrom.')
        elif not all([i > 0 for i in box_size]):
            raise ValueError('Error: box dimensions are required to be positive.')
        elif spacing <= 0:
            raise ValueError('Error: spacing should be greater than zero.')

        x, y, z = center
        a, b, c = box_size

        self._vina.compute_vina_maps(x, y, z, a, b, c, spacing, force_even_voxels)

        self._center = center
        self._box_size = box_size
        self._spacing = spacing
        self._voxels = np.ceil(np.array(box_size) / self._spacing).astype(np.int)

        # Necessary step to know if we can write maps or not later
        if force_even_voxels:
            self._voxels[0] += int(self._voxels[0] % 2 == 1)
            self._voxels[1] += int(self._voxels[1] % 2 == 1)
            self._voxels[2] += int(self._voxels[2] % 2 == 1)

    def load_maps(self, map_prefix_filename):
        """Load vina or ad4 affinity maps.

        Args:
            map_prefix_filename (str): affinity map prefix filename

        """
        if not glob.glob('%s.*.map' % map_prefix_filename):
            raise RuntimeError('Error: Cannot find affinity maps with %s' % map_prefix_filename)

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
        if self._center is None:
            raise RuntimeError('Error: no affinity maps were defined.')
        elif not overwrite:
            existing_maps = glob.glob('%s.*.map' % map_prefix_filename)
            if existing_maps:
                raise RuntimeError('Error: Cannot overwrite existing affinity maps (%s)' % existing_maps)

        if any(self._voxels % 2 == 1):
            error_msg = 'Error: Can\'t write maps. Number of voxels (NPTS/NELEMENTS) is odd: %s.' % self._voxels
            error_msg += ' Set \'force_even_voxels\' to True when computing vina maps.'
            raise RuntimeError(error_msg)

        self._vina.write_maps(map_prefix_filename, gpf_filename, fld_filename, receptor_filename)
    
    def write_pose(self, pdbqt_filename, remarks='', overwrite=False):
        """Write pose (after randomize or optimize).

        Args:
            pdbqt_filename (str): output PDBQT filename
            remarks (str): REMARKS to add in the PDBQT filename
            overwrite (bool): allow overwriting (default: false)

        """
        if not utils.check_file_writable(pdbqt_filename):
            raise RuntimeError('Error: Cannot write pose at %s.' % pdbqt_filename)
        elif not overwrite:
            if os.path.exists(pdbqt_filename):
                raise RuntimeError('Error: Cannot overwrite %s, already exists.' % pdbqt_filename)

        self._vina.write_pose(pdbqt_filename, remarks)
    
    def write_poses(self, pdbqt_filename, n_poses=9, energy_range=3.0, overwrite=False):
        """Write poses from docking.

        Args:
            pdbqt_filename (str): PDBQT file containing poses found
            n_pose (int): number of poses to write (default: 9)
            energy_range (float): maximum energy difference from best pose (default: 3.0 kcal/mol)
            overwrite (bool): allow overwriting (default: false)

        """
        if not utils.check_file_writable(pdbqt_filename):
            raise RuntimeError(
                'Error: Cannot write docking results at %s.' % pdbqt_filename)
        elif not overwrite:
            if os.path.exists(pdbqt_filename):
                raise RuntimeError(
                    'Error: Cannot overwrite %s, already exists.' % pdbqt_filename)
        elif n_poses <= 0:
            raise ValueError(
                'Error: number of poses written must be greater than zero.')
        elif energy_range <= 0:
            raise ValueError('Error: energy range must be greater than zero.')

        self._vina.write_poses(pdbqt_filename, n_poses, energy_range)
    
    def poses(self, n_poses=9, energy_range=3.0, coordinates_only=False):
        """Get poses from docking.

        Args:
            n_pose (int): number of poses to retrieve (default: 9)
            energy_range (float): maximum energy difference from best pose (default: 3.0 kcal/mol)
            coordinates_only (bool): return coordinates for each pose only

        Returns:
            str/ndarray: PDBQT file string or Array of coordinates of poses if coordinates_only=True

        """
        if n_poses <= 0:
            raise ValueError('Error: number of poses written must be greater than zero.')
        elif energy_range <= 0:
            raise ValueError('Error: energy range must be greater than zero.')
        
        if coordinates_only:
            # Dirty hack to get the coordinates from C++, it is not advised to have vector<vector<vector<double>>>
            coordinates = np.array(self._vina.get_poses_coordinates(n_poses, energy_range))
            coordinates = coordinates.reshape((coordinates.shape[0], coordinates.shape[1] // 3, 3))
            coordinates = np.around(coordinates, decimals=3)
            return coordinates
        else:
            return self._vina.get_poses(n_poses, energy_range)
    
    def energies(self, n_poses=9, energy_range=3.0):
        """Get pose energies from docking.

        Args:
            n_pose (int): number of poses to retrieve (default: 9)
            energy_range (float): maximum energy difference from best pose (default: 3.0 kcal/mol)
        
        Returns:
            ndarray: Array of energies from each pose (rows=poses, columns=energies) 
            
            Vina/Vinardo FF:
                columns=[total, inter, intra, torsions, intra best pose]
            
            AutoDock FF:
                columns=[total, inter, intra, torsions, -intra]

        """
        if n_poses <= 0:
            raise ValueError('Error: number of poses written must be greater than zero.')
        elif energy_range <= 0:
            raise ValueError('Error: energy range must be greater than zero.')

        return np.around(self._vina.get_poses_energies(n_poses, energy_range), decimals=3)

    def randomize(self):
        """Randomize the input ligand conformation."""
        self._vina.randomize()
    
    def score(self):
        """Score current pose.

        Returns:
            nadarray: Array of energies from current pose.

            Vina/Vinardo FF:
                columns=[total, lig_inter, flex_inter, other_inter, flex_intra, lig_intra, torsions, lig_intra best pose]
            
            AutoDock FF:
                columns=[total, lig_inter, flex_inter, other_inter, flex_intra, lig_intra, torsions, -lig_intra]

        """
        # It does not make sense to report energies with a precision higher than 3
        # since the coordinates precision is 3.
        energies = np.around(self._vina.score(), decimals=3)
        return energies

    def optimize(self, max_steps=0):
        """Quick local BFGS energy optimization.

        Args:
            max_steps (int): Maximum number of local minimization steps (default: 0). When max_steps is equal to 0,
                             the maximum number of steps will be equal to (25 + num_movable_atoms) / 3).

        Returns:
            nadarray: Array of energies from optimized pose.

            Vina/Vinardo FF:
                columns=[total, lig_inter, flex_inter, other_inter, flex_intra, lig_intra, torsions, lig_intra best pose]
            
            AutoDock FF:
                columns=[total, lig_inter, flex_inter, other_inter, flex_intra, lig_intra, torsions, -lig_intra]

        """
        if max_steps < 0:
            raise ValueError('Error: max_steps cannot be negative.')

        # It does not make sense to report energies with a precision higher than 3
        # since the coordinates precision is 3.
        energies = np.around(self._vina.optimize(max_steps), decimals=3)
        return energies

    def dock(self, exhaustiveness=8, n_poses=20, min_rmsd=1.0, max_evals=0):
        """Docking: global search optimization.

        Args:
            exhaustiveness (int): Number of MC run (default: 8)
            n_poses (int): number of pose to generate (default: 20)
            min_rmsd (float): minimum RMSD difference between poses (default: 1.0 Ansgtrom)
            max_evals (int): Maximum number of evaluation (default: 0; use heuristics rules)

        """
        if exhaustiveness <= 0:
            raise ValueError('Error: exhaustiveness must be 1 or greater.')
        elif n_poses <= 0:
            raise ValueError('Error: number of poses to generate must be greater than zero.')
        elif min_rmsd <= 0:
            raise ValueError('Error: minimal RMSD must be greater than zero.')
        elif max_evals < 0:
            raise ValueError('Error: maximum evaluations must be positive.')

        self._vina.global_search(exhaustiveness, n_poses, min_rmsd, max_evals)
