#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Vina CLI
#

import argparse
import os

import numpy as np

from .vina import Vina


def cmd_lineparser():
    parser = argparse.ArgumentParser(description='AutoDock-Vina 1.2.0 (Python CLI)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Input receptor, flex and ligand
    parser.add_argument('-r', '--receptor', dest='receptor', default=None,
                        type=str, action='store', help='rigid part of the receptor (PDBQT)')
    parser.add_argument('-f', '--flex', dest='flex', default=None,
                        type=str, action='store', help='flexible side chains, if any (PDBQT)')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-l', '--ligand', dest='ligands', default=None,
                       nargs='+', action='store', help='ligand (PDBQT)')
    group.add_argument('-b', '--batch', dest='batch', default=None,
                       nargs='+', action='store', help='ligand (PDBQT)')
    # Scoring function
    parser.add_argument('-s', '--scoring', dest='sf_name', default='vina',
                        type=str, choices=['ad4', 'vina', 'vinardo'], action='store',
                        help='scoring function (ad4, vina or vinardo)')
    # Dimensions maps
    group_dim = parser.add_argument_group('Search space options')
    group_dim.add_argument('-m', '--maps', dest='maps', default=None,
                        action='store', help='affinity maps for ad4, vina or vinardo scoring function')
    group_dim.add_argument('--center_x', dest='center_x', default=None,
                        type=float, action='store', help='X coordinate of the center (Angstrom)')
    group_dim.add_argument('--center_y', dest='center_y', default=None,
                        type=float, action='store', help='Y coordinate of the center (Angstrom)')
    group_dim.add_argument('--center_z', dest='center_z', default=None,
                        type=float, action='store', help='Z coordinate of the center (Angstrom)')
    group_dim.add_argument('--size_x', dest='size_x', default=None,
                        type=float, action='store', help='size in the X dimension (Angstrom)')
    group_dim.add_argument('--size_y', dest='size_y', default=None,
                        type=float, action='store', help='size in the Y dimension (Angstrom)')
    group_dim.add_argument('--size_z', dest='size_z', default=None,
                        type=float, action='store', help='size in the Z dimension (Angstrom)')
    # Actions
    group_actions = parser.add_argument_group('Action options')
    group_actions.add_argument('--randomize_only', dest='randomize_only', default=False,
                        action='store_true', help='randomize input, attempting to avoid clashes')
    group_actions.add_argument('--score_only', dest='score_only', default=False,
                        action='store_true', help='score only - needs maps')
    group_actions.add_argument('--local_only', dest='local_only', default=False,
                        action='store_true', help='do local search only')
    # Output
    group_out = parser.add_argument_group('Output options')
    group_out.add_argument('-o', '--out', dest='out', default=None,
                        type=str, action='store',
                        help='output models (PDBQT), the default is chosen based on the ligand file name')
    group_out.add_argument('-d', '--dir', dest='dir', default=None,
                        type=str, action='store',
                        help='output directory for batch mode')
    group_out.add_argument('--write_maps', dest='write_maps', default=None,
                        type=str, action='store',
                        help='output filename (directory + prefix name) for maps')
    # Extra arguments
    group_extra = parser.add_argument_group('Extra options')
    group_extra.add_argument('--cpu', dest='cpu', default=0,
                        type=int, action='store',
                        help='the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)')
    group_extra.add_argument('--seed', dest='seed', default=0,
                        type=int, action='store',
                        help='explicit random seed')
    group_extra.add_argument('--exhaustiveness', dest='exhaustiveness', default=8,
                        type=int, action='store',
                        help='exhaustiveness of the global search (roughly proportional to time): 1+')
    group_extra.add_argument('--max_evals', dest='max_evals', default=0,
                        type=int, action='store',
                        help='number of evaluations in each MC run (if zero, which is the default, the number of MC steps is based on heuristics)')
    group_extra.add_argument('--num_modes', dest='num_modes', default=9,
                        type=int, action='store',
                        help='maximum number of binding modes to generate')
    group_extra.add_argument('--min_rmsd', dest='min_rmsd', default=1.0,
                        type=float, action='store',
                        help='minimum RMSD between output poses')
    group_extra.add_argument('--energy_range', dest='energy_range', default=3.0,
                        type=float, action='store',
                        help='maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)')
    group_extra.add_argument('--spacing', dest='grid_spacing', default=0.375,
                        type=float, action='store',
                        help='grid spacing (Angstrom)')
    group_extra.add_argument('--verbosity', dest='verbosity', default=1,
                        type=int, action='store',
                        help='verbosity (0=no output, 1=normal, 2=verbose)')
    # Scoring function weight
    group_weight = parser.add_argument_group('Scoring function weights options')
    # Vina
    group_weight.add_argument('--weight_gauss1', dest='weight_gauss1', default=-0.035579,
                             type=float, action='store', help='gauss_1 weight')
    group_weight.add_argument('--weight_gauss2', dest='weight_gauss2', default=-0.005156,
                             type=float, action='store', help='gauss_2 weight')
    group_weight.add_argument('--weight_repulsion', dest='weight_repulsion', default=0.840245,
                             type=float, action='store', help='repulsion weight')
    group_weight.add_argument('--weight_hydrophobic', dest='weight_hydrophobic', default=-0.035069,
                             type=float, action='store', help='hydrophobic weight')
    group_weight.add_argument('--weight_hydrogen', dest='weight_hydrogen', default=-0.587439,
                             type=float, action='store', help='Hydrogen bond weights')
    group_weight.add_argument('--weight_rot', dest='weight_rot', default=0.05846,
                             type=float, action='store', help='N_rot weight')
    # Vinardo
    group_weight.add_argument('--weight_vinardo_gauss1', dest='weight_vinardo_gauss1', default=-0.045,
                              type=float, action='store', help='vinardo gauss_1 weight')
    group_weight.add_argument('--weight_vinardo_repulsion', dest='weight_vinardo_repulsion', default=0.8,
                              type=float, action='store', help='vinardo repulsion weight')
    group_weight.add_argument('--weight_vinardo_hydrophobic', dest='weight_vinardo_hydrophobic', default=-0.035,
                              type=float, action='store', help='vinardo hydrophobic weight')
    group_weight.add_argument('--weight_vinardo_hydrogen', dest='weight_vinardo_hydrogen', default=-0.600,
                              type=float, action='store', help='vinardo Hydrogen bond weights')
    group_weight.add_argument('--weight_vinardo_rot', dest='weight_vinardo_rot', default=0.05846,
                              type=float, action='store', help='vinardo N_rot weight')
    # AutoDock 4
    group_weight.add_argument('--weight_ad4_vdw', dest='weight_ad4_vdw', default=0.1662,
                             type=float, action='store', help='ad4_vdw weight')
    group_weight.add_argument('--weight_ad4_hb', dest='weight_ad4_hb', default=0.1209,
                             type=float, action='store', help='ad4_hb weight')
    group_weight.add_argument('--weight_ad4_elec', dest='weight_ad4_elec', default=0.1406,
                             type=float, action='store', help='ad4_electstrostatic weight')
    group_weight.add_argument('--weight_ad4_desolv', dest='weight_ad4_desolv', default=0.1322,
                             type=float, action='store', help='ad4_desolvation weight')
    group_weight.add_argument('--weight_ad4_rot', dest='weight_ad4_rot', default=0.2983,
                             type=float, action='store', help='ad4_rot weight')
    # Macrocycle
    group_weight.add_argument('--weight_glue', dest='weight_glue', default=50.0,
                             type=float, action='store', help='macrocycle glue weight')

    args = parser.parse_args()
    
    # We need the rigid receptor or the maps, not both
    if args.receptor is not None and args.maps is not None:
        parser.error('ERROR: Cannot specify both receptor and affinity maps at the same time, --flex argument is allowed with receptor or maps.')
    
    if args.maps is None:
        if not all(v is not None for v in [args.center_x, args.center_y, args.center_z]):
            parser.error('ERROR: The center of the grid box was not defined correctly (X %s Y %s Z %s).' \
                         % (args.center_x, args.center_y, args.center_z))
        if not all(v is not None for v in [args.size_x, args.size_y, args.size_z]):
            parser.error('ERROR: The size of the grid box was not defined correctly (X %s Y %s Z %s).' \
                         % (args.size_x, args.size_y, args.size_z))

    # Check that arguments are compatible with the scoring function
    if args.sf_name == 'vina' or args.sf_name == 'vinardo':
        if args.receptor is None and args.maps is None:
            parser.error('ERROR: The receptor or affinity maps must be specified.')
    elif args.sf_name == 'ad4':
        if args.receptor is not None:
            parser.error('ERROR: No receptor allowed, only --flex argument with the AD4 scoring function.')
        if args.maps is None:
            parser.error('ERROR: Affinity maps are missing.')
    
    # Check that arguments are compatible with the ligand/batch mode
    if args.batch is None and args.ligands is None:
        parser.error('ERROR: Missing ligand(s).')

    if args.ligands is not None:
        if args.dir is not None:
            parser.error('ERROR: In --ligand mode, --dir argument cannot be specified.')
        if args.out is None and not args.score_only:
            if len(args.ligands) == 1:
                args.out = '%s_out.pdbqt' % os.path.splitext(args.ligands[0])[0]
            else:
                parser.error('ERROR: Output name must be defined when docking simultaneously multiple ligands.')
    elif args.batch is not None:
        if args.dir is None:
            parser.error('ERROR: Need to specify an output directory for batch mode.')
        else:
            if not os.path.isdir(args.dir):
                parser.error('ERROR: Directory %s does not exist.' % args.args)

    return args


def main():
    args = cmd_lineparser()

    if args.verbosity > 0:
        print('Scoring function : %s' % args.sf_name)
        if args.receptor is not None:
            print('Rigid receptor: %s' % args.receptor)
        if args.flex is not None:
            print('Flex receptor: %s ' % flex_name)

        if len(args.ligands) == 1:
                print('Ligand: %s' % args.ligands[0])
        elif len(args.ligands) > 1:
            print('Ligands:')
            for ligand in args.ligands:
                print('  - %s' % ligand)
        elif len(args.batch) > 1:
                print('Ligands (batch mode): %d molecules.' % len(args.batch))

        if args.maps is None:
            print('Center: X %.3f Y %.3f Z %.3f' % (args.center_x, args.center_y, args.center_z))
            print('Size: X %4d Y %4d Z %4d' % (args.size_x, args.size_y, args.size_z))
            print('Grid space: %f' % args.grid_spacing)

        print('Exhaustiveness: %d' % args.exhaustiveness)
        print('CPU: %d' % args.cpu)
        if args.seed != 0:
            print('Seed: %d' % args.seed)
        print('Verbosity: %d' % args.verbosity)

    v = Vina(args.sf_name, args.cpu, args.seed, args.verbosity)

    # rigid_name variable can be ignored for AD4
    if args.receptor is not None or args.flex is not None:
        v.set_receptor(args.receptor, args.flex)

    # Technically we don't have to initialize weights,
    # because they are initialized during the Vina object creation with the default weights
    # but we still do it in case the user decided to change them
    if args.sf_name == 'vina':
        v.set_weights([args.weight_gauss1, args.weight_gauss2, args.weight_repulsion,
                       args.weight_hydrophobic, args.weight_hydrogen, args.weight_glue,
                       args.weight_rot])
    elif args.sf_name == 'vinardo':
        v.set_weights([args.weight_vinardo_gauss1, args.weight_vinardo_repulsion,
                       args.weight_vinardo_hydrophobic, args.weight_vinardo_hydrogen, args.weight_vinardo_glue,
                       args.weight_vinardo_rot])
    else:
        v.set_weights([args.weight_ad4_vdw, args.weight_ad4_hb, args.weight_ad4_elec,
                       args.weight_ad4_dsolv, args.weight_glue, args.weight_ad4_rot])
        v.load_maps(args.maps)

        # It works, but why would you do this?!
        if args.write_maps is not None:
            v.write_maps(args.write_maps)

    if args.ligands is not None:
        v.set_ligand_from_file(args.ligands)

        if args.sf_name == 'vina' or args.sf_name == 'vinardo':
            if args.maps is not None:
                v.load_maps(args.maps);
            else:
                # Will compute maps only for Vina atom types in the ligand(s)
                box_center = (args.center_x, args.center_y, args.center_z)
                box_size = (args.size_x, args.size_y, args.size_z)
                v.compute_vina_maps(box_center, box_size, args.grid_spacing)

                if args.write_maps is not None:
                    v.write_maps(args.write_maps)

        if args.randomize_only:
            v.randomize()
            v.write_pose(args.out)
        elif args.score_only:
            v.score()
        elif args.local_only:
            v.optimize()
            v.write_pose(args.name)
        else:
            v.dock(args.exhaustiveness, args.num_modes, args.min_rmsd, args.max_evals)
            v.write_poses(args.out, args.num_modes, args.energy_range)
    else:
        if args.sf_name == 'vina' or args.sf_name == 'vinardo':
            if args.maps is not None:
                v.load_maps(args.maps);
            else:
                # Will compute maps for all the Vina atom types
                box_center = (args.center_x, args.center_y, args.center_z)
                box_size = (args.size_x, args.size_y, args.size_z)
                v.compute_vina_maps(box_center, box_size, args.grid_spacing)

                if args.write_maps is not None:
                    v.write_maps(args.write_maps)

        for ligand in args.batch:
            v.set_ligand_from_file(ligand)

            molecule_name = os.path.splitext(os.path.basename(ligand))[0]
            out_name = '%s%s%s_out.pdbqt' % (args.dir, os.path.sep(), molecule_name)

            if args.randomize_only:
                v.randomize()
                v.write_pose(out_name)
            elif args.score_only:
                v.score()
            elif args.local_only:
                v.optimize()
                v.write_pose(out_name)
            else:
                v.dock(args.exhaustiveness, args.num_modes, args.min_rmsd, args.max_evals)
                v.write_poses(args.out_name, args.num_modes, args.energy_range)


if __name__ == '__main__':
    main()
