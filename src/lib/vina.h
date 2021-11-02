/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
		   The Olson Lab, 
		   The Scripps Research Institute

*/

#ifndef VINA_H
#define VINA_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include <boost/algorithm/string.hpp>
//#include <openbabel/mol.h>
#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "conf.h"
#include "model.h"
#include "common.h"
#include "cache.h"
#include "non_cache.h"
#include "ad4cache.h"
#include "quasi_newton.h"
#include "coords.h" // add_to_output_container
#include "utils.h"
#include "scoring_function.h"
#include "precalculate.h"


class Vina {
public:
	// Constructor
	Vina(const std::string &sf_name="vina", int cpu=0, int seed=0, int verbosity=1, bool no_refine=false, std::function<void(double)>* progress_callback = NULL) {
		m_verbosity = verbosity;
		m_receptor_initialized = false;
		m_ligand_initialized = false;
		m_map_initialized = false;
		m_seed = generate_seed(seed);
		m_no_refine = no_refine;
		m_progress_callback = progress_callback;

		// Look for the number of cpu
		if (cpu <= 0) {
			unsigned num_cpus = boost::thread::hardware_concurrency();

			if (num_cpus > 0) {
				m_cpu = num_cpus;
			} else {
				std::cerr << "WARNING: Could not determined the number of concurrent thread supported on this machine. ";
				std::cerr << "You might need to set it manually using cpu argument or fix the issue.\n";
				exit(EXIT_FAILURE);
			}
		} else {
			m_cpu = cpu;
		}

		if (sf_name.compare("vina") == 0) {
			m_sf_choice = SF_VINA;
			set_vina_weights();
		} else if (sf_name.compare("vinardo") == 0) {
			m_sf_choice = SF_VINARDO;
			set_vinardo_weights();
		} else if (sf_name.compare("ad4") == 0) {
			m_sf_choice = SF_AD42;
			set_ad4_weights();
		} else {
			std::cerr << "ERROR: Scoring function " << sf_name << " not implemented (choices: vina, vinardo or ad4)\n";
			exit (EXIT_FAILURE);
		}
	}
	// Destructor
	~Vina();

	void cite();
	int seed() { return m_seed; }
	void set_receptor(const std::string &rigid_name=std::string(), const std::string &flex_name=std::string());
	void set_ligand_from_string(const std::string &ligand_string);
	void set_ligand_from_string(const std::vector<std::string> &ligand_string);
	void set_ligand_from_file(const std::string& ligand_name);
	void set_ligand_from_file(const std::vector<std::string>& ligand_name);
	//void set_ligand(OpenBabel::OBMol* mol);
	//void set_ligand(std::vector<OpenBabel::OBMol*> mol);
	void set_vina_weights(double weight_gauss1=-0.035579, double weight_gauss2=-0.005156,
						       double weight_repulsion=0.840245, double weight_hydrophobic=-0.035069,
						       double weight_hydrogen=-0.587439, double weight_glue=50,
						       double weight_rot=0.05846);
	void set_vinardo_weights(double weight_gauss1=-0.045,
							       double weight_repulsion=0.8, double weight_hydrophobic=-0.035,
							       double weight_hydrogen=-0.600, double weight_glue=50,
							       double weight_rot=0.05846);
	void set_ad4_weights(double weight_ad4_vdw=0.1662, double weight_ad4_hb=0.1209,
						      double weight_ad4_elec=0.1406, double weight_ad4_dsolv=0.1322,
						      double weight_glue=50, double weight_ad4_rot=0.2983);
	std::vector<double> grid_dimensions_from_ligand(double buffer_size=4);
	void compute_vina_maps(double center_x, double center_y, double center_z,
								  double size_x, double size_y, double size_z,
								  double granularity=0.5, bool force_even_voxels=false);
	void load_maps(std::string maps);
	void randomize(const int max_steps=10000);
	std::vector<double> score();
	std::vector<double> optimize(const int max_steps=0);
	void global_search(const int exhaustiveness=8, const int n_poses=20, const double min_rmsd=1.0, const int max_evals=0);
	std::string get_poses(int how_many=9, double energy_range=3.0);
	std::vector< std::vector<double> > get_poses_coordinates(int how_many=9, double energy_range=3.0);
	std::vector< std::vector<double> > get_poses_energies(int how_many=9, double energy_range=3.0);
	void write_pose(const std::string &output_name, const std::string &remark = std::string());
	void write_poses(const std::string &output_name, int how_many=9, double energy_range=3.0);
	void write_maps(const std::string& map_prefix="receptor", const std::string& gpf_filename="NULL",
					    const std::string& fld_filename="NULL", const std::string& receptor_filename="NULL");
	void show_score(const std::vector<double> energies);

private:
	// model and poses
	model m_receptor;
	model m_model;
	output_container m_poses;
	//OpenBabel::OBMol m_mol;
	bool m_receptor_initialized;
	bool m_ligand_initialized;
	// scoring function
	scoring_function_choice m_sf_choice;
	flv m_weights;
	ScoringFunction m_scoring_function;
	precalculate_byatom m_precalculated_byatom;
	precalculate m_precalculated_sf;
	// maps
	cache m_grid;
	ad4cache m_ad4grid;
	non_cache m_non_cache;
	bool m_map_initialized;
	// global search
	int m_cpu;
	int m_seed;
	// others
	int m_verbosity;
	bool m_no_refine;
	std::function<void(double)>* m_progress_callback;

	std::string vina_remarks(output_type& pose, fl lb, fl ub);
	output_container remove_redundant(const output_container& in, fl min_rmsd);

	void set_forcefield();
	std::vector<double> score(double intramolecular_energy);
	std::vector<double> optimize(output_type& out, const int max_steps=0);
	int generate_seed(const int seed=0);
};

#endif
