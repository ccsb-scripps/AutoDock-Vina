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
#include "ad4cache.h"
#include "parse_error.h"
#include "quasi_newton.h"
#include "tee.h"
#include "coords.h" // add_to_output_container
#include "utils.h"
#include "scoring_function.h"
#include "precalculate.h"


class Vina {
public:
    // Constructor
    Vina(const std::string& sf_name="vina", int exhaustiveness=8, int cpu=0, int seed=0, bool no_cache=false, int verbosity=1) {
        m_exhaustiveness = exhaustiveness;
        m_verbosity = verbosity;
        m_receptor_initialized = false;
        m_ligand_initialized = false;
        m_map_initialized = false;
        m_seed = generate_seed(seed);

        if (m_exhaustiveness < 1) {
            std::cerr << "ERROR: Exhaustiveness must be 1 or greater";
            exit (EXIT_FAILURE);
        }

        // Look for the number of cpu
        if (cpu == 0) {
            unsigned num_cpus = boost::thread::hardware_concurrency();

            if (num_cpus > 0) {
                m_cpu = num_cpus;
            } else {
                m_cpu = 1;
            }
            
        } else if (cpu < 0) {
            m_cpu = 1;
            std::cerr << "WARNING: Number of CPUs set to a value lower than 0, it was automatically set to 1 per default.\n";
        } else {
            m_cpu = cpu;
        }

        if (exhaustiveness < cpu) {
            std::cerr << "WARNING: At low exhaustiveness, it may be impossible to utilize all CPUs.\n";
        }

        if (sf_name.compare("vina") == 0) {
            m_sf_choice = SF_VINA;
            set_vina_weights();
        } else if (sf_name.compare("ad4") == 0) {
            m_sf_choice = SF_AD42;
            set_ad4_weights();
        } else {
            std::cerr << "ERROR: Scoring function " << sf_name << " not implemented (choices: vina or ad4)\n";
            exit (EXIT_FAILURE);
        }
    }
    // Destructor
    ~Vina();

    void cite();
    void set_receptor(const std::string& rigid_name);
    void set_receptor(const std::string& rigid_name, const std::string& flex_name);
    void set_ligand(const std::string& ligand_name);
    void set_ligand(const std::vector<std::string>& ligand_name);
    //void set_ligand(OpenBabel::OBMol* mol);
    //void set_ligand(std::vector<OpenBabel::OBMol*> mol);
    void set_vina_weights(double weight_gauss1=-0.035579,  double weight_gauss2=-0.005156, 
                          double weight_repulsion=0.840245, double weight_hydrophobic=-0.035069, 
                          double weight_hydrogen=-0.587439, double weight_glue=50,
                          double weight_rot=0.05846);
    void set_ad4_weights(double weight_ad4_vdw =0.1662, double weight_ad4_hb=0.1209, 
                         double weight_ad4_elec=0.1406, double weight_ad4_dsolv=0.1322, 
                         double weight_glue=50, double weight_ad4_rot =0.2983);
    void compute_vina_maps(double center_x, double center_y, double center_z, double size_x, double size_y, double size_z, double granularity=0.375);
    void load_ad4_maps(std::string ad4_maps);
    void randomize(const int max_steps=10000);
    std::vector<double> score();
    std::vector<double> optimize(const int max_steps=0);
    void global_search(const int n_poses=20, const double min_rmsd=1.0);
    void write_results(const std::string& output_name, int how_many=9, double energy_range=3.0);
    void write_pose(const std::string& output_name, const std::string& remark=std::string());
    void write_maps(const std::string& map_prefix="receptor", const std::string& gpf_filename="NULL",
                    const std::string& fld_filename="NULL", const std::string& receptor_filename="NULL");

private:
    //OpenBabel::OBMol m_mol;
    // model and poses
    model m_receptor;
    model m_model;
    output_container m_poses;
    bool m_receptor_initialized;
    bool m_ligand_initialized;
    // scoring function
    scoring_function_choice m_sf_choice;
    flv m_weights;
    ScoringFunction m_scoring_function;
    precalculate_byatom m_precalculated_byatom;
    // maps
    cache m_grid;
    ad4cache m_ad4grid;
    grid_dims m_gd;
    vec m_corner1;
    vec m_corner2;
    bool m_map_initialized;
    // global search
    int m_cpu;
    int m_seed;
    int m_exhaustiveness;
    // others
    int m_verbosity;

    std::string vina_remark(fl e, fl lb, fl ub);
    std::string vina_remarks(output_type& pose, fl lb, fl ub);
    output_container remove_redundant(const output_container& in, fl min_rmsd);

    void set_forcefield();
    void set_vina_box(double center_x, double center_y, double center_z, double size_x, double size_y, double size_z, double granularity=0.375);
    std::vector<double> score(double intramolecular_energy);
    std::vector<double> optimize(output_type& out, const int max_steps=0);
    int generate_seed(const int seed=0);
    void show_score(const std::vector<double> energies);
};

#endif
