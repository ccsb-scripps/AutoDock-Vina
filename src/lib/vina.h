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
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "conf.h"
#include "model.h"
#include "common.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "current_weights.h"
#include "quasi_newton.h"
#include "tee.h"
#include "coords.h" // add_to_output_container
#include "utils.h"


class Vina {
public:
    Vina(bool score_only = false, bool local_only = false, bool randomize_only = false, int cpu = 0, 
         int seed = 0, int exhaustiveness = 8, sz num_modes = 9, fl energy_range = 3.0, 
         bool no_cache = false, int verbosity = 2) {
        m_score_only = score_only;
        m_local_only = local_only;
        m_randomize_only = randomize_only;
        m_exhaustiveness = exhaustiveness;
        m_num_modes = num_modes;
        m_energy_range = energy_range;
        m_no_cache = no_cache;
        m_verbosity = verbosity;

        // Generate seed
        if (seed == 0) {
            m_seed = auto_seed();
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
        } else {
            m_cpu = cpu;
        }

        // Default VINA weigths
        fl weight_gauss1      = -0.035579;
        fl weight_gauss2      = -0.005156;
        fl weight_repulsion   =  0.840245;
        fl weight_hydrophobic = -0.035069;
        fl weight_hydrogen    = -0.587439;
        fl weight_rot         =  0.05846;
        m_weights.push_back(weight_gauss1);
        m_weights.push_back(weight_gauss2);
        m_weights.push_back(weight_repulsion);
        m_weights.push_back(weight_hydrophobic);
        m_weights.push_back(weight_hydrogen);
        m_weights.push_back(5 * weight_rot / 0.1 - 1);
    }

    void set_receptor(const std::string& rigid_name);
    void set_receptor(const std::string& rigid_name, const std::string& flex_name);
    void set_ligand(const std::string& ligand_name);
    void set_ligand(const std::vector<std::string>& ligand_name);
    void set_box(double center_x, double center_y, double center_z, int size_x, int size_y, int size_z, double granularity);
    void set_weights(flv weigths);
    void run(const std::string& out_name);

private:
    model m_m;
    grid_dims m_gd;
    flv m_weights;
    int m_cpu; 
    int m_seed; 
    int m_verbosity;
    int m_exhaustiveness;
    bool m_local_only; 
    bool m_score_only;
    bool m_randomize_only;
    bool m_no_cache;
    sz m_num_modes;
    fl m_energy_range;
    tee m_log;

    void write_all_output(model& m, const output_container& out, sz how_many, 
                          const std::string& output_name, const std::vector<std::string>& remarks);
    void do_randomization(model& m, const std::string& out_name, const vec& corner1, 
                          const vec& corner2, int seed, int verbosity, tee& log);
    void refine_structure(model& m, const precalculate& prec, non_cache& nc, output_type& out, 
                          const vec& cap, sz max_step);
    std::string vina_remark(fl e, fl lb, fl ub);
    output_container remove_redundant(const output_container& in, fl min_rmsd);
    void do_search(model& m, const scoring_function& sf, const precalculate& prec, 
                   const igrid& ig, const precalculate& prec_widened, const igrid& ig_widened, non_cache& nc, // nc.slope is changed
                   const std::string& out_name, const vec& corner1, const vec& corner2, const parallel_mc& par, 
                   fl energy_range, sz num_modes, int seed, int verbosity, bool score_only, bool local_only, 
                   tee& log, const terms& t, const flv& weights);
};


#endif