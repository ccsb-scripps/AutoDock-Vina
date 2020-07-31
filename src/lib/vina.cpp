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

#include "vina.h"


int Vina::generate_seed(const int seed) {
    // Seed generator, if the global seed (m_seed) was defined to 0
    // it seems that we want to generate random seed, otherwise it means
    // that we want to use a particular seed.
    if (seed == 0) {
        return auto_seed();
    } else {
        return seed;
    }
}

std::string Vina::vina_remark(fl e, fl lb, fl ub) {
    std::ostringstream remark;

    remark.setf(std::ios::fixed, std::ios::floatfield);
    remark.setf(std::ios::showpoint);

    remark << "REMARK VINA RESULT: " 
           << std::setw(9) << std::setprecision(1) << e
           << "  " << std::setw(9) << std::setprecision(3) << lb
           << "  " << std::setw(9) << std::setprecision(3) << ub
           << '\n';

    return remark.str();
}

output_container Vina::remove_redundant(const output_container& in, fl min_rmsd) {
    output_container tmp;
    VINA_FOR_IN(i, in)
        add_to_output_container(tmp, in[i], min_rmsd, in.size());
    return tmp;
}

void Vina::set_receptor(const std::string& rigid_name) {
    // Read the receptor PDBQT file
    VINA_CHECK(!rigid_name.empty());

    m_receptor = parse_receptor_pdbqt(rigid_name);
    m_receptor_initialized = true;
}

void Vina::set_receptor(const std::string& rigid_name, const std::string& flex_name) {
    // Read the receptor PDBQT file
    VINA_CHECK(!rigid_name.empty());

    if (!flex_name.empty()) {
        m_receptor = parse_receptor_pdbqt(rigid_name, flex_name);
    } else {
        m_receptor = parse_receptor_pdbqt(rigid_name);
    }
    m_receptor_initialized = true;
}

void Vina::set_ligand(const std::string& ligand_name) {
    // Read ligand PDBQT file and add it to the model
    VINA_CHECK(m_receptor_initialized); // m_model
    
    // Replace current model with receptor and reinitialize poses
    m_model = m_receptor;
    output_container m_poses;

    m_model.append(parse_ligand_pdbqt(ligand_name));

    // Store in Vina object
    m_ligand_initialized = true;
}

void Vina::set_ligand(const std::vector<std::string>& ligand_name) {
    // Read ligand PDBQT files and add them to the model
    VINA_CHECK(!ligand_name.empty());
    VINA_CHECK(m_receptor_initialized); // m_model

    // Replace current model with receptor and reinitialize poses
    m_model = m_receptor;
    output_container m_poses;

    VINA_RANGE(i, 0, ligand_name.size())
        m_model.append(parse_ligand_pdbqt(ligand_name[i]));

    // Store in Vina object
    m_ligand_initialized = true;
}

/*
void Vina::set_ligand(OpenBabel::OBMol* mol) {
    // Read ligand PDBQT file and add it to the model
    VINA_CHECK(m_receptor_initialized); // m_model
    
    // Replace current model with receptor and reinitialize poses
    m_model = m_receptor;
    output_container m_poses;

    m_model.append(parse_ligand_pdbqt(mol));

    // Store in Vina object
    m_ligand_initialized = true;
}
*/

/*
void Vina::set_ligand(std::vector<OpenBabel::OBMol*> mol) {
    // Read ligand PDBQT files and add them to the model
    VINA_CHECK(!mol.empty());
    VINA_CHECK(m_receptor_initialized); // m_model

    // Replace current model with receptor and reinitialize poses
    m_model = m_receptor;
    output_container m_poses;

    VINA_RANGE(i, 0, mol.size())
        m_model.append(parse_ligand_pdbqt(mol[i]));

    // Store in Vina object
    m_ligand_initialized = true;
}
*/

void Vina::set_box(double center_x, double center_y, double center_z, int size_x, int size_y, int size_z, double granularity) {
    // Setup the search box
    try {
        if (size_x <= 0 || size_y <= 0 || size_z <= 0) {
            throw "Search space dimensions should be positive";
        }
    } catch (const char* e) {
        m_log << "Exception: " << e << "\n";
        exit (EXIT_FAILURE);
    }

    if (size_x * size_y * size_z > 27e3) {
        m_log << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
    }

    grid_dims gd;
    vec span(size_x, size_y, size_z);
    vec center(center_x, center_y, center_z);

    VINA_FOR_IN(i, gd) {
        gd[i].n = sz(std::ceil(span[i] / granularity));
        fl real_span = granularity * gd[i].n;
        gd[i].begin = center[i] - real_span / 2;
        gd[i].end = gd[i].begin + real_span;
    }

    const vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
    const vec corner2(gd[0].end,   gd[1].end,   gd[2].end);

    // Store in Vina object
    m_gd = gd;
    m_corner1 = corner1;
    m_corner2 = corner2;
    m_box_initialized = true;
}

void Vina::set_vina_weights(const double weight_gauss1, const double weight_gauss2, const double weight_repulsion, 
                            const double weight_hydrophobic, const double weight_hydrogen, const double weight_rot) {
    flv weights;
    weights.push_back(weight_gauss1);
    weights.push_back(weight_gauss2);
    weights.push_back(weight_repulsion);
    weights.push_back(weight_hydrophobic);
    weights.push_back(weight_hydrogen);
    weights.push_back(5 * weight_rot / 0.1 - 1);
    VINA_CHECK(weights.size() == 6);

    // Store in Vina object
    m_weights = weights;
}

void Vina::compute_vina_grid() {
    // Setup the search box
    // Check first that the receptor was added
    // And the box and the ff were defined
    VINA_CHECK(m_receptor_initialized); // m_model
    VINA_CHECK(m_box_initialized); // m_gd

    const fl slope = 1e6; // FIXME: too large? used to be 100
    const szv atom_types_needed = m_model.get_movable_atom_types(m_precalculated_sf.atom_typing_used());

    cache grid("scoring_function_version001", m_gd, slope, atom_type::XS);

    doing(m_verbosity, "Computing grid", m_log);
    grid.populate(m_model, m_precalculated_sf, atom_types_needed);
    done(m_verbosity, m_log);

    // Store in Vina object
    m_grid = grid;
    m_grid_initialized = true;
}

void Vina::write_pose(const std::string& output_name, const std::string& remark) {
    std::ostringstream format_remark;
    format_remark.setf(std::ios::fixed, std::ios::floatfield);
    format_remark.setf(std::ios::showpoint);

    // Add REMARK keyword to be PDB valid
    if(!remark.empty()){
        format_remark << "REMARK " << remark << " \n";
    }

    ofile f(make_path(output_name));
    m_model.write_structure(f, format_remark.str());
}

void Vina::write_results(const std::string& output_name, int how_many, double energy_range) {
    std::string remark;
    model best_model = m_model;
    boost::optional<model> ref;
    int n = 0;
    double best_energy = 0;

    if (how_many < 0) {
        m_log << "WARNING: Number of pose to write set to < 0. Automatically adjusted to default value: 9.\n";
        how_many = 9;
    }

    if (energy_range < 0) {
        m_log << "WARNING: Energy range set to < 0 kcal/mol. Automatically adjusted to default value: 3. kcal/mol.\n";
        energy_range = 3.0;
    }

    if(!m_poses.empty()) {
        m_log.setf(std::ios::fixed, std::ios::floatfield);
        m_log.setf(std::ios::showpoint);
        m_log << '\n';
        m_log << "mode |   affinity | dist from best mode\n";
        m_log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
        m_log << "-----+------------+----------+----------\n";

        // Get the best conf and its energy
        best_model.set(m_poses[0].c);
        best_energy = m_poses[0].e;

        // Open output file
        ofile f(make_path(output_name));

        VINA_FOR_IN(i, m_poses) {
            /* Stop if:
                - We wrote the number of conf asked
                - If there is no conf to write
                - The energy of the current conf is superior than best_energy + energy_range
            */
            if(n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range) break; // check energy_range sanity FIXME

            // Push the current pose to model
            m_model.set(m_poses[i].c);

            // Get RMSD between current pose and best_model
            const model& r = ref ? ref.get() : best_model;
            const fl lb = m_model.rmsd_lower_bound(r);
            const fl ub = m_model.rmsd_upper_bound(r);

            m_log << std::setw(4) << i+1 << "    " << std::setw(9) << std::setprecision(1) << m_poses[i].e; // intermolecular_energies[i];
            m_log << "  " << std::setw(9) << std::setprecision(3) << lb << "  " << std::setw(9) << std::setprecision(3) << ub; // FIXME need user-readable error messages in case of failures
            m_log.endl();

            // Write conf
            remark = vina_remark(m_poses[i].e, lb, ub);
            m_model.write_model(f, n + 1, remark);

            n++;
        }

        // Push back the best conf in model
        m_model.set(m_poses[0].c);

    } else {
        if (m_verbosity > 1) {
            m_log << "WARNING: Could not find any conformations. No conformations were written.";
            m_log.endl();
        }
    }
}

void Vina::randomize(const int max_steps) {
    // Randomize ligand/flex residues conformation
    // Check the box was defined
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_box_initialized); // m_gd, m_corner1, m_corner2

    conf c;
    int seed = generate_seed();
    double penalty = 0;
    double best_clash_penalty = 0;
    std::stringstream sstm;

    rng generator(static_cast<rng::result_type>(seed));

    // It's okay to take the initial conf since we will randomize it
    conf init_conf = m_model.get_initial_conf();
    conf best_conf = init_conf;

    sstm << "Randomize conformation (random seed: " << seed << ")";
    doing(m_verbosity, sstm.str(), m_log);
    VINA_FOR(i, max_steps) {
        c = init_conf;
        c.randomize(m_corner1, m_corner2, generator);
        penalty = m_model.clash_penalty();

        if (i == 0 || penalty < best_clash_penalty) {
            best_conf = c;
            best_clash_penalty = penalty;
        }
    }
    done(m_verbosity, m_log);

    m_model.set(best_conf);

    if (m_verbosity > 1) {
        m_log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
        m_log.endl();
    }
}

void Vina::score_robust() {
    // Check first that the receptor was added
    // And the box and the ff were defined
    VINA_CHECK(m_ligand_initialized); // m_model

    flv term_values;
    double e = 0;
    double e2 = 0;
    double intramolecular_energy = 0;
    const vec authentic_v(1000, 1000, 1000);
    naive_non_cache nnc(&m_precalculated_sf); // for out of grid issues

    intramolecular_energy = m_model.eval_intramolecular(m_precalculated_sf, authentic_v);
    e = m_model.eval_adjusted(m_scoring_function, m_precalculated_sf, nnc, authentic_v, intramolecular_energy);
    m_log << "Intramolecular: " << std::fixed << std::setprecision(5) << intramolecular_energy << " (kcal/mol)";
    m_log.endl();
    m_log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
    m_log.endl();

    term_values = m_t.evale_robust(m_model);
    VINA_CHECK(term_values.size() == 5);
    
    m_log << "Intermolecular contributions to the terms, before weighting:\n";
    m_log << std::setprecision(5);
    m_log << "    gauss 1     : " << term_values[0] << '\n';
    m_log << "    gauss 2     : " << term_values[1] << '\n';
    m_log << "    repulsion   : " << term_values[2] << '\n';
    m_log << "    hydrophobic : " << term_values[3] << '\n';
    m_log << "    Hydrogen    : " << term_values[4] << '\n';

    VINA_CHECK(m_weights.size() == term_values.size() + 1);
    
    VINA_FOR_IN(i, term_values)
        e2 += term_values[i] * m_weights[i];

    e2 = m_scoring_function.conf_independent(m_model, e2);

    m_log << "Affinity: " << std::fixed << std::setprecision(5) << e2 << " (kcal/mol)";
    m_log.endl();
    
    if(e < 100 && std::abs(e2 - e) > 0.05) {
        m_log << "WARNING: the individual terms are inconsisent with the\n";
        m_log << "WARNING: affinity. Consider reporting this as a bug:\n";
        m_log << "WARNING: http://vina.scripps.edu/manual.html#bugs\n";
    }
}

double Vina::score() {
    // Score the current conf in the model
    // Check if ff and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model

    double e = 0;
    double intramolecular_energy = 0;
    double intermolecular_energy = 0;
    const vec authentic_v(1000, 1000, 1000);
    naive_non_cache nnc(&m_precalculated_sf); // for out of grid issues

    intramolecular_energy = m_model.eval_intramolecular(m_precalculated_sf, authentic_v);
    e = m_model.eval_adjusted(m_scoring_function, m_precalculated_sf, nnc, authentic_v, intramolecular_energy);
    
    if (m_verbosity > 1) {
        m_log << "Intramolecular (eval_intramolecular): " << std::fixed << std::setprecision(5) << intramolecular_energy << " (kcal/mol)";
        m_log.endl();
        m_log << "Energy (eval_adjusted): " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
        m_log.endl();
    }

    return e;
}

void Vina::optimize(output_type& out, int max_steps) {
    // Local optimization of the ligand conf
    // Check if ff, box and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_box_initialized); // m_gd

    change g(m_model.get_size());
    quasi_newton quasi_newton_par;
    const fl slope = 1e6;
    const vec authentic_v(1000, 1000, 1000);
    non_cache nc(m_model, m_gd, &m_precalculated_sf, slope);

    // Define the number minimization steps based on the number moving atoms
    if(max_steps == 0) {
        max_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
    }
    quasi_newton_par.max_steps = max_steps;

    if (m_verbosity > 1) {
        m_log << "Before local optimization:";
        m_log.endl();
        score();
    }

    doing(m_verbosity, "Performing local search", m_log);
    VINA_FOR(p, 5) {
        nc.slope = 100 * std::pow(10.0, 2.0 * p);
        quasi_newton_par(m_model, m_precalculated_sf, nc, out, g, authentic_v);
        
        if(nc.within(m_model))
            break;
    }
    done(m_verbosity, m_log);

    out.coords = m_model.get_heavy_atom_movable_coords();

    if(!nc.within(m_model)) {
        m_log << "WARNING: not all movable atoms are within the search space\n";
        out.e = max_fl;
    }

    if (m_verbosity > 1) {
        // Get energy of the new pose
        m_log << "After local optimization:";
        m_log.endl();
        score();
    }
}

void Vina::optimize(int max_steps) {
    // Local optimization of the ligand conf
    // Check if ff, box and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_box_initialized); // m_gd

    // We have to find a way to get rid of this out thing...
    // And also get the current conf and not the initial conf
    double e = 0;
    conf c = m_model.get_initial_conf();
    output_type out(c, e);

    optimize(out, max_steps);
}

void Vina::global_search(const int n_poses, const double min_rmsd) {
    // Vina search (Monte-carlo and local optimization)
    // Check if ff, box and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_grid_initialized); // m_grid

    int seed = generate_seed();
    double e = 0;
    double intramolecular_energy = 0;
    const vec authentic_v(1000, 1000, 1000);
    output_container poses;
    std::stringstream sstm;

    rng generator(static_cast<rng::result_type>(seed));

    // We have to fix scoring function mess by declaring it in the Vina object
    everything t;
    weighted_terms scoring_function(&t, m_weights);

    // Setup Monte-Carlo search
    parallel_mc m_parallelmc;
    sz heuristic = m_model.num_movable_atoms() + 10 * m_model.get_size().num_degrees_of_freedom();
    m_parallelmc.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
    m_parallelmc.mc.ssd_par.evals = unsigned((25 + m_model.num_movable_atoms()) / 3);
    m_parallelmc.mc.min_rmsd = min_rmsd;
    m_parallelmc.mc.num_saved_mins = n_poses;
    m_parallelmc.mc.hunt_cap = vec(10, 10, 10);
    m_parallelmc.num_tasks = m_exhaustiveness;
    m_parallelmc.num_threads = m_cpu;
    m_parallelmc.display_progress = (m_verbosity > 1);

    sstm << "Performing search (random seed: " << seed << ")";
    doing(m_verbosity, sstm.str(), m_log);
    m_parallelmc(m_model, poses, m_precalculated_sf, m_grid, m_precalculated_sf, m_grid, m_corner1, m_corner2, generator);
    done(m_verbosity, m_log);

    doing(m_verbosity, "Refining results", m_log);
    VINA_FOR_IN(i, poses) {
        optimize(poses[i], m_parallelmc.mc.ssd_par.evals);
    }
    done(m_verbosity, m_log);

    doing(m_verbosity, "Removing duplicates", m_log);
    m_poses = remove_redundant(poses, min_rmsd);
    done(m_verbosity, m_log);

    if(!poses.empty()) {        
        VINA_FOR_IN(i, poses) {
            m_model.set(poses[i].c);
            poses[i].e = score();
        }
        
        // the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
        poses.sort();
        // we update model with the best conf
        m_model.set(poses[0].c);
    } else {
        if (m_verbosity > 1) {
            m_log << "WARNING: Could not find any conformations completely within the search space.\n"
              << "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
            m_log.endl();
        }
    }

    // Store results in Vina object
    m_poses = poses;
}


Vina::~Vina() {
    //OpenBabel::OBMol m_mol;
    // model and poses
    model m_receptor;
    model m_model;
    cache m_grid;
    output_container m_poses;
    bool m_receptor_initialized;
    bool m_ligand_initialized;
    // scoring function
    everything m_t;
    flv m_weights;
    non_cache m_nc;
    weighted_terms m_scoring_function;
    precalculate m_precalculated_sf;
    // maps
    grid_dims m_gd;
    vec m_corner1;
    vec m_corner2;
    bool m_no_cache;
    bool m_box_initialized;
    bool m_grid_initialized;
    // global search
    parallel_mc m_parallelmc;
    int m_cpu;
    int m_seed;
    int m_exhaustiveness;
    // others
    int m_verbosity;
    tee m_log;
}