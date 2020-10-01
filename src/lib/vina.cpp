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
#include "scoring_function.h"
#include "precalculate.h"


void Vina::cite() {
    const std::string cite_message = "\
#################################################################\n\
# If you used AutoDock Vina in your work, please cite:          #\n\
#                                                               #\n\
# O. Trott, A. J. Olson,                                        #\n\
# AutoDock Vina: improving the speed and accuracy of docking    #\n\
# with a new scoring function, efficient optimization and       #\n\
# multithreading, Journal of Computational Chemistry 31 (2010)  #\n\
# 455-461                                                       #\n\
#                                                               #\n\
# DOI 10.1002/jcc.21334                                         #\n\
#                                                               #\n\
# Please see http://vina.scripps.edu for more information.      #\n\
#################################################################\n";

    std::cout << cite_message << '\n';
}

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

std::string Vina::vina_remarks(output_type& pose, fl lb, fl ub) {
    std::ostringstream remark;

    remark.setf(std::ios::fixed, std::ios::floatfield);
    remark.setf(std::ios::showpoint);

    remark << "REMARK VINA RESULT: " 
           << std::setw(9) << std::setprecision(1) << pose.e
           << "  " << std::setw(9) << std::setprecision(3) << lb
           << "  " << std::setw(9) << std::setprecision(3) << ub
           << '\n';
    
    remark << "REMARK ITER + INTRA:     " << std::setw(12) << std::setprecision(4) << pose.total << "\n";
    remark << "REMARK INTER:            " << std::setw(12) << std::setprecision(4) << pose.inter << "\n";
    remark << "REMARK INTRA:            " << std::setw(12) << std::setprecision(4) << pose.intra << "\n";
    if(m_sf_choice==SF_AD42) remark << "REMARK CONF_INDEPENDENT: " << std::setw(12) << std::setprecision(4) << pose.conf_independent << "\n";
    remark << "REMARK UNBOUND:          " << std::setw(12) << std::setprecision(4) << pose.unbound << "\n";

    return remark.str();
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

    m_receptor = parse_receptor_pdbqt(rigid_name, m_scoring_function.get_atom_typing());
    m_model = m_receptor;
    m_receptor_initialized = true;
}

void Vina::set_receptor(const std::string& rigid_name, const std::string& flex_name) {
    // Read the receptor PDBQT file
    VINA_CHECK(!rigid_name.empty());

    if (!flex_name.empty()) {
        m_receptor = parse_receptor_pdbqt(rigid_name, flex_name, m_scoring_function.get_atom_typing());
    } else {
        m_receptor = parse_receptor_pdbqt(rigid_name, m_scoring_function.get_atom_typing());
    }

    m_model = m_receptor;
    m_receptor_initialized = true;
}

void Vina::set_ligand(const std::string& ligand_name) {
    // Read ligand PDBQT file and add it to the model
    VINA_CHECK(m_receptor_initialized); // m_model
    
    // Replace current model with receptor and reinitialize poses
    m_model = m_receptor;
    output_container poses;
    // ... and add ligand to the model
    m_model.append(parse_ligand_pdbqt(ligand_name, m_scoring_function.get_atom_typing()));
    m_model.about();

    // Because we precalculate ligand atoms interactions
    precalculate_byatom precalculated_byatom(m_scoring_function, m_model);

    // Store in Vina object
    m_poses = poses;
    m_precalculated_byatom = precalculated_byatom;
    m_ligand_initialized = true;
}

void Vina::set_ligand(const std::vector<std::string>& ligand_name) {
    // Read ligand PDBQT files and add them to the model
    VINA_CHECK(!ligand_name.empty());
    VINA_CHECK(m_receptor_initialized); // m_model

    // Replace current model with receptor and reinitialize poses
    m_model = m_receptor;
    output_container poses;

    VINA_RANGE(i, 0, ligand_name.size())
        m_model.append(parse_ligand_pdbqt(ligand_name[i], m_scoring_function.get_atom_typing()));

    // Because we precalculate ligand atoms interactions
    precalculate_byatom precalculated_byatom(m_scoring_function, m_model);

    // Store in Vina object
    m_poses = poses;
    m_precalculated_byatom = precalculated_byatom;
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

void Vina::set_vina_weights(double weight_gauss1, double weight_gauss2, double weight_repulsion, 
                            double weight_hydrophobic, double weight_hydrogen, double weight_glue,
                            double weight_rot) {
    flv weights;

    if (m_sf_choice == SF_VINA) {
        weights.push_back(weight_gauss1);
        weights.push_back(weight_gauss2);
        weights.push_back(weight_repulsion);
        weights.push_back(weight_hydrophobic);
        weights.push_back(weight_hydrogen);
        weights.push_back(weight_glue);
        weights.push_back(5 * weight_rot / 0.1 - 1);

        // Store in Vina object
        m_weights = weights;

        // Since we set (different) weights, we automatically initialize the forcefield
        set_forcefield();
    }
}

void Vina::set_ad4_weights(double weight_ad4_vdw , double weight_ad4_hb, 
                           double weight_ad4_elec, double weight_ad4_dsolv, 
                           double weight_glue, double weight_ad4_rot) {
    flv weights;

    if (m_sf_choice == SF_AD42) {
        weights.push_back(weight_ad4_vdw);
        weights.push_back(weight_ad4_hb);
        weights.push_back(weight_ad4_elec);
        weights.push_back(weight_ad4_dsolv);
        weights.push_back(weight_glue);
        weights.push_back(weight_ad4_rot);

        // Store in Vina object
        m_weights = weights;

        // Since we set (different) weights, we automatically initialize the forcefield
        set_forcefield();
    }
}

void Vina::set_forcefield() 
{
    ScoringFunction scoring_function(m_sf_choice, m_weights);
    // Store in Vina object
    m_scoring_function = scoring_function;
}

void Vina::set_vina_box(double center_x, double center_y, double center_z, double size_x, double size_y, double size_z, double granularity) {
    // Setup the search box
    if (size_x <= 0 || size_y <= 0 || size_z <= 0) {
        std::cerr << "ERROR: Search space dimensions should be positive.\n";
        exit(EXIT_FAILURE);
    }

    if (size_x * size_y * size_z > 27e3) {
        std::cerr << "WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)\n";
    }

    grid_dims gd;
    vec span(size_x, size_y, size_z);
    vec center(center_x, center_y, center_z);

    VINA_FOR_IN(i, gd) {
        // Make sure that the number of grid points is an odd number
        sz n = std::ceil(span[i] / granularity);
        gd[i].n = std::ceil(n / 2) * 2 + 1;
        fl real_span = granularity * (gd[i].n - 1) / 2.;
        gd[i].begin = center[i] - real_span;
        gd[i].end = center[i] + real_span;
    }

    const vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
    const vec corner2(gd[0].end,   gd[1].end,   gd[2].end);

    // Store in Vina object
    m_gd = gd;
    m_corner1 = corner1;
    m_corner2 = corner2;
}

void Vina::compute_vina_maps(double center_x, double center_y, double center_z, double size_x, double size_y, double size_z, double granularity) {
    // Setup the search box
    // Check first that the receptor was added
    // And the box and the ff were defined
    VINA_CHECK(m_receptor_initialized); // m_model

    const fl slope = 1e6; // FIXME: too large? used to be 100
    szv atom_types;
    atom_type::t atom_typing = m_scoring_function.get_atom_typing();

    if (m_ligand_initialized)
        atom_types = m_model.get_movable_atom_types(atom_typing);
    else
        atom_types = m_scoring_function.get_atom_types();

    doing("Computing Vina grid", m_verbosity, 0);
    // Initialize the vina box
    set_vina_box(center_x, center_y, center_z, size_x, size_y, size_z, granularity);
    precalculate precalculated_sf(m_scoring_function);
    cache grid(slope);
    grid.populate(m_model, precalculated_sf, m_gd, atom_types);
    done(m_verbosity, 0);

    // Store in Vina object
    m_grid = grid;
    m_map_initialized = true;
}

void Vina::load_maps(std::string maps) {
    const fl slope = 1e6; // FIXME: too large? used to be 100
    grid_dims gd;

    if (m_sf_choice == SF_VINA) {
        doing("Reading Vina maps", m_verbosity, 0);
        cache grid(slope);
        gd = grid.read(maps);
        done(m_verbosity, 0);
        m_grid = grid;
    } else {
        doing("Reading AD4.2 maps", m_verbosity, 0);
        ad4cache grid(slope);
        gd = grid.read(maps);
        done(m_verbosity, 0);
        m_ad4grid = grid;
    }

    // Store in Vina object
    const vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
    const vec corner2(gd[0].end, gd[1].end, gd[2].end);
    m_gd = gd;
    m_corner1 = corner1;
    m_corner2 = corner2;
    m_map_initialized = true;
}

void Vina::write_maps(const std::string& map_prefix, const std::string& gpf_filename,
                      const std::string& fld_filename, const std::string& receptor_filename) {
    VINA_CHECK(m_map_initialized); // // m_gd, m_corner1, m_corner2, m_grid/ad4grid

    szv atom_types;
    atom_type::t atom_typing = m_scoring_function.get_atom_typing();

    if (m_ligand_initialized)
        atom_types = m_model.get_movable_atom_types(atom_typing);
    else
        atom_types = m_scoring_function.get_atom_types();

    if (m_sf_choice == SF_VINA) {
        doing("Writing Vina maps", m_verbosity, 0);
        m_grid.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename);
        done(m_verbosity, 0);
    } else {
        // Add electrostatics and desolvation maps
        atom_types.push_back(AD_TYPE_SIZE);
        atom_types.push_back(AD_TYPE_SIZE + 1);
        doing("Writing AD4.2 maps", m_verbosity, 0);
        m_ad4grid.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename);
        done(m_verbosity, 0);
    }
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
    std::string remarks;
    model best_model = m_model;
    boost::optional<model> ref;
    int n = 0;
    double best_energy = 0;

    if (how_many < 0) {
        std::cerr << "WARNING: Number of pose to write set to < 0. Automatically adjusted to default value: 9.\n";
        how_many = 9;
    }

    if (energy_range < 0) {
        std::cerr << "WARNING: Energy range set to < 0 kcal/mol. Automatically adjusted to default value: 3. kcal/mol.\n";
        energy_range = 3.0;
    }

    if (!m_poses.empty()) {
        if (m_verbosity > 0) {
            std::cout << '\n';
            std::cout << "mode |   affinity | dist from best mode\n";
            std::cout << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
            std::cout << "-----+------------+----------+----------\n";
        }

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
            if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range) break; // check energy_range sanity FIXME

            // Push the current pose to model
            m_model.set(m_poses[i].c);

            // Get RMSD between current pose and best_model
            const model& r = ref ? ref.get() : best_model;
            const fl lb = m_model.rmsd_lower_bound(r);
            const fl ub = m_model.rmsd_upper_bound(r);

            if (m_verbosity > 0) {
                std::cout << std::setw(4) << i+1 << "    " << std::setw(9) << std::setprecision(1) << m_poses[i].e; // intermolecular_energies[i];
                std::cout << "  " << std::setw(9) << std::setprecision(3) << lb << "  " << std::setw(9) << std::setprecision(3) << ub << "\n"; // FIXME need user-readable error messages in case of failures
            }

            // Write conf
            remarks = vina_remarks(m_poses[i], lb, ub);
            m_model.write_model(f, n + 1, remarks);

            n++;
        }

        // Push back the best conf in model
        m_model.set(m_poses[0].c);
    } else {
        std::cerr << "WARNING: Could not find any conformations. No conformations were written.\n";
    }
}

void Vina::randomize(const int max_steps) {
    // Randomize ligand/flex residues conformation
    // Check the box was defined
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_map_initialized); // // m_gd, m_corner1, m_corner2, m_grid/ad4grid

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
    doing(sstm.str(), m_verbosity, 0);
    VINA_FOR(i, max_steps) {
        c = init_conf;
        c.randomize(m_corner1, m_corner2, generator);
        penalty = m_model.clash_penalty();

        if (i == 0 || penalty < best_clash_penalty) {
            best_conf = c;
            best_clash_penalty = penalty;
        }
    }
    done(m_verbosity, 0);

    m_model.set(best_conf);

    if (m_verbosity > 1) {
        std::cout << "Clash penalty: " << best_clash_penalty << "\n";
    }
}

void Vina::show_score(const std::vector<double> energies) {
    std::cout << "Estimated Free Energy of Binding   : " << std::fixed << std::setprecision(5) << energies[0] << " (kcal/mol) [=(1)+(2)+(3)+(4)]\n";
    std::cout << "(1) Final Intermolecular Energy    : " << std::fixed << std::setprecision(5) << energies[1] + energies[2] + energies[3] << " (kcal/mol)\n";
    std::cout << "    Ligand                         : " << std::fixed << std::setprecision(5) << energies[1] << " (kcal/mol)\n";
    std::cout << "    Flexible side chains           : " << std::fixed << std::setprecision(5) << energies[2] << " (kcal/mol)\n";
    std::cout << "    Other                          : " << std::fixed << std::setprecision(5) << energies[3] << " (kcal/mol)\n";
    std::cout << "(2) Final Total Internal Energy    : " << std::fixed << std::setprecision(5) << energies[4] << " (kcal/mol)\n";
    std::cout << "(3) Torsional Free Energy          : " << std::fixed << std::setprecision(5) << energies[5] << " (kcal/mol)\n";
    if (m_sf_choice == SF_VINA) {
        std::cout << "(4) Unbound System's Energy        : " << std::fixed << std::setprecision(5) << energies[6] << " (kcal/mol)\n";
    } else {
        std::cout << "(4) Unbound System's Energy [=(2)] : " << std::fixed << std::setprecision(5) << energies[6] << " (kcal/mol)\n";
    }
}

std::vector<double> Vina::score(double intramolecular_energy) {
    // Score the current conf in the model
    // Check if ff and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_map_initialized); // // m_gd, m_corner1, m_corner2, m_grid/ad4grid

    double total = 0;
    double lig_grids = 0;
    double flex_grids = 0;
    double lig_intra = 0;
    double conf_independent = 0;
    double other_pairs = 0; // TODO
    const vec authentic_v(1000, 1000, 1000);
    std::vector<double> energies;

    if (m_sf_choice == SF_VINA) {
        lig_grids = m_grid.eval(m_model, authentic_v[1]); // [1]
        flex_grids = m_grid.eval_intra(m_model, authentic_v[1]); // [1]
        other_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1] TODO missing other pairs
        lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2]
        total = m_scoring_function.conf_independent(m_model, lig_grids + flex_grids + other_pairs + lig_intra - intramolecular_energy); // we pass intermolecular energy, got to remove flex-flex
        // We want to know how much penalty conf_independent function added to the energy
        conf_independent = total - (lig_grids + flex_grids + other_pairs + lig_intra - intramolecular_energy);
    } else {
        lig_grids = m_ad4grid.eval(m_model, authentic_v[1]); // [1]
        flex_grids = m_ad4grid.eval_intra(m_model, authentic_v[1]); // [1]
        other_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1]
        lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2]
        conf_independent = m_scoring_function.conf_independent(m_model, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
        total = lig_grids + flex_grids + other_pairs + conf_independent; // (+ lig_intra - lig_intra) TODO verify that all contributions from other_pairs are valid, got to remove flex-flex
    }

    energies.push_back(total);
    energies.push_back(lig_grids);
    energies.push_back(flex_grids);
    energies.push_back(other_pairs);
    energies.push_back(lig_intra);
    energies.push_back(conf_independent);

    if (m_sf_choice == SF_VINA) {
        energies.push_back(intramolecular_energy);
    } else {
        energies.push_back(-lig_intra);
    }

    return energies;
}

std::vector<double> Vina::score() {
    // Score the current conf in the model
    // Check if ff and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_map_initialized); // // m_gd, m_corner1, m_corner2, m_grid/ad4grid

    double intramolecular_energy = 0;
    const vec authentic_v(1000, 1000, 1000);

    if(m_sf_choice == SF_VINA) {
        intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_grid, authentic_v);
    }

    std::vector<double> energies = score(intramolecular_energy);
    return energies;
}

std::vector<double> Vina::optimize(output_type& out, int max_steps) {
    // Local optimization of the ligand conf
    // Check if ff, box and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_map_initialized); // m_gd, m_corner1, m_corner2, m_grid/ad4grid

    change g(m_model.get_size());
    quasi_newton quasi_newton_par;
    const fl slope = 1e6;
    const vec authentic_v(1000, 1000, 1000);
    std::vector<double> energies_before_opt;
    std::vector<double> energies_after_opt;
    int evalcount = 0;

    // Define the number minimization steps based on the number moving atoms
    if (max_steps == 0) {
        max_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
        if (m_verbosity > 1)
            std::cout << "Number of local optimization steps: " << max_steps << "\n";
    }
    quasi_newton_par.max_steps = max_steps;

    if (m_verbosity > 0) {
        energies_before_opt = score();
    }

    doing("Performing local search", m_verbosity, 0);
    VINA_FOR(p, 5) {
        if (m_sf_choice == SF_VINA) quasi_newton_par(m_model, m_precalculated_byatom, m_grid,    out, g, authentic_v, evalcount);
        if (m_sf_choice == SF_AD42) quasi_newton_par(m_model, m_precalculated_byatom, m_ad4grid, out, g, authentic_v, evalcount);
    }
    done(m_verbosity, 0);

    energies_after_opt = score();
    out.coords = m_model.get_heavy_atom_movable_coords();

    if (m_verbosity > 1) {
        std::cout << "Before local optimization:";
        show_score(energies_before_opt);
        std::cout << "After local optimization:";
        show_score(energies_after_opt);
    }

    return energies_after_opt;
}

std::vector<double> Vina::optimize(int max_steps) {
    // Local optimization of the ligand conf
    // Check if ff, box and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_map_initialized); // // m_gd, m_corner1, m_corner2, m_grid/ad4grid

    // We have to find a way to get rid of this out thing...
    // And also get the current conf and not the initial conf
    double e = 0;
    conf c = m_model.get_initial_conf();
    output_type out(c, e);

    std::vector<double> energies = optimize(out, max_steps);

    return energies;
}

void Vina::global_search(const int exhaustiveness, const int n_poses, const double min_rmsd, const int max_evals)
{
    // Vina search (Monte-carlo and local optimization)
    // Check if ff, box and ligand were initialized
    VINA_CHECK(m_ligand_initialized); // m_model
    VINA_CHECK(m_map_initialized); // // m_gd, m_corner1, m_corner2, m_grid/ad4grid

    if (exhaustiveness < 1) {
        std::cerr << "ERROR: Exhaustiveness must be 1 or greater";
        exit(EXIT_FAILURE);
    }

    if (exhaustiveness < m_cpu) {
        std::cerr << "WARNING: At low exhaustiveness, it may be impossible to utilize all CPUs.\n";
    }

    int seed = generate_seed();
    double e = 0;
    double intramolecular_energy = 0;
    const vec authentic_v(1000, 1000, 1000);
    output_container poses;
    std::stringstream sstm;
    rng generator(static_cast<rng::result_type>(seed));

    // Setup Monte-Carlo search
    parallel_mc parallelmc;
    sz heuristic = m_model.num_movable_atoms() + 10 * m_model.get_size().num_degrees_of_freedom();
    parallelmc.mc.global_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
    parallelmc.mc.local_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
    parallelmc.mc.max_evals = max_evals;
    parallelmc.mc.min_rmsd = min_rmsd;
    parallelmc.mc.num_saved_mins = n_poses;
    parallelmc.mc.hunt_cap = vec(10, 10, 10);
    parallelmc.num_tasks = exhaustiveness;
    parallelmc.num_threads = m_cpu;
    parallelmc.display_progress = (m_verbosity > 0);

    sstm << "Performing search (random seed: " << seed << ")";
    doing(sstm.str(), m_verbosity, 1);
    if (m_sf_choice == SF_VINA) {
        parallelmc(m_model, poses, m_precalculated_byatom,    m_grid, m_corner1, m_corner2, generator);
    } else {
        parallelmc(m_model, poses, m_precalculated_byatom, m_ad4grid, m_corner1, m_corner2, generator);
    }
    done(m_verbosity, 1);

    m_poses = remove_redundant(poses, min_rmsd);

    if(!poses.empty()) {        
        // the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
        poses.sort();
        // we update model with the best conf
        m_model.set(poses[0].c);

        // For the Vina scoring function, we take the intramolecular energy from the best pose
        if (m_sf_choice == SF_VINA) 
            intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_grid, authentic_v);

        VINA_FOR_IN(i, poses) {
            if (m_verbosity > 1)
                std::cout << "ENERGY FROM SEARCH: " << poses[i].e << "\n";

            m_model.set(poses[i].c);

            // For AD42 intramolecular_energy is equal to 0
            std::vector<double> energies = score(intramolecular_energy);
            // Store energy components in current pose
            poses[i].e = energies[0];
            poses[i].inter = energies[1];
            poses[i].intra = energies[4];
            poses[i].total = energies[1] + energies[2] + energies[3] + energies[4]; // cost function for optimization
            poses[i].conf_independent = energies[5];
            poses[i].unbound = energies[6];

            if (m_verbosity > 1)
                show_score(energies);
        }
    } else {
        std::cerr << "WARNING: Could not find any conformations completely within the search space.\n";
        std::cerr << "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.\n";
    }

    // Store results in Vina object
    m_poses = poses;
}

Vina::~Vina() {
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
    // others
    int m_verbosity;
}
