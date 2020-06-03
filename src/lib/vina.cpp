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


void Vina::set_receptor(const std::string& rigid_name) {
    // Read the receptor PDBQT file
    m_m = parse_receptor_pdbqt(rigid_name);
}

void Vina::set_receptor(const std::string& rigid_name, const std::string& flex_name) {
    // Read the receptor PDBQT file
    m_m = parse_receptor_pdbqt(rigid_name, flex_name);
}

void Vina::set_ligand(const std::string& ligand_name) {
    // Read ligand PDBQT file and add it to the model
    m_m.append(parse_ligand_pdbqt(ligand_name));
}

void Vina::set_ligand(const std::vector<std::string>& ligand_name) {
    // Read ligand PDBQT files and add them to the model
    VINA_CHECK(!ligand_name.empty());

    VINA_RANGE(i, 0, ligand_name.size())
        m_m.append(parse_ligand_pdbqt(ligand_name[i]));
}

void Vina::set_box(double center_x, double center_y, double center_z, int size_x, int size_y, int size_z, double granularity = 0.375) {
    // Setup the search box
    vec span(size_x, size_y, size_z);
    vec center(center_x, center_y, center_z);

    VINA_FOR_IN(i, m_gd) {
        m_gd[i].n = sz(std::ceil(span[i] / granularity));
        fl real_span = granularity * m_gd[i].n;
        m_gd[i].begin = center[i] - real_span / 2;
        m_gd[i].end = m_gd[i].begin + real_span;
    }
}

void Vina::set_weights(flv weights) {
    VINA_CHECK(weights.size() == 6);
    m_weights = weights;
}

void Vina::write_all_output(model& m, const output_container& out, sz how_many, const std::string& output_name,
                            const std::vector<std::string>& remarks) {
  if(out.size() < how_many)
    how_many = out.size();
  VINA_CHECK(how_many <= remarks.size());
  ofile f(make_path(output_name));
  VINA_FOR(i, how_many) {
    m.set(out[i].c);
    m.write_model(f, i+1, remarks[i]); // so that model numbers start with 1
  }
}

void Vina::do_randomization(model& m,
            const std::string& out_name,
            const vec& corner1, const vec& corner2, int seed, int verbosity, tee& log) {
  conf init_conf = m.get_initial_conf();
  rng generator(static_cast<rng::result_type>(seed));
  if(verbosity > 1) {
    log << "Using random seed: " << seed;
    log.endl();
  }
  const sz attempts = 10000;
  conf best_conf = init_conf;
  fl best_clash_penalty = 0;
  VINA_FOR(i, attempts) {
    conf c = init_conf;
    c.randomize(corner1, corner2, generator);
    m.set(c);
    fl penalty = m.clash_penalty();
    if(i == 0 || penalty < best_clash_penalty) {
      best_conf = c;
      best_clash_penalty = penalty;
    }
  }
  m.set(best_conf);
  if(verbosity > 1) {
    log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
    log.endl();
  }
  m.write_structure(make_path(out_name));
}

void Vina::refine_structure(model& m, const precalculate& prec, non_cache& nc, output_type& out, const vec& cap, sz max_steps = 1000) {
  change g(m.get_size());
  quasi_newton quasi_newton_par;
  quasi_newton_par.max_steps = max_steps;
  const fl slope_orig = nc.slope;
  VINA_FOR(p, 5) {
    nc.slope = 100 * std::pow(10.0, 2.0*p);
    quasi_newton_par(m, prec, nc, out, g, cap);
    m.set(out.c); // just to be sure
    if(nc.within(m))
      break;
  }
  out.coords = m.get_heavy_atom_movable_coords();
  if(!nc.within(m))
    out.e = max_fl;
  nc.slope = slope_orig;
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

void Vina::do_search(model& m, const scoring_function& sf, const precalculate& prec, 
                     const igrid& ig, const precalculate& prec_widened, const igrid& ig_widened, non_cache& nc, // nc.slope is changed
                     const std::string& out_name, const vec& corner1, const vec& corner2, const parallel_mc& par, 
                     fl energy_range, sz num_modes, int seed, int verbosity, bool score_only, bool local_only, 
                     tee& log, const terms& t, const flv& weights) {
  conf_size s = m.get_size();
  conf c = m.get_initial_conf();
  fl e = max_fl;
  boost::optional<model> ref;
  const vec authentic_v(1000, 1000, 1000);

  if(score_only) {
    fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, c);
    naive_non_cache nnc(&prec); // for out of grid issues
    e = m.eval_adjusted(sf, prec, nnc, authentic_v, c, intramolecular_energy);
    log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
    log.endl();
    flv term_values = t.evale_robust(m);
    VINA_CHECK(term_values.size() == 5);
    log << "Intermolecular contributions to the terms, before weighting:\n";
    log << std::setprecision(5);
    log << "    gauss 1     : " << term_values[0] << '\n';
    log << "    gauss 2     : " << term_values[1] << '\n';
    log << "    repulsion   : " << term_values[2] << '\n';
    log << "    hydrophobic : " << term_values[3] << '\n';
    log << "    Hydrogen    : " << term_values[4] << '\n';
    VINA_CHECK(weights.size() == term_values.size() + 1);
    fl e2 = 0;
    VINA_FOR_IN(i, term_values)
      e2 += term_values[i] * weights[i];
    e2 = sf.conf_independent(m, e2);
    if(e < 100 && std::abs(e2 - e) > 0.05) {
      log << "WARNING: the individual terms are inconsisent with the\n";
      log << "WARNING: affinity. Consider reporting this as a bug:\n";
      log << "WARNING: http://vina.scripps.edu/manual.html#bugs\n";
    }
  }
  else if(local_only) {
    output_type out(c, e);
    doing(verbosity, "Performing local search", log);
    refine_structure(m, prec, nc, out, authentic_v, par.mc.ssd_par.evals);
    done(verbosity, log);
    fl intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out.c);
    e = m.eval_adjusted(sf, prec, nc, authentic_v, out.c, intramolecular_energy);

    log << "Affinity: " << std::fixed << std::setprecision(5) << e << " (kcal/mol)";
    log.endl();
    if(!nc.within(m))
      log << "WARNING: not all movable atoms are within the search space\n";

    doing(verbosity, "Writing output", log);
    output_container out_cont;
    out_cont.push_back(new output_type(out));
    std::vector<std::string> remarks(1, vina_remark(e, 0, 0));
    write_all_output(m, out_cont, 1, out_name, remarks); // how_many == 1
    done(verbosity, log);
  }
  else {
    rng generator(static_cast<rng::result_type>(seed));
    log << "Using random seed: " << seed;
    log.endl();
    output_container out_cont;
    doing(verbosity, "Performing search", log);
    par(m, out_cont, prec, ig, prec_widened, ig_widened, corner1, corner2, generator);
    done(verbosity, log);

    doing(verbosity, "Refining results", log);
    VINA_FOR_IN(i, out_cont)
      refine_structure(m, prec, nc, out_cont[i], authentic_v, par.mc.ssd_par.evals);

    if(!out_cont.empty()) {
      out_cont.sort();
      const fl best_mode_intramolecular_energy = m.eval_intramolecular(prec, authentic_v, out_cont[0].c);
      VINA_FOR_IN(i, out_cont)
        if(not_max(out_cont[i].e))
          out_cont[i].e = m.eval_adjusted(sf, prec, nc, authentic_v, out_cont[i].c, best_mode_intramolecular_energy); 
      // the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
      out_cont.sort();
    }

    const fl out_min_rmsd = 1;
    out_cont = remove_redundant(out_cont, out_min_rmsd);

    done(verbosity, log);

    log.setf(std::ios::fixed, std::ios::floatfield);
    log.setf(std::ios::showpoint);
    log << '\n';
    log << "mode |   affinity | dist from best mode\n";
    log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
    log << "-----+------------+----------+----------\n";

    model best_mode_model = m;
    if(!out_cont.empty())
      best_mode_model.set(out_cont.front().c);

    sz how_many = 0;
    std::vector<std::string> remarks;
    VINA_FOR_IN(i, out_cont) {
      if(how_many >= num_modes || !not_max(out_cont[i].e) || out_cont[i].e > out_cont[0].e + energy_range) break; // check energy_range sanity FIXME
      ++how_many;
      log << std::setw(4) << i+1
        << "    " << std::setw(9) << std::setprecision(1) << out_cont[i].e; // intermolecular_energies[i];
      m.set(out_cont[i].c);
      const model& r = ref ? ref.get() : best_mode_model;
      const fl lb = m.rmsd_lower_bound(r);
      const fl ub = m.rmsd_upper_bound(r);
      log << "  " << std::setw(9) << std::setprecision(3) << lb
          << "  " << std::setw(9) << std::setprecision(3) << ub; // FIXME need user-readable error messages in case of failures

      remarks.push_back(vina_remark(out_cont[i].e, lb, ub));
      log.endl();
    }
    doing(verbosity, "Writing output", log);
    write_all_output(m, out_cont, how_many, out_name, remarks);
    done(verbosity, log);

    if(how_many < 1) {
      log << "WARNING: Could not find any conformations completely within the search space.\n"
        << "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
      log.endl();
    }
  }
}

void Vina::run(const std::string& out_name) {
  doing(m_verbosity, "Setting up the scoring function", m_log);

  everything t;

  weighted_terms wt(&t, m_weights);
  precalculate prec(wt);
  const fl left  = 0.25; 
  const fl right = 0.25;
  const fl slope = 1e6; // FIXME: too large? used to be 100
  precalculate prec_widened(prec); prec_widened.widen(left, right);

  done(m_verbosity, m_log);

  vec corner1(m_gd[0].begin, m_gd[1].begin, m_gd[2].begin);
  vec corner2(m_gd[0].end,   m_gd[1].end,   m_gd[2].end);

  parallel_mc par;
  sz heuristic = m_m.num_movable_atoms() + 10 * m_m.get_size().num_degrees_of_freedom();
  par.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
  par.mc.ssd_par.evals = unsigned((25 + m_m.num_movable_atoms()) / 3);
  par.mc.min_rmsd = 1.0;
  par.mc.num_saved_mins = 20;
  par.mc.hunt_cap = vec(10, 10, 10);
  par.num_tasks = m_exhaustiveness;
  par.num_threads = m_cpu;
  par.display_progress = (m_verbosity > 1);

  if(m_randomize_only) {
    do_randomization(m_m, out_name, corner1, corner2, m_seed, m_verbosity, m_log);
  } else {
    non_cache nc        (m_m, m_gd, &prec,         slope); // if gd has 0 n's, this will not constrain anything
    non_cache nc_widened(m_m, m_gd, &prec_widened, slope); // if gd has 0 n's, this will not constrain anything
    
    if(m_no_cache) {
      do_search(m_m, wt, prec, nc, prec_widened, nc_widened, nc,
                out_name,
                corner1, corner2,
                par, m_energy_range, m_num_modes,
                m_seed, m_verbosity, m_score_only, m_local_only, m_log, t, m_weights);
    }
    else {
      bool cache_needed = !(m_score_only || m_randomize_only || m_local_only);
      if(cache_needed) doing(m_verbosity, "Analyzing the binding site", m_log);
      cache c("scoring_function_version001", m_gd, slope, atom_type::XS);
      if(cache_needed) c.populate(m_m, prec, m_m.get_movable_atom_types(prec.atom_typing_used()));
      if(cache_needed) done(m_verbosity, m_log);
      do_search(m_m, wt, prec, c, prec, c, nc,
                out_name,
                corner1, corner2,
                par, m_energy_range, m_num_modes,
                m_seed, m_verbosity, m_score_only, m_local_only, m_log, t, m_weights);
    }
  }
}
