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

#include "monte_carlo.h"
#include "coords.h"
#include "mutate.h"
#include "quasi_newton.h"

output_type monte_carlo::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	output_container tmp;
	search_stats stats;
	this->operator()(m, tmp, stats, p, ig, p_widened, ig_widened, corner1, corner2, increment_me, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

void monte_carlo::single_run(model& m, output_type& out, search_stats& stats, const precalculate& p, const igrid& ig, rng& generator) const {
    int evalcount = 0;
    printf("monte_carlo::single_run\n");
	conf_size s = m.get_size();
	change g(s);
	vec authentic_v(1000, 1000, 1000);
	out.e = max_fl;
	output_type current(out);
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	VINA_U_FOR(step, num_steps) {
		output_type candidate(current.c, max_fl);
		mutate_conf(candidate.c, m, mutation_amplitude, generator);
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap, evalcount,
                stats.bfgs_reject, stats.bfgs_accept,
                stats.hist_bfgs, stats.hist_linesearch);
		if(step == 0 || metropolis_accept(current.e, candidate.e, temperature, generator)) {
			quasi_newton_par(m, p, ig, candidate, g, authentic_v, evalcount,
                stats.bfgs_reject, stats.bfgs_accept,
                stats.hist_bfgs, stats.hist_linesearch);
			current = candidate;
			if(current.e < out.e)
				out = current;
		}
	}
	quasi_newton_par(m, p, ig, out, g, authentic_v, evalcount,
                stats.bfgs_reject, stats.bfgs_accept,
                stats.hist_bfgs, stats.hist_linesearch);
}

void monte_carlo::many_runs(model& m, output_container& out, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const { // This is never used / untested
	conf_size s = m.get_size();
	VINA_FOR(run, num_runs) {
		output_type tmp(s, 0);
        search_stats stats;
		tmp.c.randomize(corner1, corner2, generator);
		single_run(m, tmp, stats, p, ig, generator);
		out.push_back(new output_type(tmp));
	}
	out.sort();

    printf("monte_carlo::many_runs\n");
}

output_type monte_carlo::many_runs(model& m, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const {
	output_container tmp;
	many_runs(m, tmp, p, ig, corner1, corner2, num_runs, generator);
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}


// out is sorted
void monte_carlo::operator()(model& m, output_container& out, search_stats& stats, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
    int evalcount = 0;
    int verypromising = 0;
    stats.done_steps = 0;
    stats.done_evals = 0;
    stats.count_rejected = 0;
    stats.count_accepted = 0;
    stats.bfgs_reject = 0;
    stats.bfgs_accept = 0;
    stats.nr_promising = 0;
    VINA_FOR(i, 10) stats.hist_linesearch.push_back(0);
    VINA_FOR(i, ssd_par.evals) stats.hist_bfgs.push_back(0);
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	conf_size s = m.get_size();
	change g(s);
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);
	fl best_e = max_fl;
	quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
	VINA_U_FOR(step, num_steps) {

		if(increment_me)
			++(*increment_me);

        if ((evalcap > 0) & (evalcount > evalcap)) {
            //printf("BREAKING evalcap!\n");
            break;
        } 

        stats.done_steps += 1;

		output_type candidate = tmp;
		mutate_conf(candidate.c, m, mutation_amplitude, generator);
		quasi_newton_par(m, p, ig, candidate, g, hunt_cap, evalcount,
                stats.bfgs_reject, stats.bfgs_accept,
                stats.hist_bfgs, stats.hist_linesearch);
		if(step == 0 || metropolis_accept(tmp.e, candidate.e, temperature, generator)) {
            stats.count_accepted += 1;
			tmp = candidate;

			m.set(tmp.c); // FIXME? useless?

			// FIXME only for very promising ones
			if(tmp.e < best_e || out.size() < num_saved_mins) {
                verypromising++;
                stats.nr_promising++;
				quasi_newton_par(m, p, ig, tmp, g, authentic_v, evalcount,
                                 stats.bfgs_reject, stats.bfgs_accept,
                                 stats.hist_bfgs, stats.hist_linesearch);
				m.set(tmp.c); // FIXME? useless?
				tmp.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
				if(tmp.e < best_e)
					best_e = tmp.e;
                    stats.best_mc_step = step + 1; // 1-indexed
			}
        }
		else {
            stats.count_rejected += 1;
            //printf("Metropolis REJECT, step: %6d\n", step);
        }
	}
    stats.done_evals = evalcount;
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
