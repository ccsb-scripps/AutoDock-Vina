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

#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "ssd.h"
#include "incrementable.h"

struct search_stats {
    unsigned done_evals;
    unsigned done_steps;
    unsigned count_accepted;
    unsigned count_rejected;
    unsigned bfgs_reject;
    unsigned bfgs_accept;
    unsigned best_mc_step; // MC step at which best was found (1-indexed)
    unsigned nr_promising; // number of re-minimized best solutions
    std::vector<unsigned> hist_linesearch; // line search trials 
    std::vector<unsigned> hist_bfgs; // quasi-newton steps
};


struct monte_carlo {
	unsigned evalcap;
	unsigned num_steps;
	fl temperature;
	vec hunt_cap;
	fl min_rmsd;
	sz num_saved_mins;
	fl mutation_amplitude;
	ssd ssd_par;
	monte_carlo() : evalcap(0), num_steps(2500), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2) {} // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  num_steps = 50*lig_atoms = 2500

	output_type operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	output_type many_runs(model& m, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const;

	void single_run(model& m, output_type& out, search_stats& stats, const precalculate& p, const igrid& ig, rng& generator) const;
	// out is sorted
	void operator()(model& m, output_container& out, search_stats& stats, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	void many_runs(model& m, output_container& out, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const;

};

#endif

