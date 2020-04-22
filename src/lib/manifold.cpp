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

#include "manifold.h"
#include "recent_history.h"
#include "coords.h"
#include "quasi_newton.h"

const vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl

output_type manifold::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

const conf* select_rs(const output_container& mf, sz rstart, rng& generator) { // clean up
	if(!mf.empty() && mf.size() >= rstart) { 
		sz i = random_sz(0, mf.size()-1, generator);
		return &(mf[i].c);
	}
	return NULL;
}

bool determine_iu(const conf& c, const output_container& mf, std::vector<bool>& internal_too_close, fl exclusion) { // clean up
	bool tmp = true;
	VINA_FOR_IN(i, mf) {
		bool t = c.internal_too_close(mf[i].c, exclusion);
		internal_too_close[i] = t;
		if(t) 
			tmp = false;
	}
	return tmp;
}

bool conf_is_legal(const conf& c, const output_container& mf, const std::vector<bool>& internal_too_close, const scale& exclusion) {
	assert(mf.size() == internal_too_close.size());
	VINA_FOR_IN(i, mf)
		if(internal_too_close[i] && c.external_too_close(mf[i].c, exclusion))
			return false;
	return true;
}

conf generate_external(const conf& internal_conf, const output_container& mf, const std::vector<bool>& internal_too_close, bool uniq, const scale& spread, const scale& exclusion, fl rp, const conf* rs, unsigned num_attempts, bool& failed, rng& generator) {
	// FIXME if one ligand with no side chains, don't use the same (pos, orient) more than once 
	failed = false;
	VINA_U_FOR(attempt, num_attempts) {
		conf tmp(internal_conf);
		tmp.generate_external(spread, rp, rs, generator); // CHECKME
		if(uniq || conf_is_legal(tmp, mf, internal_too_close, exclusion))
			return tmp;
		if(attempt + 1 >= num_attempts) {
			failed = true;
			return tmp;
		}
	}
	VINA_CHECK(false); 
	return internal_conf; // placating the compiler
}

fl extrapolate_cap(fl from, fl to, fl progress) {
	if(from < epsilon_fl) return from;
	return from * std::exp(progress * std::log(to/from));
}

vec extrapolate_cap(const vec& from, const vec& to, fl progress) {
	vec tmp;
	VINA_FOR_IN(i, tmp)
		tmp[i] = extrapolate_cap(from[i], to[i], progress);
	return tmp;
}

void manifold_phase(sz phase, output_type& out, model& m, const output_container& mf, const precalculate& p, const igrid& ig, fl corner2corner, fl init_manifold_factor, recent_history& e_internal_stats, const manifold& par, rng& generator) { // out.first is starting conf on input
	out.e = max_fl; // FIXME ? back to above?

	const fl rp = 1 - std::exp(std::log(1-par.max_prob) * phase/par.num_phases); // max_prob had better be < 1

	fl e_internal_cost = par.relative_pair_cost * m.num_internal_pairs();
	fl e_external_cost = par.relative_pair_cost * m.num_other_pairs() + m.num_movable_atoms();
	const sz rigid_trials = 1 +  fl_to_sz(1 + par.cost_factor * e_internal_cost / e_external_cost, 100);
	unsigned num_steps = (std::max)(unsigned(1), par.num_steps);

	VINA_U_FOR(step, num_steps) {
		const fl manifold_factor = init_manifold_factor * std::exp(- par.manifold_lambda * step / par.num_steps);
		vec hunt_cap; hunt_cap = extrapolate_cap(par.hunt_cap, authentic_v, fl(step) / num_steps);

		const scale spread(corner2corner * manifold_factor,
			               pi            * manifold_factor,
						   pi            * manifold_factor);

		sz rstart = fl_to_sz(par.rstart_fraction * par.num_phases, par.num_phases);
		const conf* rs = select_rs(mf, rstart, generator);

		output_type candidate = out;
		candidate.c.generate_internal(spread.torsion, rp, rs, generator); // CHECKME

		std::vector<bool> internal_too_close(mf.size());
		bool uniq = determine_iu(candidate.c, mf, internal_too_close, par.exclusion.torsion); // internal_too_close is output

		m.seti(candidate.c);

		VINA_CHECK(rigid_trials > 0);
		bool candidate_is_set = false;
		VINA_FOR(i, rigid_trials) {
			bool failed = false;
			output_type candidate2(generate_external(candidate.c, mf, internal_too_close, uniq, spread, par.exclusion, rp, rs, par.num_attempts, failed, generator), 0);
			if(failed && 
				(candidate_is_set || i+1 < rigid_trials)) // candidate is set or a chance to set it remains
				continue;
			m.sete(candidate2.c);
			candidate2.e = m.evale(p, ig, hunt_cap);
			if(!candidate_is_set || candidate2.e < candidate.e) {
				candidate_is_set = true;
				candidate = candidate2; // this works because internal conf is the same in both // FIXME inefficient, all this
			}
		}
		VINA_CHECK(candidate_is_set);
		if(true) {
			fl e_internal = (3 * step < 2 * num_steps) ? 0 : m.evali(p, hunt_cap);
			e_internal_stats.add(e_internal);
			candidate.e += e_internal;
			if(step == 0 || candidate.e < out.e)
				out = candidate;
		}
	}
}

void manifold::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const {
	manifold tuning_par;
	tuning_par.hunt_cap = authentic_v;

	fl corner2corner = 0.5 * (corner2 - corner1).norm();
	conf_size s = m.get_size();
	change g(s);

	sz np = (std::max)(sz(1), num_phases);
	recent_history e_internal_stats(0, 10, 10);


	VINA_FOR(phase, np) {
		output_type tmp(s, 0);
		tmp.c.randomize(corner1, corner2, generator);
		manifold_phase(phase, tmp, m, out, p_widened, ig_widened, corner2corner, 1, e_internal_stats, *this, generator); // 1 - initial manifold
		if(use_ssd)
			ssd_par(m, p, ig, tmp, g, authentic_v);
		else {
			quasi_newton quasi_newton_par; quasi_newton_par.max_steps = ssd_par.evals;
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);
		}
		m.set(tmp.c);
		tmp.coords = m.get_heavy_atom_movable_coords();
		VINA_CHECK(tmp.coords.size() > 0);
		std::pair<sz, fl> closest_rmsd = find_closest(tmp.coords, out);
		if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
			if(tmp.e < out[closest_rmsd.first].e) { // the new one is better, apparently
				out[closest_rmsd.first] = tmp; // FIXME? slow
			}
		}
		else { // nothing similar
			out.push_back(new output_type(tmp));
		}
	}
	// final tunings
	const output_container null_array;
	VINA_CHECK(!out.empty());
	out.sort();
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
