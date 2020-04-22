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

#ifndef VINA_MANIFOLD_H
#define VINA_MANIFOLD_H

// This is not used  by Vina proper any more, but "design" sort of needs it

#include "ssd.h"

struct manifold {
	ssd ssd_par;
	sz num_phases;
	unsigned num_steps;
	unsigned num_attempts;
	sz num_final_tunings;

	//fl manifold_mid_factor;
	fl manifold_lambda;
	fl cost_factor;
	fl max_prob;
	fl rstart_fraction;
	fl min_rmsd;

	vec hunt_cap;

	scale exclusion;

	fl relative_pair_cost; // the cost of calculating  pairwise energy / interaction of atom with grid

	bool use_ssd;

	void print() const { 
		std::cout << "{ssd: "; ssd_par.print(); 
		std::cout << "}\n" << "num_phases=" << num_phases 
			              << ",\nnum_steps=" << num_steps 
						  << ",\nnum_attempts=" << num_attempts 
						  << ",\nnum_final_tunings=" << num_final_tunings 
						  << ",\nmanifold_lambda=" << manifold_lambda
						  << ",\ncost_factor=" << cost_factor
						  << ",\nmax_prob=" << max_prob
						  << ",\nrstart_fraction=" << rstart_fraction
						  << ",\nmin_rmsd=" << min_rmsd
						  << ",\nhunt_cap=" << hunt_cap[0] << " " << hunt_cap[1] << " " << hunt_cap[2]
						  << ",\nexclusion=" << exclusion.position << " " << exclusion.orientation << " " << exclusion.torsion
					      << ",\nrelative_pair_cost=" << relative_pair_cost
						  << ",\nuse_ssd=" << use_ssd << std::endl;
	}
	manifold() 
		: num_phases(800), num_steps(50), num_attempts(10), num_final_tunings(1),
		  /*manifold_mid_factor(0.1), */ manifold_lambda(0), cost_factor(0.5), max_prob(0.9), rstart_fraction(0.25), min_rmsd(2),
		  //hunt_cap(2, 0.1, 2),
		  hunt_cap(0.01, 0.001, 0.01),
		  exclusion(1, pi/24, pi/12),
	      relative_pair_cost(0.5),
	      use_ssd(false) {}

	output_type operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const;

	// out is sorted
	void operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const;

};

#endif
