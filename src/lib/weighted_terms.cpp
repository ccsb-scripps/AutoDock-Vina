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

#include "weighted_terms.h"

weighted_terms::weighted_terms(const terms* t, const flv& weights) : t(t), weights(weights), cutoff_(0) { // does not own t
	VINA_CHECK(t->distance_additive_terms.num_enabled() == 0);
	VINA_CHECK(t->         additive_terms.num_enabled() == 0);
	VINA_CHECK(t->   intermolecular_terms.num_enabled() == 0);
	VINA_FOR_IN(i, t->usable_terms)
		if(t->usable_terms.enabled[i]) {
			if(enabled_usable_terms.empty())
				atom_typing_used_ = t->usable_terms[i].atom_typing_used;
			else
				VINA_CHECK(atom_typing_used_ == t->usable_terms[i].atom_typing_used);

			enabled_usable_terms.push_back(i);
			cutoff_ = (std::max)(cutoff_, t->usable_terms[i].cutoff);
		}
}

fl weighted_terms::eval(sz t1, sz t2, fl r) const { // intentionally not checking for cutoff
	fl acc = 0;
	VINA_FOR_IN(i, enabled_usable_terms) 
		acc += weights[i] * t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r);
	return acc;
}
fl weighted_terms::conf_independent(const model& m, fl e) const {
	flv::const_iterator it = weights.begin() + enabled_usable_terms.size();
	conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?
	fl tmp = t->eval_conf_independent(in, e, it);
	assert(it == weights.end());
	return tmp;
}
