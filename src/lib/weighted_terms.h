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

#ifndef VINA_WEIGHTED_TERMS_H
#define VINA_WEIGHTED_TERMS_H

#include "terms.h"

struct weighted_terms : public scoring_function {
	weighted_terms(const terms* t, const flv& weights); // does not own t
	atom_type::t atom_typing_used() const { return atom_typing_used_; }
	fl cutoff() const { return cutoff_; }
	fl eval(sz t1, sz t2, fl r) const; // intentionally not checking for cutoff
	fl conf_independent(const model& m, fl e) const;
private:
	weighted_terms() {}
	const terms* t;
	flv weights;
	atom_type::t atom_typing_used_;
	fl cutoff_;
	szv enabled_usable_terms;
};

#endif
