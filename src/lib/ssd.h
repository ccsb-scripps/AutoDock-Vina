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

#ifndef VINA_SSD_H
#define VINA_SSD_H

#include "model.h"

struct ssd {
	unsigned evals;
	fl initial_factor;
	fl min_factor;
	fl up;
	fl down;
	void print() const { std::cout << "evals=" << evals << ", initial_factor=" << initial_factor << ", min_factor=" << min_factor << ", up=" << up << ", down=" << down; }
	ssd() : evals(300), initial_factor(1e-4), min_factor(1e-6), up(1.6), down(0.5) {}
	// clean up
	void operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v) const; // g must have correct size
};

#endif
