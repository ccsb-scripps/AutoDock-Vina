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

#ifndef VINA_QUASI_NEWTON_H
#define VINA_QUASI_NEWTON_H

#include "model.h"

struct quasi_newton {
    unsigned max_steps;
    fl average_required_improvement;
    quasi_newton() : max_steps(1000), average_required_improvement(0.0) {}
    // clean up
    void operator()(model& m, const precalculate_byatom& p, const igrid& ig, output_type& out, change& g, const vec& v, int& evalcount) const; // g must have correct size
};

#endif
