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

#include "quasi_newton.h"
#include "bfgs.h"


struct quasi_newton_aux {
    model* m;
    const precalculate_byatom* p;
    const igrid* ig;
    const vec v;

    quasi_newton_aux(model* m_, const precalculate_byatom* p_, const igrid* ig_, const vec& v_) : m(m_), p(p_), ig(ig_), v(v_) {}
    
    fl operator()(const conf& c, change& g) {
        // Before evaluating conf, we have to update model
        m->set(c);
        const fl tmp = m->eval_deriv(*p, *ig, v, g);
        return tmp;
    }
};

void quasi_newton::operator()(model& m, const precalculate_byatom& p, const igrid& ig, output_type& out, change& g, const vec& v, int& evalcount) const { // g must have correct size
    quasi_newton_aux aux(&m, &p, &ig, v);

    fl res = bfgs(aux, out.c, g, max_steps, average_required_improvement, 10, evalcount);

    // Update model a last time after optimization
    m.set(out.c);
    out.e = res;
}
