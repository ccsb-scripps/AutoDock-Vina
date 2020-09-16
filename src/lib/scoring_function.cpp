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

#include "scoring_function.h"
#include "model.h"
#include "atom_type.h"


fl ScoringFunction::eval(atom& a, atom& b, fl r) const { // intentionally not checking for cutoff
    fl acc = 0;

    VINA_FOR (i, m_num_potentials) {
        acc += m_weights[i] * m_potentials[i]->eval(a, b, r);
    }

    return acc;
}

fl ScoringFunction::eval(sz t1, sz t2, fl r) const { // intentionally not checking for cutoff
    fl acc = 0;

    VINA_FOR (i, m_num_potentials) {
        acc += m_weights[i] * m_potentials[i]->eval(t1, t2, r);
    }

    return acc;
}

fl ScoringFunction::conf_independent(const model& m, fl e) const {
    // Iterator for weights
    flv::const_iterator it = m_weights.begin() + m_num_potentials;
    conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?

    VINA_FOR (i, m_num_conf_independents) {
        // We don't accumulate energy. Why? I don't know...
        e = m_conf_independents[i]->eval(in, e, it);
    }

    assert(it == m_weights.end());
    return e;
}

szv ScoringFunction::get_atom_types() const {
    szv tmp;

    VINA_FOR(i, num_atom_types(m_atom_typing)) {
      tmp.push_back(i);
    }

    return tmp;
}