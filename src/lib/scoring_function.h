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

#ifndef VINA_SCORING_FUNCTION_H
#define VINA_SCORING_FUNCTION_H

#include <stdlib.h>
#include <list>
#include "atom.h"
#include "conf_independent.h"
#include "potentials.h"
#include "common.h"


//Forward declaration
struct model;

enum scoring_function_choice {SF_VINA, SF_AD42, SF_VINARDO};

class ScoringFunction {
public:
    ScoringFunction() { }
    ScoringFunction(const scoring_function_choice sf_choice, const flv& weights);
    ~ScoringFunction() { }
    fl eval(atom& a, atom& b, fl r) const; // intentionally not checking for cutoff
    fl eval(sz t1, sz t2, fl r) const;
    fl conf_independent(const model& m, fl e) const;
    fl get_cutoff() const { return m_cutoff; }
    fl get_max_cutoff() const { return m_max_cutoff; }
    atom_type::t get_atom_typing() const { return m_atom_typing; }
    szv get_atom_types() const;
    sz get_num_atom_types() const { return num_atom_types(m_atom_typing); }
    flv get_weights() const { return m_weights; }

private:
    std::vector<Potential*> m_potentials;
    std::vector<ConfIndependent*> m_conf_independents;
    flv m_weights;
    fl m_cutoff;
    fl m_max_cutoff;
    int m_num_potentials;
    int m_num_conf_independents;
    atom_type::t m_atom_typing;
};

#endif
