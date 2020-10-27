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
    ScoringFunction(const scoring_function_choice sf_choice, const flv& weights){
        switch (sf_choice)
        {
            case SF_VINA:
            {
                m_potentials.push_back(new vina_gaussian(0, 0.5, 8.0));
                m_potentials.push_back(new vina_gaussian(3, 2.0, 8.0));
                m_potentials.push_back(new vina_repulsion(0.0, 8.0));
                m_potentials.push_back(new vina_hydrophobic(0.5, 1.5, 8.0));
                m_potentials.push_back(new vina_non_dir_h_bond(-0.7, 0, 8.0));
                m_potentials.push_back(new linearattraction(20.0));
                m_conf_independents.push_back(new num_tors_div());
                m_atom_typing = atom_type::XS;
                m_cutoff = 8.0;
                m_max_cutoff = 20.0;
                break;
            }
            case SF_VINARDO:
            {
                m_potentials.push_back(new vinardo_gaussian(0, 0.8, 8.0));
                m_potentials.push_back(new vinardo_repulsion(0, 8.0));
                m_potentials.push_back(new vinardo_hydrophobic(0, 2.5, 8.0));
                m_potentials.push_back(new vinardo_non_dir_h_bond(-0.6, 0, 8.0));
                m_potentials.push_back(new linearattraction(20.0));
                m_conf_independents.push_back(new num_tors_div());
                m_atom_typing = atom_type::XS;
                m_cutoff = 8.0;
                m_max_cutoff = 20.0;
                break;
            }
            case SF_AD42:
            {
                m_potentials.push_back(new ad4_vdw(0.5, 100000, 8.0));
                m_potentials.push_back(new ad4_hb(0.5, 100000, 8.0));
                m_potentials.push_back(new ad4_electrostatic(100, 20.48));
                m_potentials.push_back(new ad4_solvation(3.6, 0.01097, true, 20.48));
                m_potentials.push_back(new linearattraction(20.0));
                m_conf_independents.push_back(new ad4_tors_add());
                m_atom_typing = atom_type::AD;
                m_cutoff = 20.48;
                m_max_cutoff = 20.48;
                break;
            }
            default:
            {
                std::cout << "INSIDE everything::everything()   sfchoice = " << sf_choice << "\n";
                VINA_CHECK(false);
                break;
            }
        }

        m_num_potentials = m_potentials.size();
        m_num_conf_independents = m_conf_independents.size();
        m_weights = weights;
    };
    ~ScoringFunction() { }
    fl eval(atom& a, atom& b, fl r) const{ // intentionally not checking for cutoff
        fl acc = 0;
        VINA_FOR (i, m_num_potentials)
        {
            acc += m_weights[i] * m_potentials[i]->eval(a, b, r);
        }
        return acc;
    };
    fl eval(sz t1, sz t2, fl r) const{
        fl acc = 0;
        VINA_FOR (i, m_num_potentials)
        {
            acc += m_weights[i] * m_potentials[i]->eval(t1, t2, r);
        }
        return acc;
    };
    fl conf_independent(const model& m, fl e) const{
        // Iterator for weights
        flv::const_iterator it = m_weights.begin() + m_num_potentials;
        conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?
        VINA_FOR (i, m_num_conf_independents)
        {
            // We don't accumulate energy. Why? I don't know...
            e = m_conf_independents[i]->eval(in, e, it);
        }
        assert(it == m_weights.end());
        return e;
    };
    fl get_cutoff() const { return m_cutoff; }
    fl get_max_cutoff() const { return m_max_cutoff; }
    atom_type::t get_atom_typing() const { return m_atom_typing; }
    szv get_atom_types() const{
        szv tmp;
        VINA_FOR(i, num_atom_types(m_atom_typing))
        {
          tmp.push_back(i);
        }
        return tmp;
    };
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
