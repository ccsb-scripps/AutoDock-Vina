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
    ScoringFunction(const scoring_function_choice sf_choice, const flv& weights, const atom_type::t atom_types) {
        switch(sf_choice) {
            case SF_VINA: {
                m_cutoff = 8.0;
                m_potentials = {new vina_gaussian(0, 0.5),
                                new vina_gaussian(3, 2.0),
                                new vina_repulsion(0.0),
                                new vina_hydrophobic(0.5, 1.5),
                                new vina_non_dir_h_bond(-0.7, 0)};
                m_conf_independents = {new num_tors_div()};
                break;
            }
            case SF_VINARDO: {
                std::cout << "\n\nVinardo scoring function is not implemented yet.\n\nAborting.\n\n";
                VINA_CHECK(false);
                break;
            }
            case SF_AD42: {
                m_cutoff = 20.48;
                m_potentials = {new ad4_vdw(0.5, 100000),
                                new ad4_hb(0.5, 100000),
                                new ad4_electrostatic(100),
                                new ad4_solvation(3.6, 0.01097, true)};
                m_conf_independents = {new ad4_tors_add()};
                break;
            }
            default: {
                std::cout << "INSIDE everything::everything()   sfchoice = " << sf_choice << "\n";
                VINA_CHECK(false);
                break;
            }
        }

        m_num_potentials = m_potentials.size();
        m_num_conf_independents = m_conf_independents.size();
        m_weights = weights;
        m_atom_types = atom_types;
    }

    ~ScoringFunction() { }

    fl eval(atom& a, atom& b, fl r) const; // intentionally not checking for cutoff
    fl eval(sz t1, sz t2, fl r) const;
    fl conf_independent(const model& m, fl e) const;
    fl get_cutoff() const { return m_cutoff; }
    atom_type::t get_atom_types() const { return m_atom_types; }

private:
    std::vector<Potential*> m_potentials;
    std::vector<ConfIndependent*> m_conf_independents;
    flv m_weights;
    fl m_cutoff;
    int m_num_potentials;
    int m_num_conf_independents;
    atom_type::t m_atom_types;
};

#endif
