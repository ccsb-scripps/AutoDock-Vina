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

#ifndef VINA_CONF_INDEPENDENT_H
#define VINA_CONF_INDEPENDENT_H

#include <stdlib.h>
#include "atom.h"
#include "common.h"


// Forward declaration
struct model;

class conf_independent_inputs {
public:
    fl torsdof; // from TORSDOF keyword in pdbqt file
    fl num_tors;
    fl num_rotors;
    fl num_heavy_atoms;
    fl num_hydrophobic_atoms;
    fl ligand_max_num_h_bonds;
    fl num_ligands;
    fl ligand_lengths_sum;

    operator flv() const;
    conf_independent_inputs();
    conf_independent_inputs(const model& m);
    std::vector<std::string> get_names() const;
private:
    unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i) const; // FIXME? - could be static, but I don't feel like declaring function friends
    unsigned atom_rotors(const model& m, const atom_index& i) const; // the number of rotatable bonds to heavy ligand atoms
};

// Conf independent
class ConfIndependent {
public:
    virtual ~ConfIndependent() { }
    virtual fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) { return 0; };
};

// Vina
class num_tors_sqr : public ConfIndependent {
public:
    num_tors_sqr() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class num_tors_sqrt : public ConfIndependent {
public:
    num_tors_sqrt() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class num_tors_div : public ConfIndependent {
public:
    num_tors_div() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class ligand_length : public ConfIndependent {
public:
    ligand_length() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class num_ligands : public ConfIndependent {
public:
    num_ligands() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class num_heavy_atoms_div : public ConfIndependent {
public:
    num_heavy_atoms_div() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class num_heavy_atoms : public ConfIndependent {
public:
    num_heavy_atoms() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

class num_hydrophobic_atoms : public ConfIndependent {
public:
    num_hydrophobic_atoms() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

// AD42
class ad4_tors_add : public ConfIndependent {
public:
    ad4_tors_add() { }
    sz size() const { return 1; }
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

#endif
