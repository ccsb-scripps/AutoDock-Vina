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

#ifndef VINA_POTENTIALS_H
#define VINA_POTENTIALS_H

#include "atom.h"
#include "int_pow.h"


class Potential {
public:
    virtual ~Potential() { }
    virtual fl eval(const atom& a, const atom& b, fl r) { return 0; };
    virtual fl eval(sz t1, sz t2, fl r) { return 0; };
    virtual fl get_cutoff() { return 0; }
};

// Vina
class vina_gaussian : public Potential {
public:
    vina_gaussian(fl offset_, fl width_, fl cutoff_) : offset(offset_), width(width_), cutoff(cutoff_) { }
    //~vina_gaussian() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override;
    fl get_cutoff() override { return cutoff; }
private:
    fl offset; // added to optimal distance
    fl width;
    fl cutoff;
    
    fl gauss(fl x) const;
};

class vina_repulsion : public Potential {
public:
    vina_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) { }
    //~vina_repulsion() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override;
    fl get_cutoff() override { return cutoff; }
private:
    fl offset; // added to vdw
    fl cutoff;
};

class vina_hydrophobic : public Potential {
public:
    vina_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }
    //~vina_hydrophobic() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override;
    fl get_cutoff() override { return cutoff; }
private:
    fl good;
    fl bad;
    fl cutoff;
};

class vina_non_dir_h_bond : public Potential {
public:
    vina_non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }
    //~vina_non_dir_h_bond() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override;
    fl get_cutoff() override { return cutoff; }
private:
    fl good;
    fl bad;
    fl cutoff;
};

// AD42
class ad4_electrostatic : public Potential {
public:
    ad4_electrostatic(fl cap_, fl cutoff_) : cap(cap_), cutoff(cutoff_) { }
    //~ad4_electrostatic() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override { return 0; }
    fl get_cutoff() override { return cutoff; }
private:
    fl cap;
    fl cutoff;
};

class ad4_solvation : public Potential {
public:
    ad4_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_), cutoff(cutoff_) { }
    //~ad4_solvation() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override { return 0; }
    fl get_cutoff() override { return cutoff; }
private:
    fl desolvation_sigma;
    fl solvation_q;
    bool charge_dependent;
    fl cutoff;

    fl volume(const atom_type& a) const;
    fl solvation_parameter(const atom_type& a) const;
};

class ad4_vdw : public Potential {
public:
    ad4_vdw(fl smoothing_, fl cap_, fl cutoff_) : smoothing(smoothing_), cap(cap_), cutoff(cutoff_) { }
    //~ad4_vdw() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override { return 0; }
    fl get_cutoff() override { return cutoff; }
private:
    fl smoothing;
    fl cap;
    fl cutoff;
};

class ad4_hb : public Potential {
public:
    ad4_hb(fl smoothing_, fl cap_, fl cutoff_) : smoothing(smoothing_), cap(cap_), cutoff(cutoff_) { }
    //~ad4_hb() { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override { return 0; }
    fl get_cutoff() override { return cutoff; }
private:
    fl smoothing;
    fl cap;
    fl cutoff;
};

// Macrocycle - Vina and AD42
class linearattraction : public Potential {
public:
    linearattraction(fl cutoff_): cutoff(cutoff_) { }
    fl eval(const atom& a, const atom& b, fl r) override;
    fl eval(sz t1, sz t2, fl r) override;
    fl get_cutoff() override { return cutoff; }
private:
    fl cutoff;
};

#endif