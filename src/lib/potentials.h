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

// Vina common functions
inline fl slope_step(fl x_bad, fl x_good, fl x) {
    if (x_bad < x_good) {
        if (x <= x_bad) return 0;
        if (x >= x_good) return 1;
    }
    else {
        if (x >= x_bad) return 0;
        if (x <= x_good) return 1;
    }
    return (x - x_bad) / (x_good - x_bad);
}

inline bool is_glue_type(sz xs_t) {
    if ((xs_t==XS_TYPE_G0) || (xs_t==XS_TYPE_G1) || (xs_t==XS_TYPE_G2) || (xs_t==XS_TYPE_G3)) return true;
    return false;}

inline fl optimal_distance(sz xs_t1, sz xs_t2) {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) return 0.0; // G0, G1, G2 or G3
    return xs_radius(xs_t1) + xs_radius(xs_t2);
}

inline fl smooth_div(fl x, fl y) {
    if (std::abs(x) < epsilon_fl) return 0;
    if (std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
    return x / y;
}

// Vinardo common functions
inline fl optimal_distance_vinardo(sz xs_t1, sz xs_t2) {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) return 0.0; // G0, G1, G2 or G3
    return xs_vinardo_radius(xs_t1) + xs_vinardo_radius(xs_t2);
}

// AD42 common functions
inline fl smoothen(fl r, fl rij, fl smoothing) {
    fl out;
    smoothing *= 0.5;

    if (r > rij + smoothing)
        out = r - smoothing;
    else if(r < rij - smoothing)
        out = r + smoothing;
    else
        out = rij;

    return out;
}

inline fl ad4_hb_eps(sz& a) {
    if (a < AD_TYPE_SIZE) return ad_type_property(a).hb_depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

inline fl ad4_hb_radius(sz& t) {
    if (t < AD_TYPE_SIZE) return ad_type_property(t).hb_radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

inline fl ad4_vdw_eps(sz& a) {
    if(a < AD_TYPE_SIZE) return ad_type_property(a).depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

inline fl ad4_vdw_radius(sz& t) {
    if(t < AD_TYPE_SIZE) return ad_type_property(t).radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

// Macrocycle - Vina and AD42
inline bool is_glued(sz xs_t1, sz xs_t2) {
    return (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_H_CG0) ||
       (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_P_CG0) ||
       (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_H_CG0) ||
       (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_P_CG0) ||

       (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_H_CG1) ||
       (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_P_CG1) ||
       (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_H_CG1) ||
       (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_P_CG1) ||

       (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_H_CG2) ||
       (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_P_CG2) ||
       (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_H_CG2) ||
       (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_P_CG2) ||

       (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_H_CG3) ||
       (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_P_CG3) ||
       (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_H_CG3) ||
       (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_P_CG3);
}

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
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        return gauss(r - (optimal_distance(a.xs, b.xs) + offset)); // hard-coded to XS atom type
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        return gauss(r - (optimal_distance(t1, t2) + offset)); // hard-coded to XS atom type
    };
    fl get_cutoff() { return cutoff; }
private:
    fl offset; // added to optimal distance
    fl width;
    fl cutoff;

    fl gauss(fl x) {
        return std::exp(-sqr(x / width));
    };
};

class vina_repulsion : public Potential {
public:
    vina_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) { }
    //~vina_repulsion() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        fl d = r - (optimal_distance(a.xs, b.xs) + offset); // hard-coded to XS atom type
        if (d > 0.0)
            return 0.0;
        return d * d;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl d = r - (optimal_distance(t1, t2) + offset); // hard-coded to XS atom type
        if(d > 0.0)
            return 0.0;
        return d*d;
    };
    fl get_cutoff() { return cutoff; }
private:
    fl offset; // added to vdw
    fl cutoff;
};

class vina_hydrophobic : public Potential {
public:
    vina_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }
    //~vina_hydrophobic() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
        else return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        else
            return 0.0;
    };
    fl get_cutoff() { return cutoff; }
private:
    fl good;
    fl bad;
    fl cutoff;
};

class vina_non_dir_h_bond : public Potential {
public:
    vina_non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }
    //~vina_non_dir_h_bond() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_h_bond_possible(a.xs, b.xs))
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
        return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_h_bond_possible(t1, t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        return 0.0;
    };
    fl get_cutoff() { return cutoff; }
private:
    fl good;
    fl bad;
    fl cutoff;
};

// Vinardo
class vinardo_gaussian : public Potential {
public:
    vinardo_gaussian(fl offset_, fl width_, fl cutoff_) : offset(offset_), width(width_), cutoff(cutoff_) { }
    //~vina_gaussian() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        return gauss(r - (optimal_distance_vinardo(a.xs, b.xs) + offset)); // hard-coded to XS atom type
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        return gauss(r - (optimal_distance_vinardo(t1, t2) + offset)); // hard-coded to XS atom type
    };
    fl get_cutoff() { return cutoff; }
private:
    fl offset; // added to optimal distance
    fl width;
    fl cutoff;

    fl gauss(fl x) {
        return std::exp(-sqr(x / width));
    };
};

class vinardo_repulsion : public Potential {
public:
    vinardo_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) {}
    //~vina_repulsion() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        fl d = r - (optimal_distance_vinardo(a.xs, b.xs) + offset); // hard-coded to XS atom type
        if (d > 0.0)
            return 0.0;
        return d * d;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl d = r - (optimal_distance_vinardo(t1, t2) + offset); // hard-coded to XS atom type
        if(d > 0.0)
            return 0.0;
        return d*d;
    };
    fl get_cutoff() { return cutoff; }
private:
    fl offset; // added to vdw
    fl cutoff;
};

class vinardo_hydrophobic : public Potential {
public:
    vinardo_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }
    //~vina_hydrophobic() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs));
        else return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2));
        else
            return 0.0;
    };
    fl get_cutoff() { return cutoff; }
private:
    fl good;
    fl bad;
    fl cutoff;
};

class vinardo_non_dir_h_bond : public Potential {
public:
    vinardo_non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }
    //~vina_non_dir_h_bond() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_h_bond_possible(a.xs, b.xs))
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs));
        return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_h_bond_possible(t1, t2))
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2));
        return 0.0;
    };
    fl get_cutoff() { return cutoff; }
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
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl q1q2 = a.charge * b.charge * 332.0;
        fl B = 78.4 + 8.5525;
        fl lB = -B * 0.003627;
        fl diel = -8.5525 + (B / (1 + 7.7839 * std::exp(lB * r)));
        if (r < epsilon_fl)
            return q1q2 * cap / diel;
        else {
            return q1q2 * (std::min)(cap, 1.0 / (r * diel));
        }
    };
    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }
private:
    fl cap;
    fl cutoff;
};

class ad4_solvation : public Potential {
public:
    ad4_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_), cutoff(cutoff_) { }
    //~ad4_solvation() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl q1 = a.charge;
        fl q2 = b.charge;
        VINA_CHECK(not_max(q1));
        VINA_CHECK(not_max(q2));
        fl solv1 = solvation_parameter(a);
        fl solv2 = solvation_parameter(b);
        fl volume1 = volume(a);
        fl volume2 = volume(b);
        fl my_solv = charge_dependent ? solvation_q : 0;
        fl tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 + 
                  (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-0.5 * sqr(r / desolvation_sigma));
        VINA_CHECK(not_max(tmp));
        return tmp;
    };
    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }
private:
    fl desolvation_sigma;
    fl solvation_q;
    bool charge_dependent;
    fl cutoff;

    fl volume(const atom_type& a) const {
        if (a.ad < AD_TYPE_SIZE)
            return ad_type_property(a.ad).volume;
        else if (a.xs < XS_TYPE_SIZE)
            return 4.0 * pi / 3.0 * int_pow<3>(xs_radius(a.xs));
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };
    fl solvation_parameter(const atom_type& a) const {
        if (a.ad < AD_TYPE_SIZE)
            return ad_type_property(a.ad).solvation;
        else if (a.xs == XS_TYPE_Met_D)
            return metal_solvation_parameter;
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };
};

class ad4_vdw : public Potential {
public:
    ad4_vdw(fl smoothing_, fl cap_, fl cutoff_) : smoothing(smoothing_), cap(cap_), cutoff(cutoff_) { }
    //~ad4_vdw() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        sz t1 = a.ad;
        sz t2 = b.ad;
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        fl vdw_depth = std::sqrt(ad4_vdw_eps(t1) * ad4_vdw_eps(t2));
        if (hb_depth < 0) return 0.0; // interaction is hb, not vdw.
        r = smoothen(r, vdw_rij, smoothing);
        fl c_12 = int_pow<12>(vdw_rij) * vdw_depth;
        fl c_6  = int_pow<6>(vdw_rij)  * vdw_depth * 2.0;
        fl r6   = int_pow<6>(r);
        fl r12  = int_pow<12>(r);
        if(r12 > epsilon_fl && r6 > epsilon_fl)
            return (std::min)(cap, c_12 / r12 - c_6 / r6);
        else
            return cap;
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };
    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }
private:
    fl smoothing;
    fl cap;
    fl cutoff;
};

class ad4_hb : public Potential {
public:
    ad4_hb(fl smoothing_, fl cap_, fl cutoff_) : smoothing(smoothing_), cap(cap_), cutoff(cutoff_) { }
    //~ad4_hb() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        sz t1 = a.ad;
        sz t2 = b.ad;
        fl hb_rij = ad4_hb_radius(t1) + ad4_hb_radius(t2);
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        if (hb_depth >= 0)
            return 0.0; // interaction is vdw, not hb.
        r = smoothen(r, hb_rij, smoothing);
        fl c_12 = int_pow<12>(hb_rij) * -hb_depth * 10 / 2.0;
        fl c_10 = int_pow<10>(hb_rij) * -hb_depth * 12 / 2.0;
        fl r10  = int_pow<10>(r);
        fl r12  = int_pow<12>(r);
        if (r12 > epsilon_fl && r10 > epsilon_fl)
            return (std::min)(cap, c_12 / r12 - c_10 / r10);
        else
            return cap;
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };
    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }
private:
    fl smoothing;
    fl cap;
    fl cutoff;
};

// Macrocycle - Vina and AD42
class linearattraction : public Potential {
public:
    linearattraction(fl cutoff_): cutoff(cutoff_) { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if (is_glued(a.xs, b.xs))
            return r;
        else 
            return 0.0;
    };
    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if (is_glued(t1, t2))
            return r;
        else
            return 0.0;
    };
    fl get_cutoff() { return cutoff; }
private:
    fl cutoff;
};

#endif
