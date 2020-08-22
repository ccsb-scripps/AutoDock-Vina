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

#include "everything.h"
#include "int_pow.h"

inline fl gaussian(fl x, fl width) {
	return std::exp(-sqr(x/width));
}

// distance_additive terms

template<unsigned i>
struct electrostatic : public usable {
	fl cap;
	electrostatic(fl cap_, fl cutoff_) : usable(cutoff_), cap(cap_) {
		name = std::string("electrostatic(i=") + to_string(i) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
		fl tmp = int_pow<i>(r);
		fl q1q2 = a.charge * b.charge * 332.;
        fl B = 78.4 + 8.5525;
        fl lB = -B * 0.003627;
        fl diel = -8.5525 + (B / (1 + 7.7839 * std::exp(lB * r)));
		if(tmp < epsilon_fl) return q1q2 * cap / diel;
		else {
            //std::cout << diel << " dist=" << r << " q1=" << a.charge << " q2=" << b.charge <<  " q1q2=" << q1q2 << " tmp=" << tmp << " energy=" << q1q2 * (std::min)(cap, 1/(int_pow<i>(r) * diel)) << "\n";
            return q1q2 * (std::min)(cap, 1/(int_pow<i>(r) * diel));
        }
	}
    fl eval(sz t1, sz t2, fl r) const {return 0;}
};

fl solvation_parameter(const atom_type& a) {
	if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).solvation;
	else if(a.xs == XS_TYPE_Met_D) return metal_solvation_parameter;
	VINA_CHECK(false); 
	return 0; // placating the compiler
}

fl volume(const atom_type& a) {
	if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).volume;
	else if(a.xs < XS_TYPE_SIZE) return 4*pi / 3 * int_pow<3>(xs_radius(a.xs));
	VINA_CHECK(false);
	return 0; // placating the compiler
}


struct ad4_solvation : public usable {
	fl desolvation_sigma;
	fl solvation_q;
	bool charge_dependent;
	ad4_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : usable(cutoff_), solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_) {
		name = std::string("ad4_solvation(d-sigma=") + to_string(desolvation_sigma) + ", s/q=" + to_string(solvation_q) + ", q=" + to_string(charge_dependent) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
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
			    (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-0.5 * sqr(r/(desolvation_sigma)));

		VINA_CHECK(not_max(tmp));
		return tmp;
	}
    fl eval(sz t1, sz t2, fl r) const {return 0;}

};

inline fl optimal_distance(sz xs_t1, sz xs_t2) {
	return xs_radius(xs_t1) + xs_radius(xs_t2);
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

inline fl ad4_hb_eps(sz& a) {
    if(a < AD_TYPE_SIZE) return ad_type_property(a).hb_depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

inline fl ad4_hb_radius(sz& t) {
    if(t < AD_TYPE_SIZE) return ad_type_property(t).hb_radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}
struct gauss : public usable {
    fl optimal_r;
	fl offset; // added to optimal distance
	fl width;
	gauss(fl offset_, fl width_, fl cutoff_) : usable(cutoff_), offset(offset_), width(width_) {
		name = std::string("gauss(o=") + to_string(offset) + ", w=" + to_string(width) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
		return gaussian(r - (optimal_distance(a.xs, b.xs) + offset), width); // hard-coded to XS atom type
	}
	fl eval(sz t1, sz t2, fl r) const {
		return gaussian(r - (optimal_distance(t1, t2) + offset), width); // hard-coded to XS atom type
	}
};

struct repulsion : public usable {
	fl offset; // added to vdw
	repulsion(fl offset_, fl cutoff_) : usable(cutoff_), offset(offset_) {
		name = std::string("repulsion(o=") + to_string(offset) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
		fl d = r - (optimal_distance(a.xs, b.xs) + offset); // hard-coded to XS atom type
		if(d > 0) 
			return 0;
		return d*d;
	}
	fl eval(sz t1, sz t2, fl r) const {
		fl d = r - (optimal_distance(t1, t2) + offset); // hard-coded to XS atom type
		if(d > 0) 
			return 0;
		return d*d;
	}

};

inline fl slope_step(fl x_bad, fl x_good, fl x) {
	if(x_bad < x_good) {
		if(x <= x_bad) return 0;
		if(x >= x_good) return 1;
	}
	else {
		if(x >= x_bad) return 0;
		if(x <= x_good) return 1;
	}
	return (x - x_bad) / (x_good - x_bad);
}

struct hydrophobic : public usable {
	fl good;
	fl bad;
	hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
		name = "hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
		if(xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
			return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
		else return 0;
	}
	fl eval(sz t1, sz t2, fl r) const {
		if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		else return 0;
	}

};

struct non_hydrophobic : public usable {
	fl good;
	fl bad;
	non_hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
		name = "non_hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const {
		if(!xs_is_hydrophobic(t1) && !xs_is_hydrophobic(t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		else return 0;
	}
};

template<unsigned n, unsigned m>
void find_vdw_coefficients(fl position, fl depth, fl& c_n, fl& c_m) {
	BOOST_STATIC_ASSERT(n != m); 
	c_n = int_pow<n>(position) * depth * m / (fl(n)-fl(m));
	c_m = int_pow<m>(position) * depth * n / (fl(m)-fl(n));
}


template<unsigned i, unsigned j>
struct vdw : public usable {
	fl smoothing;
	fl cap;
	vdw(fl smoothing_, fl cap_, fl cutoff_) 
		: usable(cutoff_), smoothing(smoothing_), cap(cap_) {
		name = "vdw(i=" + to_string(i) + ", j=" + to_string(j) + ", s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const {
		fl d0 = optimal_distance(t1, t2);
		fl depth = 1; 
		fl c_i = 0;
		fl c_j = 0;
		find_vdw_coefficients<i, j>(d0, depth, c_i, c_j);
		if     (r > d0 + smoothing) r -= smoothing;
		else if(r < d0 - smoothing) r += smoothing;
		else r = d0;

		fl r_i = int_pow<i>(r);
		fl r_j = int_pow<j>(r);
		if(r_i > epsilon_fl && r_j > epsilon_fl)
			return (std::min)(cap, c_i / r_i + c_j / r_j);
		else 
			return cap;
	}
};

fl smoothen(fl r, fl rij, fl smoothing) {
    fl out;
    smoothing *= 0.5;
    if     (r > rij + smoothing) out = r - smoothing;
    else if(r < rij - smoothing) out = r + smoothing;
    else out = rij;
    return out;
}

struct ad4_vdw : public usable {
	fl smoothing;
	fl cap;
	ad4_vdw(fl smoothing_, fl cap_, fl cutoff_) 
		: usable(cutoff_), smoothing(smoothing_), cap(cap_) {
		name = "ad4_vdw(s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
        //std::cout << "we are here!!!\n";
        sz t1 = a.ad; 
        sz t2 = b.ad; 
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
		fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
		fl vdw_depth = std::sqrt(ad4_vdw_eps(t1) * ad4_vdw_eps(t2));
        
        if (hb_depth < 0) return 0; // interaction is hb, not vdw.

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
        return 0; // placating the compiler
	}
    fl eval(sz t1, sz t2, fl r) const {return 0;}
};

struct ad4_hb : public usable {
	fl smoothing;
	fl cap;
	ad4_hb(fl smoothing_, fl cap_, fl cutoff_) 
		: usable(cutoff_), smoothing(smoothing_), cap(cap_) {
		name = "ad4_hb(s=" + to_string(smoothing) + ", ^=" + to_string(cap) + ", c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
        sz t1 = a.ad; 
        sz t2 = b.ad; 
        //std::cout << "HERE " << t1 << "\n";
        fl hb_rij = ad4_hb_radius(t1) + ad4_hb_radius(t2);
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        //std::cout << "hb_eps " << ad4_hb_eps(t1) << " " << ad4_hb_eps(t2) << " r=" << r << " t1=" << t1 << ", t2= " << t2 << "\n";
		fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        
        if (hb_depth >= 0) return 0; // interaction is vdw, not hb.

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
        return 0; // placating the compiler
	}
    fl eval(sz t1, sz t2, fl r) const {return 0;}
};


struct non_dir_h_bond : public usable {
	fl good;
	fl bad;
	non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
		name = std::string("non_dir_h_bond(g=") + to_string(good) + ", b=" + to_string(bad) + ")";
	}
	fl eval(const atom& a, const atom& b, fl r) const {
		if(xs_h_bond_possible(a.xs, b.xs))
			return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
		return 0;
	}
	fl eval(sz t1, sz t2, fl r) const {
		if(xs_h_bond_possible(t1, t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		return 0;
	}

};

inline fl read_iterator(flv::const_iterator& i) {
	fl x = *i; 
	++i;
	return x;
}

fl smooth_div(fl x, fl y) {
	if(std::abs(x) < epsilon_fl) return 0;
	if(std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
	return x / y;
}

struct num_tors_add : public conf_independent {
	num_tors_add() { name = "num_tors_add"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		//fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		fl w = 0.1 * (read_iterator(i) + 1); // FIXME?
        std::cout << "TORS ADD, w=" << w << " num_tors=" << in.num_tors << "\n";
		return x + 1 + w * in.num_tors/5.0;
	}
};

struct num_tors_sqr : public conf_independent {
	num_tors_sqr() { name = "num_tors_sqr"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		return x + w * sqr(fl(in.num_tors)) / 5;
	}
};

struct num_tors_sqrt : public conf_independent {
	num_tors_sqrt() { name = "num_tors_sqrt"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		return x + w * std::sqrt(fl(in.num_tors)) / sqrt(5.0);
	}
};

struct num_tors_div : public conf_independent {
	num_tors_div() { name = "num_tors_div"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.1 * (read_iterator(i) + 1); // w is in [0..0.2]
		return smooth_div(x, 1 + w * in.num_tors/5.0);
	}
};

struct ligand_length : public conf_independent {
	ligand_length() { name = "ligand_length"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = read_iterator(i);
		return x + w * in.ligand_lengths_sum;
	}
};

struct num_ligands : public conf_independent {
	num_ligands() { name = "num_ligands"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 1 * read_iterator(i); // w is in [-1.. 1]
		return x + w * in.num_ligands;
	}
};

struct num_heavy_atoms_div : public conf_independent {
	num_heavy_atoms_div() { name = "num_heavy_atoms_div"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.05 * read_iterator(i); 
		return smooth_div(x, 1 + w * in.num_heavy_atoms); 
	}
};

struct num_heavy_atoms : public conf_independent {
	num_heavy_atoms() { name = "num_heavy_atoms"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.05 * read_iterator(i); 
		return x + w * in.num_heavy_atoms;
	}
};

struct num_hydrophobic_atoms : public conf_independent {
	num_hydrophobic_atoms() { name = "num_hydrophobic_atoms"; }
	sz size() const { return 1; }
	fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) const {
		fl w = 0.05 * read_iterator(i); 
		return x + w * in.num_hydrophobic_atoms;
	}
};

everything::everything(scoring_function_choice sfchoice) { // enabled according to design.out227

    switch(sfchoice) {
        case SF_VINA:
	        add(1, new gauss(0, 0.5, 8.0)); // offset, width, cutoff // WEIGHT: -0.035579
	        add(1, new gauss(3, 2.0, 8.0)); // offset, width, cutoff // WEIGHT: -0.005156
	        add(1, new repulsion( 0.0, 8.0)); // offset, cutoff // WEIGHT:  0.840245
	        add(1, new hydrophobic(0.5, 1.5, 8.0)); // good, bad, cutoff // WEIGHT:  -0.035069
	        add(1, new non_dir_h_bond(-0.7, 0, 8.0)); // good, bad, cutoff // WEIGHT:  -0.587439
	        add(1, new num_tors_div()); // WEIGHT: 1.923 -- FIXME too close to limit?
            break;
        case SF_VINARDO:
            std::cout << "\n\nVinardo scoring function is not implemented yet.\n\nAborting.\n\n";
            VINA_CHECK(false);
            break;
        case SF_AD42: 
	        add(1, new ad4_vdw( 0.5, 100000, 8.0)); // smoothing, cap, cutoff
	        add(1, new ad4_hb(  0.5, 100000, 8.0)); // smoothing, cap, cutoff
	        add(1, new electrostatic<1>(100, 20.48)); // exponent, cap, cutoff
	        add(1, new ad4_solvation(3.6, 0.01097,  true, 20.48)); // desolvation_sigma, solvation_q, charge_dependent, cutoff
	        add(1, new num_tors_add());

            break;
        default:
            std::cout << "INSIDE everything::everything()   sfchoice = " << sfchoice << "\n";
            VINA_CHECK(false);
            break;
    }

	//add(d, new ad4_solvation(3.6, 0.01097, false, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff

	//add(d, new electrostatic<2>(100, cutoff)); // cap, cutoff
    //
	//add(d, new ad4_hb<10,  12>( 0.5, 100, cutoff)); // smoothing, cap, cutoff

	//add(d, new gauss(0,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(0.5, 0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1.5, 0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2.5, 0.3, cutoff)); // offset, width, cutoff

	//add(1, new gauss(0, 0.5, cutoff)); // offset, width, cutoff // WEIGHT: -0.035579
	//add(d, new gauss(1, 0.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.5, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 0.7, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 0.7, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.7, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 0.9, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(4, 1.5, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 2.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 2.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 2.0, cutoff)); // offset, width, cutoff
	//add(1, new gauss(3, 2.0, cutoff)); // offset, width, cutoff // WEIGHT: -0.005156
	//add(d, new gauss(4, 2.0, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(4, 3.0, cutoff)); // offset, width, cutoff

	//add(d, new repulsion( 0.4, cutoff)); // offset, cutoff
	//add(d, new repulsion( 0.2, cutoff)); // offset, cutoff
	//add(1, new repulsion( 0.0, cutoff)); // offset, cutoff // WEIGHT:  0.840245
	//add(d, new repulsion(-0.2, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.4, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.6, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.8, cutoff)); // offset, cutoff
	//add(d, new repulsion(-1.0, cutoff)); // offset, cutoff

	//add(d, new hydrophobic(0.5, 1, cutoff)); // good, bad, cutoff
	//add(1, new hydrophobic(0.5, 1.5, cutoff)); // good, bad, cutoff // WEIGHT:  -0.035069
	//add(d, new hydrophobic(0.5, 2, cutoff)); // good, bad, cutoff
	//add(d, new hydrophobic(0.5, 3, cutoff)); // good, bad, cutoff

	//add(1, new non_hydrophobic(0.5, 1.5, cutoff));

	//add(d, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff

	//add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff // WEIGHT:  -0.587439
	//add(d, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff
	//add(d, new non_dir_h_bond(-0.7, 0.2, cutoff)); // good, bad, cutoff
	//add(d, new non_dir_h_bond(-0.7, 0.4, cutoff)); // good, bad, cutoff
	// additive

	// conf-independent
	//add(d, new num_ligands());

	//add(1, new num_tors_div()); // WEIGHT: 1.923 -- FIXME too close to limit?
	//add(d, new num_heavy_atoms_div());
	//add(d, new num_heavy_atoms());
	//add(1, new num_tors_add());
	//add(d, new num_tors_sqr());
	//add(d, new num_tors_sqrt());
	//add(d, new num_hydrophobic_atoms());
	///add(1, new ligand_length());

	//add(d, new num_tors(100, 100, false)); // cap, past_cap, heavy_only
	//add(1, new num_tors(100, 100,  true)); // cap, past_cap, heavy_only
	//add(d, new num_tors(  2,   1,  true)); // cap, past_cap, heavy_only
	//add(d, new num_heavy_atoms());
	//add(d, new ligand_max_num_h_bonds());
	//add(1, new num_ligands());
}
