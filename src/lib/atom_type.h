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

#ifndef VINA_ATOM_TYPE_H
#define VINA_ATOM_TYPE_H

#include "atom_constants.h"
#include "triangular_matrix_index.h"

struct atom_type {
	enum t {EL, AD, XS, SY};
	sz el, ad, xs, sy;
	atom_type() : el(EL_TYPE_SIZE), ad(AD_TYPE_SIZE), xs(XS_TYPE_SIZE), sy(SY_TYPE_SIZE) {}
	sz get(t atom_typing_used) const {
		switch(atom_typing_used) {
			case EL: return el;
			case AD: return ad;
			case XS: return xs;
			case SY: return sy;
			default: assert(false); return max_sz;
		}
	}
	bool is_hydrogen() const {
		return ad_is_hydrogen(ad);
	}
	bool is_heteroatom() const {
		return ad_is_heteroatom(ad) || xs == XS_TYPE_Met_D;
	}
	bool acceptable_type() const {
		return ad < AD_TYPE_SIZE || xs == XS_TYPE_Met_D;
	}
	void assign_el() {
		el = ad_type_to_el_type(ad);
		if(ad == AD_TYPE_SIZE && xs == XS_TYPE_Met_D)
			el = EL_TYPE_Met;
	}
	bool same_element(const atom_type& a) const { // does not distinguish metals or unassigned types
		return el == a.el;
	}
	fl covalent_radius() const {
		if(ad < AD_TYPE_SIZE)        return ad_type_property(ad).covalent_radius;
		else if(xs == XS_TYPE_Met_D) return metal_covalent_radius;
		VINA_CHECK(false);           return 0; // never happens - placating the compiler
	}
	fl optimal_covalent_bond_length(const atom_type& x) const {
		return covalent_radius() + x.covalent_radius();
	}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & el;
		ar & ad;
		ar & xs;
		ar & sy;
	}
};

inline sz num_atom_types(atom_type::t atom_typing_used) {
	switch(atom_typing_used) {
		case atom_type::EL: return EL_TYPE_SIZE;
		case atom_type::AD: return AD_TYPE_SIZE;
		case atom_type::XS: return XS_TYPE_SIZE;
		case atom_type::SY: return SY_TYPE_SIZE;
		default: assert(false); return max_sz;
	}
}

inline sz get_type_pair_index(atom_type::t atom_typing_used, const atom_type& a, const atom_type& b) { // throws error if any arg is unassigned in the given typing scheme
	sz n = num_atom_types(atom_typing_used);

	sz i = a.get(atom_typing_used); VINA_CHECK(i < n);
	sz j = b.get(atom_typing_used); VINA_CHECK(j < n);

	if(i <= j) return triangular_matrix_index(n, i, j);
	else       return triangular_matrix_index(n, j, i);
}

#endif
