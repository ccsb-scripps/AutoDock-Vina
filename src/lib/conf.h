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

#ifndef VINA_CONF_H
#define VINA_CONF_H

#include <boost/ptr_container/ptr_vector.hpp> // typedef output_container

#include "quaternion.h"
#include "random.h"

struct scale {
	fl position;
	fl orientation;
	fl torsion;
	scale(fl position_, fl orientation_, fl torsion_) : position(position_), orientation(orientation_), torsion(torsion_) {}
};

struct conf_size {
	szv ligands;
	szv flex;
	sz num_degrees_of_freedom() const {
		return sum(ligands) + sum(flex) + 6 * ligands.size();
	}
};

inline void torsions_set_to_null(flv& torsions) {
	VINA_FOR_IN(i, torsions)
		torsions[i] = 0;
}

inline void torsions_increment(flv& torsions, const flv& c, fl factor) { // new torsions are normalized
	VINA_FOR_IN(i, torsions) {
		torsions[i] += normalized_angle(factor * c[i]);
		normalize_angle(torsions[i]);
	}
}

inline void torsions_randomize(flv& torsions, rng& generator) {
	VINA_FOR_IN(i, torsions)
		torsions[i] = random_fl(-pi, pi, generator);
}

inline bool torsions_too_close(const flv& torsions1, const flv& torsions2, fl cutoff) {
	assert(torsions1.size() == torsions2.size());
	VINA_FOR_IN(i, torsions1)
		if(std::abs(normalized_angle(torsions1[i] - torsions2[i])) > cutoff) 
			return false;
	return true;
}

inline void torsions_generate(flv& torsions, fl spread, fl rp, const flv* rs, rng& generator) {
	assert(!rs || rs->size() == torsions.size()); // if present, rs should be the same size as torsions
	VINA_FOR_IN(i, torsions)
		if(rs && random_fl(0, 1, generator) < rp)
			torsions[i] = (*rs)[i];
		else
			torsions[i] += random_fl(-spread, spread, generator);
}

struct rigid_change {
	vec position;
	vec orientation;
	rigid_change() : position(0, 0, 0), orientation(0, 0, 0) {}
	void print() const {
		::print(position);
		::print(orientation);
	}
};

struct rigid_conf {
	vec position;
	qt orientation;
	rigid_conf() : position(0, 0, 0), orientation(qt_identity) {}
	void set_to_null() {
		position = zero_vec;
		orientation = qt_identity;
	}
	void increment(const rigid_change& c, fl factor) {
		position += factor * c.position;
		vec rotation; rotation = factor * c.orientation;
		quaternion_increment(orientation, rotation); // orientation does not get normalized; tests show rounding errors growing very slowly
	}
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		position = random_in_box(corner1, corner2, generator);
		orientation = random_orientation(generator);
	}
	bool too_close(const rigid_conf& c, fl position_cutoff, fl orientation_cutoff) const {
		if(vec_distance_sqr(position, c.position) > sqr(position_cutoff)) return false;
		if(sqr(quaternion_difference(orientation, c.orientation)) > sqr(orientation_cutoff)) return false;
		return true;
	}
	void mutate_position(fl spread, rng& generator) {
		position += spread * random_inside_sphere(generator);
	}
	void mutate_orientation(fl spread, rng& generator) {
		vec tmp; tmp = spread * random_inside_sphere(generator);
		quaternion_increment(orientation, tmp);
	}
	void generate(fl position_spread, fl orientation_spread, fl rp, const rigid_conf* rs, rng& generator) {
		if(rs && random_fl(0, 1, generator) < rp)
			position = rs->position;
		else
			mutate_position(position_spread, generator);
		if(rs && random_fl(0, 1, generator) < rp)
			orientation = rs->orientation;
		else
			mutate_orientation(orientation_spread, generator);
	}
	void apply(const vecv& in, vecv& out, sz begin, sz end) const {
		assert(in.size() == out.size());
		const mat m = quaternion_to_r3(orientation);
		VINA_RANGE(i, begin, end)
			out[i] = m * in[i] + position;
	}
	void print() const {
		::print(position);
		::print(orientation);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & position;
		ar & orientation;
	}
};

struct ligand_change {
	rigid_change rigid;
	flv torsions;
	void print() const {
		rigid.print();
		printnl(torsions);
	}
};

struct ligand_conf {
	rigid_conf rigid;
	flv torsions;
	void set_to_null() {
		rigid.set_to_null();
		torsions_set_to_null(torsions);
	}
	void increment(const ligand_change& c, fl factor) {
		rigid.increment(c.rigid, factor);
		torsions_increment(torsions, c.torsions, factor);
	}
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		rigid.randomize(corner1, corner2, generator);
		torsions_randomize(torsions, generator);
	}
	void print() const {
		rigid.print();
		printnl(torsions);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & rigid;
		ar & torsions;
	}
};

struct residue_change {
	flv torsions;
	void print() const {
		printnl(torsions);
	}
};

struct residue_conf {
	flv torsions;
	void set_to_null() {
		torsions_set_to_null(torsions);
	}
	void increment(const residue_change& c, fl factor) {
		torsions_increment(torsions, c.torsions, factor);
	}
	void randomize(rng& generator) {
		torsions_randomize(torsions, generator);
	}
	void print() const {
		printnl(torsions);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & torsions;
	}
};

struct change {
	std::vector<ligand_change> ligands;
	std::vector<residue_change> flex;
	change(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0);
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0);
	}
	fl operator()(sz index) const { // returns by value
		VINA_FOR_IN(i, ligands) {
			const ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			const residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return 0; // shouldn't happen, placating the compiler
	}
	fl& operator()(sz index) {
		VINA_FOR_IN(i, ligands) {
			ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return ligands[0].rigid.position[0]; // shouldn't happen, placating the compiler
	}
	sz num_floats() const {
		sz tmp = 0;
		VINA_FOR_IN(i, ligands)
			tmp += 6 + ligands[i].torsions.size();
		VINA_FOR_IN(i, flex)
			tmp += flex[i].torsions.size();
		return tmp;
	}
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
};

struct conf {
	std::vector<ligand_conf> ligands;
	std::vector<residue_conf> flex;
	conf() {}
	conf(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0); // FIXME?
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0); // FIXME?
	}
	void set_to_null() {
		VINA_FOR_IN(i, ligands)
			ligands[i].set_to_null();
		VINA_FOR_IN(i, flex)
			flex[i].set_to_null();
	}
	void increment(const change& c, fl factor) { // torsions get normalized, orientations do not
		VINA_FOR_IN(i, ligands)
			ligands[i].increment(c.ligands[i], factor);
		VINA_FOR_IN(i, flex)
			flex[i]   .increment(c.flex[i],    factor);
	}
	bool internal_too_close(const conf& c, fl torsions_cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!torsions_too_close(ligands[i].torsions, c.ligands[i].torsions, torsions_cutoff))
				return false;
		return true;
	}
	bool external_too_close(const conf& c, const scale& cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!ligands[i].rigid.too_close(c.ligands[i].rigid, cutoff.position, cutoff.orientation))
				return false;
		assert(flex.size() == c.flex.size());
		VINA_FOR_IN(i, flex)
			if(!torsions_too_close(flex[i].torsions, c.flex[i].torsions, cutoff.torsion))
				return false;
		return true;
	}
	bool too_close(const conf& c, const scale& cutoff) const {
		return internal_too_close(c, cutoff.torsion) &&
			   external_too_close(c, cutoff); // a more efficient implementation is possible, probably
	}
	void generate_internal(fl torsion_spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid.position.assign(0);
			ligands[i].rigid.orientation = qt_identity;
			const flv* torsions_rs = rs ? (&rs->ligands[i].torsions) : NULL;
			torsions_generate(ligands[i].torsions, torsion_spread, rp, torsions_rs, generator);
		}
	}
	void generate_external(const scale& spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			const rigid_conf* rigid_conf_rs = rs ? (&rs->ligands[i].rigid) : NULL;
			ligands[i].rigid.generate(spread.position, spread.orientation, rp, rigid_conf_rs, generator);
		}
		VINA_FOR_IN(i, flex) {
			const flv* torsions_rs = rs ? (&rs->flex[i].torsions) : NULL;
			torsions_generate(flex[i].torsions, spread.torsion, rp, torsions_rs, generator);
		}
	}
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		VINA_FOR_IN(i, ligands)
			ligands[i].randomize(corner1, corner2, generator);
		VINA_FOR_IN(i, flex)
			flex[i].randomize(generator);
	}
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & ligands;
		ar & flex;
	}
};

struct output_type {
	conf c;
	fl e;
	fl lb;
	fl ub;
	fl intra;
	fl inter;
    fl conf_independent;
    fl unbound;
    fl total;
	vecv coords;
	output_type(const conf& c_, fl e_) : c(c_), e(e_) {}
	//output_type(const conf& c_, fl e_, fl intra_, fl conf_independent_) : c(c_), e(e_), intra(intra_), conf_independent(conf_independent_) {}
};

typedef boost::ptr_vector<output_type> output_container;

inline bool operator<(const output_type& a, const output_type& b) { // for sorting output_container
	return a.e < b.e;
}

#endif
