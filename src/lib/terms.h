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

#ifndef VINA_TERMS_H
#define VINA_TERMS_H

#include <boost/ptr_container/ptr_vector.hpp> 
#include "model.h"

struct term {
	std::string name;
	virtual ~term() {}
};

struct distance_additive : public term {
	fl cutoff;
	distance_additive(fl cutoff_) : cutoff(cutoff_) {}
	virtual fl eval(const atom_base& a, const atom_base& b, fl r) const = 0;
	virtual ~distance_additive() {}
};

struct usable : public distance_additive {
	atom_type::t atom_typing_used;
	usable(fl cutoff_) : distance_additive(cutoff_), atom_typing_used(atom_type::XS) {}
	fl eval(const atom_base& a, const atom_base& b, fl r) const { // should not be overriden
		return eval(a.get(atom_typing_used), b.get(atom_typing_used), r);
	}
	virtual fl eval(sz t1, sz t2, fl r) const { return 0; } 
	virtual ~usable() {}
};

struct additive : public term {
	fl cutoff;
	additive() : cutoff(max_fl) {}
	virtual fl eval(const model& m, const atom_index& i, const atom_index& j) const = 0;
	virtual ~additive() {}
};

struct intermolecular : public term {
	virtual fl eval(const model& m) const = 0;
};

struct conf_independent_inputs {
	fl num_tors;
	fl num_rotors;
	fl num_heavy_atoms;
	fl num_hydrophobic_atoms;
	fl ligand_max_num_h_bonds;
	fl num_ligands;
	fl ligand_lengths_sum;
	operator flv() const;
	conf_independent_inputs(const model& m);
	std::vector<std::string> get_names() const;
	conf_independent_inputs();
private:
	unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i) const; // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned atom_rotors(const model& m, const atom_index& i) const; // the number of rotatable bonds to heavy ligand atoms

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & num_tors;
		ar & num_rotors;
		ar & num_heavy_atoms;
		ar & num_hydrophobic_atoms;
		ar & ligand_max_num_h_bonds;
		ar & num_ligands;
		ar & ligand_lengths_sum;
	}
};

struct conf_independent : public term {
	virtual fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const = 0;
	virtual sz size() const = 0; // how many parameters does it take
}; 

template<typename T>
struct term_set {
	std::vector<bool> enabled;
	boost::ptr_vector<T> fun; // FIXME? const T?
	void add(unsigned e, T* f) { // FIXME? const T* ?
		enabled.push_back(e > 0);
		fun.push_back(f);
	}
	sz num_enabled() const {
		sz tmp = 0;
		VINA_FOR_IN(i, enabled)
			if(enabled[i])
				++tmp;
		return tmp;
	}
	void get_names(bool enabled_only, std::vector<std::string>& out) const { // appends to "out"
		VINA_CHECK(enabled.size() == fun.size());
		VINA_FOR_IN(i, fun)
			if(!enabled_only || enabled[i])
				out.push_back(fun[i].name);
	}
	void filter(flv::const_iterator& in, flv& out) const {
		VINA_CHECK(enabled.size() == fun.size());
		VINA_FOR_IN(i, enabled) {
			if(enabled[i])
				out.push_back(*in);
			++in;
		}
	}
	fl max_cutoff() const {
		fl tmp = 0;
		VINA_FOR_IN(i, fun)
			tmp = (std::max)(tmp, fun[i].cutoff);
		return tmp;
	}
	sz size() const { return fun.size(); }
	const T& operator[](sz i) const { return fun[i]; }
};


struct factors {
	flv e; // external
	flv i; // internal
	sz size() const { return e.size() + i.size(); }
	//sz num_weights() const { return (std::max)(e.size(), i.size()); } // FIXME? compiler bug? getting warnings here
	sz num_weights() const { return (e.size() > i.size()) ? e.size() : i.size(); }
	fl eval(const flv& weights, bool include_internal) const;
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & e;
		ar & i;
	}
};

struct terms {
	term_set<distance_additive> distance_additive_terms;
	term_set<usable>            usable_terms;
	term_set<additive>          additive_terms;
	term_set<intermolecular>    intermolecular_terms;
	term_set<conf_independent>  conf_independent_terms;

	// the class takes ownership of the pointer with 'add'
	void add(unsigned e, distance_additive* p) { distance_additive_terms.add(e, p); }
	void add(unsigned e, usable* p)            {            usable_terms.add(e, p); }
	void add(unsigned e, additive* p)          {          additive_terms.add(e, p); }
	void add(unsigned e, intermolecular* p)    {    intermolecular_terms.add(e, p); }
	void add(unsigned e, conf_independent* p)  {  conf_independent_terms.add(e, p); }

	std::vector<std::string> get_names(bool enabled_only) const; // does not include conf-independent
	sz size_internal() const;
	sz size() const { return size_internal() + intermolecular_terms.size(); }
	sz size_conf_independent(bool enabled_only) const; // number of parameters does not necessarily equal the number of operators
	fl max_r_cutoff() const;
	flv evale(const model& m) const;
	flv evali(const model& m) const;
	flv evale_robust(const model& m) const;
	factors eval(const model& m) const;
	fl eval_conf_independent(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const;
	flv filter_external(const flv& v) const;
	flv filter_internal(const flv& v) const;
	factors filter(const factors& f) const;
	void display_info() const;
private:
	void eval_additive_aux(const model& m, const atom_index& i, const atom_index& j, fl r, flv& out) const; // out is added to

};

#endif
