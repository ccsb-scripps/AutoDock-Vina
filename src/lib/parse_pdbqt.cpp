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

#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <boost/utility.hpp> // for noncopyable 
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include "parse_pdbqt.h"
#include "atom_constants.h"
#include "file.h"
#include "convert_substring.h"
#include "parse_error.h"

struct stream_parse_error {
	unsigned line;
	std::string reason;
	stream_parse_error(unsigned line_, const std::string& reason_) : line(line_), reason(reason_) {}
	parse_error to_parse_error(const path& name) const {
		return parse_error(name, line, reason);
	}
};

struct parsed_atom : public atom {
	unsigned number; 
	parsed_atom(sz ad_, fl charge_, const vec& coords_, unsigned number_) : number(number_) {
		ad = ad_;
		charge = charge_;
		coords = coords_;
	}
};

void add_context(context& c, std::string& str) {
	c.push_back(parsed_line(str, boost::optional<sz>()));
}

std::string omit_whitespace(const std::string& str, sz i, sz j) {
	if(i < 1) i = 1;
	if(j < i-1) j = i-1; // i >= 1
	if(j < str.size()) j = str.size();

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	// omit trailing whitespace
	while(i <= j && std::isspace(str[j-1]))
		--j;

	VINA_CHECK(i-1 < str.size());
	VINA_CHECK(j-i+1 < str.size());

	return str.substr(i-1, j-i+1);
}

struct atom_syntax_error {
	std::string nature;
	atom_syntax_error(const std::string& nature_) : nature(nature_) {}
};

template<typename T>
T checked_convert_substring(const std::string& str, sz i, sz j, const std::string& dest_nature) {
	VINA_CHECK(i >= 1);
	VINA_CHECK(i <= j+1);
	if(j > str.size()) throw atom_syntax_error("The line is too short");

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	const std::string substr = str.substr(i-1, j-i+1);
	try {
		return boost::lexical_cast<T>(substr);
	}
	catch(...) {
		throw atom_syntax_error(std::string("\"") + substr + "\" is not a valid " + dest_nature);
	}
}

parsed_atom parse_pdbqt_atom_string(const std::string& str) {
	unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "atom number");
	vec coords(checked_convert_substring<fl>(str, 31, 38, "coordinate"),
			   checked_convert_substring<fl>(str, 39, 46, "coordinate"),
			   checked_convert_substring<fl>(str, 47, 54, "coordinate"));
	fl charge = 0;
	if(!substring_is_blank(str, 69, 76))
		charge = checked_convert_substring<fl>(str, 69, 76, "charge");
	std::string name = omit_whitespace(str, 78, 79);
	sz ad = string_to_ad_type(name);
	parsed_atom tmp(ad, charge, coords, number);
	if(is_non_ad_metal_name(name))
		tmp.xs = XS_TYPE_Met_D;
	if(tmp.acceptable_type()) 
		return tmp;
	else 
		throw atom_syntax_error(std::string("\"") + name + "\" is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.");
}

struct atom_reference {
	sz index;
	bool inflex;
	atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

struct movable_atom : public atom {
	vec relative_coords;
	movable_atom(const atom& a, const vec& relative_coords_) : atom(a) {
		relative_coords = relative_coords_;
	}
};

struct rigid {
	atomv atoms;
};

typedef std::vector<movable_atom> mav;

struct non_rigid_parsed {
	vector_mutable<ligand> ligands;
	vector_mutable<residue> flex;

	mav atoms;
	atomv inflex;

	distance_type_matrix atoms_atoms_bonds;
	matrix<distance_type> atoms_inflex_bonds;
	distance_type_matrix inflex_inflex_bonds;

	distance_type_matrix mobility_matrix() const {
		distance_type_matrix tmp(atoms_atoms_bonds);
		tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
		return tmp;
	}
};

struct parsing_struct {
	// start reading after this class
	template<typename T> // T == parsing_struct
	struct node_t {
		sz context_index;
		parsed_atom a;
		std::vector<T> ps;
		node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

		// inflex atom insertion
		void insert_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
			nr.inflex.push_back(a);
		}
		void insert_immobiles_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile_inflex(nr);
		}

		// insertion into non_rigid_parsed
		void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
			vec relative_coords; relative_coords = a.coords - frame_origin;
			c[context_index].second = nr.atoms.size();
			nr.atoms.push_back(movable_atom(a, relative_coords));
		}
		void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile(nr, c, frame_origin);
		}
	};

	typedef node_t<parsing_struct> node;
	boost::optional<sz> immobile_atom; // which of `atoms' is immobile, if any
	boost::optional<atom_reference> axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
	boost::optional<atom_reference> axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
	std::vector<node> atoms;

	void add(const parsed_atom& a, const context& c) { 
		VINA_CHECK(c.size() > 0);
		atoms.push_back(node(a, c.size()-1)); 
	}
	const vec& immobile_atom_coords() const {
		VINA_CHECK(immobile_atom);
		VINA_CHECK(immobile_atom.get() < atoms.size());
		return atoms[immobile_atom.get()].a.coords;
	}
	// inflex insertion
	void insert_immobile_inflex(non_rigid_parsed& nr) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.inflex.size(), true);
			atoms[immobile_atom.get()].insert_inflex(nr);
		}
	}

	// insertion into non_rigid_parsed
	void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.atoms.size(), false);
			atoms[immobile_atom.get()].insert(nr, c, frame_origin);
		}
	}

	bool essentially_empty() const { // no sub-branches besides immobile atom, including sub-sub-branches, etc
		VINA_FOR_IN(i, atoms) {
			if(immobile_atom && immobile_atom.get() != i)
				return false;
			const node& nd = atoms[i];
			if(!nd.ps.empty())
				return false; // FIXME : iffy
		}
		return true;
	}
};

unsigned parse_one_unsigned(const std::string& str, const std::string& start, unsigned count) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp;
	in_str >> tmp;
	if(!in_str || tmp < 0) 
		throw stream_parse_error(count, "Syntax error");
	return unsigned(tmp);
}

void parse_two_unsigneds(const std::string& str, const std::string& start, unsigned count, unsigned& first, unsigned& second) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp1, tmp2;
	in_str >> tmp1;
	in_str >> tmp2;
	if(!in_str || tmp1 < 0 || tmp2 < 0) 
		throw stream_parse_error(count, "Syntax error");
	first = unsigned(tmp1);
	second = unsigned(tmp2);
}

void parse_pdbqt_rigid(const path& name, rigid& r) {
	ifile in(name);
	unsigned count = 0;
	std::string str;
	while(std::getline(in, str)) {
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "TER")) {} // ignore 
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				r.atoms.push_back(parse_pdbqt_atom_string(str));
			}
			catch(atom_syntax_error& e) {
				throw parse_error(name, count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw parse_error(name, count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
}


void parse_pdbqt_root_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				p.add(parse_pdbqt_atom_string(str), c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "ENDROOT")) return;
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_root(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ROOT")) {
			parse_pdbqt_root_aux(in, count, p, c);
			break;
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to); // forward declaration

void parse_pdbqt_branch_aux(std::istream& in, unsigned& count, const std::string& str, parsing_struct& p, context& c) {
	unsigned first, second;
	parse_two_unsigneds(str, "BRANCH", count, first, second); 
	sz i = 0;
	for(; i < p.atoms.size(); ++i)
		if(p.atoms[i].a.number == first) {
			p.atoms[i].ps.push_back(parsing_struct());
			parse_pdbqt_branch(in, count, p.atoms[i].ps.back(), c, first, second);
			break;
		}
	if(i == p.atoms.size())
		throw stream_parse_error(count, "No atom number " + boost::lexical_cast<std::string>(first) + " in this branch");
}

void parse_pdbqt_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue) {
	parse_pdbqt_root(in, count, p, c);

	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(!residue && starts_with(str, "TORSDOF")) {
			if(torsdof) throw stream_parse_error(count, "TORSDOF can occur only once");
			torsdof = parse_one_unsigned(str, "TORSDOF", count);
		}
		else if(residue && starts_with(str, "END_RES")) return; 
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
	if(atm)
		VINA_RANGE(i, r.begin, r.end) {
			atom_reference& ar = atm.get();
			if(ar.inflex) 
				nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED; //(max_unsigned); // first index - atoms, second index - inflex
			else
				nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED; // (max_unsigned);
		}
}

void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
	if(axis_begin && axis_end) {
		atom_reference& r1 = axis_begin.get();
		atom_reference& r2 = axis_end  .get();
		if(r2.inflex) {
			VINA_CHECK(r1.inflex); // no atom-inflex rotors
			nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
		}
		else
			if(r1.inflex)
				nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (atoms, inflex)
			else
				nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
	}
}

typedef std::pair<sz, sz> axis_numbers;
typedef boost::optional<axis_numbers> axis_numbers_option;

void nr_update_matrixes(non_rigid_parsed& nr) {
	// atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)

	nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);  
	nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // first index - inflex, second index - atoms
	nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
}

template<typename B> // B == branch or main_branch or flexible_body 
void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
	b.node.begin = nr.atoms.size();
	VINA_FOR_IN(i, p.atoms) {  // postprocess atoms into 'b.node'
		parsing_struct::node& p_node = p.atoms[i];
		if(p.immobile_atom && i == p.immobile_atom.get()) {} // skip immobile_atom - it's already inserted in "THERE"
		else p_node.insert(nr, c, b.node.get_origin());
		p_node.insert_immobiles(nr, c, b.node.get_origin());
	}
	b.node.end = nr.atoms.size();

	nr_update_matrixes(nr);
	add_bonds(nr, p.axis_begin, b.node); // b.node is used as atom_range
	add_bonds(nr, p.axis_end  , b.node); // b.node is used as atom_range
	set_rotor(nr, p.axis_begin, p.axis_end);

	VINA_RANGE(i, b.node.begin, b.node.end)
		VINA_RANGE(j, i+1, b.node.end)
			nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME


	VINA_FOR_IN(i, p.atoms) { 	// postprocess children
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile already inserted // FIXME ?!
				b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch will assign begin and end
				postprocess_branch(nr, ps, c, b.children.back());
			}
		}
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
	VINA_CHECK(!p.atoms.empty());
	nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch will assign begin and end
	postprocess_branch(nr, p, c, nr.ligands.back());
	nr_update_matrixes(nr); // FIXME ?
}

void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		p_node.insert_inflex(nr);
		p_node.insert_immobiles_inflex(nr);
	}
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile atom already inserted // FIXME ?!
				nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords))); // postprocess_branch will assign begin and end
				postprocess_branch(nr, ps, c, nr.flex.back());
			}
		}
	}
	nr_update_matrixes(nr); // FIXME ?
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	unsigned count = 0;
	parsing_struct p;
	boost::optional<unsigned> torsdof;
	try {
		parse_pdbqt_aux(in, count, p, c, torsdof, false);
		if(p.atoms.empty()) 
			throw parse_error(name, count, "No atoms in the ligand");
		if(!torsdof)
			throw parse_error(name, count, "Missing TORSDOF");
		postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
	}
	catch(stream_parse_error& e) {
		throw e.to_parse_error(name);
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_residue(std::istream& in, unsigned& count, parsing_struct& p, context& c) { 
	boost::optional<unsigned> dummy;
	parse_pdbqt_aux(in, count, p, c, dummy, true);
}

void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	unsigned count = 0;
	std::string str;

	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BEGIN_RES")) {
			try {
				parsing_struct p;
				parse_pdbqt_residue(in, count, p, c);
				postprocess_residue(nr, p, c);
			}
			catch(stream_parse_error& e) {
				throw e.to_parse_error(name);
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} //ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(starts_with(str, "ENDBRANCH")) {
			unsigned first, second;
			parse_two_unsigneds(str, "ENDBRANCH", count, first, second);
			if(first != from || second != to) 
				throw stream_parse_error(count, "Inconsistent branch numbers");
			if(!p.immobile_atom) 
				throw stream_parse_error(count, "Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch");
			return;
		}
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				parsed_atom a = parse_pdbqt_atom_string(str);
				if(a.number == to)
					p.immobile_atom = p.atoms.size();
				p.add(a, c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}


//////////// new stuff //////////////////


struct pdbqt_initializer {
	model m;
	void initialize_from_rigid(const rigid& r) { // static really
		VINA_CHECK(m.grid_atoms.empty());
		m.grid_atoms = r.atoms;
	}
	void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
		VINA_CHECK(m.ligands.empty());
		VINA_CHECK(m.flex   .empty());

		m.ligands = nrp.ligands;
		m.flex    = nrp.flex;

		VINA_CHECK(m.atoms.empty());

		sz n = nrp.atoms.size() + nrp.inflex.size();
		m.atoms.reserve(n);
		m.coords.reserve(n);

		VINA_FOR_IN(i, nrp.atoms) {
			const movable_atom& a = nrp.atoms[i];
			atom b = static_cast<atom>(a);
			b.coords = a.relative_coords;
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_FOR_IN(i, nrp.inflex) {
			const atom& a = nrp.inflex[i];
			atom b = a;
			b.coords = zero_vec; // to avoid any confusion; presumably these will never be looked at
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_CHECK(m.coords.size() == n);

		m.internal_coords.resize(m.coords.size(), zero_vec); // FIXME

		m.minus_forces = m.coords;
		m.m_num_movable_atoms = nrp.atoms.size();

		if(is_ligand) {
			VINA_CHECK(m.ligands.size() == 1);
			m.ligands.front().cont = c;
		}
		else
			m.flex_context = c;

	}
	void initialize(const distance_type_matrix& mobility) {
		m.initialize(mobility);
	}
};

model parse_ligand_pdbqt  (const path& name) { // can throw parse_error
	non_rigid_parsed nrp;
	context c;
	parse_pdbqt_ligand(name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_nrp(nrp, c, true);
	tmp.initialize(nrp.mobility_matrix());
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name, const path& flex_name) { // can throw parse_error
	rigid r;
	non_rigid_parsed nrp;
	context c;
	parse_pdbqt_rigid(rigid_name, r);
	parse_pdbqt_flex(flex_name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	tmp.initialize_from_nrp(nrp, c, false);
	tmp.initialize(nrp.mobility_matrix());
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name) { // can throw parse_error
	rigid r;
	parse_pdbqt_rigid(rigid_name, r);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	distance_type_matrix mobility_matrix;
	tmp.initialize(mobility_matrix);
	return tmp.m;
}
