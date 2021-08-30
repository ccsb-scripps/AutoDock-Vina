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

#include <random>

#include "model.h"
#include "file.h"
#include "curl.h"
#include "precalculate.h"

template<typename T>
atom_range get_atom_range(const T& t) {
	atom_range tmp = t.node;
	VINA_FOR_IN(i, t.children) {
		atom_range r = get_atom_range(t.children[i]);
		if(tmp.begin > r.begin) tmp.begin = r.begin;
		if(tmp.end   < r.end  ) tmp.end   = r.end;
	}
	return tmp;
}

struct branch_metrics {
	sz length;
	sz corner2corner;
	branch_metrics() : length(0), corner2corner(0) {}
};

template<typename T>
branch_metrics get_branch_metrics(const T& t) {
	branch_metrics tmp;
	if(!t.children.empty()) {
		sz corner2corner_max = 0;
		szv lengths;
		VINA_FOR_IN(i, t.children) {
			branch_metrics res = get_branch_metrics(t.children[i]);
			if(corner2corner_max < res.corner2corner)
				corner2corner_max = res.corner2corner;
			lengths.push_back(res.length + 1); // FIXME? weird compiler warning (sz -> unsigned)
		}
		std::sort(lengths.begin(), lengths.end());

		tmp.length = lengths.back();

		tmp.corner2corner = tmp.length;
		if(lengths.size() >= 2)
			tmp.corner2corner += lengths[lengths.size() - 1];

		if(tmp.corner2corner < corner2corner_max)
			tmp.corner2corner = corner2corner_max;
	}
	return tmp;
}

sz model::ligand_longest_branch(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).length;
}

sz model::ligand_length(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

void ligand::set_range() {
	atom_range tmp = get_atom_range(*this);
	begin = tmp.begin;
	end   = tmp.end;
}

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

struct appender_info {
	sz grid_atoms_size;
	sz m_num_movable_atoms;
	sz atoms_size;

	appender_info(const model& m) : grid_atoms_size(m.grid_atoms.size()), m_num_movable_atoms(m.m_num_movable_atoms), atoms_size(m.atoms.size()) {}
};

class appender {
	appender_info a_info;
	appender_info b_info;
	sz new_grid_index(sz x) const {
		return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
	}
public:
	bool is_a;

	appender(const model& a, const model& b) : a_info(a), b_info(b), is_a(true) {}

	sz operator()(sz x) const { // transform coord index
		if(is_a) {
			if(x < a_info.m_num_movable_atoms)  return x; // a-movable unchanged
			else                                return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
		}
		else {
			if(x < b_info.m_num_movable_atoms)  return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
			else                                return x + a_info.atoms_size; // all a's spliced before b-inflex
		}
	}
	atom_index operator()(const atom_index& x) const { // transform atom_index
		atom_index tmp(x);
		if(tmp.in_grid) tmp.i = new_grid_index(tmp.i);
		else            tmp.i = operator()(tmp.i);
		return tmp;
	}

	// type-directed old -> new transformations
	void update(interacting_pair& ip) const {
		ip.a = operator()(ip.a);
		ip.b = operator()(ip.b);
	}
	void update(vec& v) const { // coordinates & forces - do nothing
	}
	void update(ligand& lig) const {
		lig.transform(*this); // ligand as an atom_range subclass
		transform_ranges(lig, *this);
		VINA_FOR_IN(i, lig.pairs)
			this->update(lig.pairs[i]);
		VINA_FOR_IN(i, lig.cont)
			this->update(lig.cont[i]); // parsed_line update, below
	}
	void update(residue& r) const {
		transform_ranges(r, *this);
	}
	void update(parsed_line& p) const {
		if(p.second)
			p.second = operator()(p.second.get());
	}
	void update(atom& a) const {
		VINA_FOR_IN(i, a.bonds) {
			bond& b = a.bonds[i];
			b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
		}
	}

	// ligands, flex, flex_context, atoms; also used for other_pairs
	template<typename T>
	void append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
		sz a_sz = a.size();
		vector_append(a, b);

		is_a = true;
		VINA_FOR(i, a_sz)
			update(a[i]);

		is_a = false;
		VINA_RANGE(i, a_sz, a.size())
			update(a[i]);
	}

	//  coords, minus_forces, atoms
	template<typename T>
	void coords_append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
		std::vector<T> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

		is_a = true;
		VINA_FOR_IN(i, a)
			update(a[i]);

		is_a = false;
		VINA_FOR_IN(i, b_copy)
			update(b_copy[i]);

		// interleave 
		typedef typename std::vector<T>::const_iterator const cci;
		cci b1 = b_copy.begin();
		cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;
		cci b3 = b_copy.end();

		a.insert(a.begin() + a_info.m_num_movable_atoms , b1 , b2);
		a.insert(a.end()                                , b2 , b3);
	}
};

void model::append(const model& m) {
	VINA_CHECK(atom_typing_used() == m.atom_typing_used());

	appender t(*this, m);

	t.append(other_pairs, m.other_pairs);
	t.append(inter_pairs, m.inter_pairs);
	t.append(glue_pairs, m.glue_pairs);

	VINA_CHECK(minus_forces.size() == coords.size());
	VINA_CHECK(m.minus_forces.size() == m.coords.size());

	t.coords_append(coords, m.coords);
	t.coords_append(minus_forces, m.minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

	t.append(ligands, m.ligands);
	t.append(flex, m.flex);
	t.append(flex_context, m.flex_context);

	t.append(grid_atoms, m.grid_atoms);
	t.coords_append(atoms, m.atoms);

	// Add interaction pairs between previously added atoms and (very likely) ligand atoms
	/* Interactions:
	- flex     - ligand   : YES (inter_pairs)
	- flex_i   - flex_j   : YES (other_pairs) but append is used mostly for adding ligand
	- ligand_i - ligand_j : YES (inter_pairs)
	- macrocycle closure interactions: NO (1-2, 1-3, 1-4)
	*/
	VINA_FOR(i, m_num_movable_atoms) {
		VINA_RANGE(j, m_num_movable_atoms, m_num_movable_atoms + m.m_num_movable_atoms) {
			if (is_closure_clash(i, j)) continue; // 1-2, 1-3 or 1-4 interaction around CG-CG bond

			const atom& a = atoms[i];
			const atom& b = atoms[j];

			sz t1 = a.get(atom_typing_used());
			sz t2 = b.get(atom_typing_used());
			sz n = num_atom_types(atom_typing_used());

			if (t1 < n && t2 < n) {
				sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
				
				if (is_glue_pair(i, j)) {
					glue_pairs.push_back(interacting_pair(type_pair_index, i, j)); // glue_i - glue_j
				} else if (is_atom_in_ligand(i) && is_atom_in_ligand(j)) {
					inter_pairs.push_back(interacting_pair(type_pair_index, i, j)); // INTER: ligand_i - ligand_j
				} else if (is_atom_in_ligand(i) || is_atom_in_ligand(j)) {
					inter_pairs.push_back(interacting_pair(type_pair_index, i, j)); // INTER: flex - ligand
				} else {
					other_pairs.push_back(interacting_pair(type_pair_index, i, j)); // INTRA: flex_i - flex_j
				}
			}
		}
	}

	m_num_movable_atoms += m.m_num_movable_atoms;
}

///////////////////  end  MODEL::APPEND /////////////////////////


/////////////////// begin MODEL::INITIALIZE /////////////////////////

atom_index model::sz_to_atom_index(sz i) const {
	if(i < grid_atoms.size()) return atom_index(i                    ,  true);
	else                      return atom_index(i - grid_atoms.size(), false);
}

distance_type model::distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const {
	if(i.in_grid && j.in_grid) return DISTANCE_FIXED;
	if(i.in_grid) return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	if(j.in_grid) return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	assert(!i.in_grid);
	assert(!j.in_grid);
	assert(i.i < atoms.size());
	assert(j.i < atoms.size());
	sz a = i.i;
	sz b = j.i;
	if(a == b) return DISTANCE_FIXED;
	return (a < b) ? mobility(a, b) : mobility(b, a);
}

const vec& model::atom_coords(const atom_index& i) const {
	return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}

fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
	return vec_distance_sqr(atom_coords(a), atom_coords(b));
}

struct bond_less { // FIXME rm!?
	bool operator()(const bond& a, const bond& b) const {
		return a.connected_atom_index.i < b.connected_atom_index.i;
	}
};


bool model::atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const { // there is an atom closer to both a and b then they are to each other and immobile relative to them
	fl r2 = distance_sqr_between(a, b);
	VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
		sz i = relevant_atoms[relevant_atoms_i];
		atom_index c = sz_to_atom_index(i);
		if(a == c || b == c) continue;
		distance_type ac = distance_type_between(mobility, a, c);
		distance_type bc = distance_type_between(mobility, b, c);
		if(ac != DISTANCE_VARIABLE &&
		   bc != DISTANCE_VARIABLE &&
		   distance_sqr_between(a, c) < r2 &&
		   distance_sqr_between(b, c) < r2)
			return true;
	}
	return false;
}

struct beads {
	fl radius_sqr;
	std::vector<std::pair<vec, szv> > data;
	beads(sz reserve_size, fl radius_sqr_) : radius_sqr(radius_sqr_) { data.reserve(reserve_size); }
	void add(sz index, const vec& coords) {
		VINA_FOR_IN(i, data) {
			if(vec_distance_sqr(coords, data[i].first) < radius_sqr) {
				data[i].second.push_back(index);
				return;
			}
		}
		// not found
		std::pair<vec, szv> tmp;
		tmp.first = coords;
		tmp.second.push_back(index);
		data.push_back(tmp);
	}
};

void model::assign_bonds(const distance_type_matrix& mobility) { // assign bonds based on relative mobility, distance and covalent length
	const fl bond_length_allowance_factor = 1.1;
	sz n = grid_atoms.size() + atoms.size();

	// construct beads
	const fl bead_radius = 15;
	beads beads_instance(n, sqr(bead_radius));
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		beads_instance.add(i, atom_coords(i_atom_index));
	}
	// assign bonds
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		const vec& i_atom_coords = atom_coords(i_atom_index);
		atom& i_atom = get_atom(i_atom_index);

		const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
		fl i_atom_covalent_radius = max_covalent_r;
		if(i_atom.ad < AD_TYPE_SIZE)
			i_atom_covalent_radius = ad_type_property(i_atom.ad).covalent_radius;

		//find relevant atoms
		szv relevant_atoms;
		const fl bead_cutoff_sqr = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r));
		VINA_FOR_IN(b, beads_instance.data) {
			if(vec_distance_sqr(beads_instance.data[b].first, i_atom_coords) > bead_cutoff_sqr) continue;
			const szv& bead_elements = beads_instance.data[b].second;
			VINA_FOR_IN(bead_elements_i, bead_elements) {
				sz j = bead_elements[bead_elements_i];
				atom_index j_atom_index = sz_to_atom_index(j);
				atom& j_atom = get_atom(j_atom_index);
				const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
				distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
				if(dt != DISTANCE_VARIABLE && i != j) {
					fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
					//if(r2 < sqr(bond_length_allowance_factor * bond_length))
					if(r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)))
						relevant_atoms.push_back(j);
				}
			}
		}
		// find bonded atoms
		VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
			sz j = relevant_atoms[relevant_atoms_i];
			if(j <= i) continue; // already considered
			atom_index j_atom_index = sz_to_atom_index(j);
			atom& j_atom = get_atom(j_atom_index);
			const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
			distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
			fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
			if(r2 < sqr(bond_length_allowance_factor * bond_length) && !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms)) {
				bool rotatable = (dt == DISTANCE_ROTOR);
				fl length = std::sqrt(r2);
				i_atom.bonds.push_back(bond(j_atom_index, length, rotatable));
				j_atom.bonds.push_back(bond(i_atom_index, length, rotatable));
			}

		}
	}
}

bool model::bonded_to_HD(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).ad == AD_TYPE_HD) 
			return true;
	}
	return false;
}

bool model::bonded_to_heteroatom(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).is_heteroatom())
			return true;
	}
	return false;
}

void model::assign_types() {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		atom& a = get_atom(ai);
		a.assign_el();
		sz& x = a.xs;

		bool acceptor   = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA); // X-Score forumaltion apparently ignores SA
		bool donor_NorO = (a.el == EL_TYPE_Met || bonded_to_HD(a));

		switch(a.el) {
			case EL_TYPE_H    : break;
			case EL_TYPE_C    :{
                if     (a.ad == AD_TYPE_CG0){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG0 : XS_TYPE_C_H_CG0;}
                else if(a.ad == AD_TYPE_CG1){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG1 : XS_TYPE_C_H_CG1;}
                else if(a.ad == AD_TYPE_CG2){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG2 : XS_TYPE_C_H_CG2;}
                else if(a.ad == AD_TYPE_CG3){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG3 : XS_TYPE_C_H_CG3;}
                else                        {x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H;}
                break;
            }
			case EL_TYPE_N    : x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
			case EL_TYPE_O    : x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
			case EL_TYPE_S    : x = XS_TYPE_S_P; break;
			case EL_TYPE_P    : x = XS_TYPE_P_P; break;
			case EL_TYPE_F    : x = XS_TYPE_F_H; break;
			case EL_TYPE_Cl   : x = XS_TYPE_Cl_H; break;
			case EL_TYPE_Br   : x = XS_TYPE_Br_H; break;
			case EL_TYPE_I    : x = XS_TYPE_I_H; break;
			case EL_TYPE_Si   : x = XS_TYPE_Si; break;
			case EL_TYPE_At   : x = XS_TYPE_At; break;
			case EL_TYPE_Met  : x = XS_TYPE_Met_D; break;
            case EL_TYPE_Dummy: {
                if      (a.ad == AD_TYPE_G0) x = XS_TYPE_G0;
                else if (a.ad == AD_TYPE_G1) x = XS_TYPE_G1;
                else if (a.ad == AD_TYPE_G2) x = XS_TYPE_G2;
                else if (a.ad == AD_TYPE_G3) x = XS_TYPE_G3;
                else if (a.ad == AD_TYPE_W)  x = XS_TYPE_SIZE; // no W atoms in XS types
                else VINA_CHECK(false);
                break;
            }
			case EL_TYPE_SIZE : break;
			default: VINA_CHECK(false);
		}
	}
}

void model::bonded_to(sz a, sz n, szv& out) const {
	if(!has(out, a)) { // not found
		out.push_back(a);
		if(n > 0) 
			VINA_FOR_IN(i, atoms[a].bonds) {
				const bond& b = atoms[a].bonds[i];
				if(!b.connected_atom_index.in_grid)
					bonded_to(b.connected_atom_index.i, n-1, out);
			}
	}
}

szv model::bonded_to(sz a, sz n) const {
	szv tmp;
	bonded_to(a, n, tmp);
	return tmp;
}

// remove 1-2, 1-3 and 1-4 interactions around CG-CG bond
bool model::is_closure_clash(sz i, sz j) const {
	sz t1 = atoms[i].get(atom_type::AD);
	sz t2 = atoms[j].get(atom_type::AD);
    if ((t1==AD_TYPE_CG0 && t2==AD_TYPE_G0) || (t2==AD_TYPE_CG0 && t1==AD_TYPE_G0) ||
        (t1==AD_TYPE_CG1 && t2==AD_TYPE_G1) || (t2==AD_TYPE_CG1 && t1==AD_TYPE_G1) ||
        (t1==AD_TYPE_CG2 && t2==AD_TYPE_G2) || (t2==AD_TYPE_CG2 && t1==AD_TYPE_G2) ||
        (t1==AD_TYPE_CG3 && t2==AD_TYPE_G3) || (t2==AD_TYPE_CG3 && t1==AD_TYPE_G3))
            return false; // it's a G-CG pair, not a clash
	szv neighbors_of_i = bonded_to(i, 1); // up to 1-2 interactions
	szv neighbors_of_j = bonded_to(j, 1); // up to 1-2 interactions
	bool i_has_CG0 = false;
	bool i_has_CG1 = false;
	bool i_has_CG2 = false;
	bool i_has_CG3 = false;
	VINA_FOR_IN(index, neighbors_of_i) {
		sz type_i = atoms[neighbors_of_i[index]].get(atom_type::AD);
			if      (type_i==AD_TYPE_CG0) i_has_CG0=true;
			else if (type_i==AD_TYPE_CG1) i_has_CG1=true;
			else if (type_i==AD_TYPE_CG2) i_has_CG2=true;
			else if (type_i==AD_TYPE_CG3) i_has_CG3=true;
	}
	VINA_FOR_IN(index, neighbors_of_j) {
		sz type_j = atoms[neighbors_of_j[index]].get(atom_type::AD);
		if ((type_j==AD_TYPE_CG0 && i_has_CG0) ||
			(type_j==AD_TYPE_CG1 && i_has_CG1) ||
			(type_j==AD_TYPE_CG2 && i_has_CG2) ||
			(type_j==AD_TYPE_CG3 && i_has_CG3))
			return true;
	}
	return false;
}

bool model::is_glue_pair(sz i, sz j) const {
	sz t1 = atoms[i].get(atom_type::AD);
	sz t2 = atoms[j].get(atom_type::AD);

	if ((t1 == AD_TYPE_CG0 && t2 == AD_TYPE_G0) || (t2 == AD_TYPE_CG0 && t1 == AD_TYPE_G0) ||
		(t1 == AD_TYPE_CG1 && t2 == AD_TYPE_G1) || (t2 == AD_TYPE_CG1 && t1 == AD_TYPE_G1) ||
		(t1 == AD_TYPE_CG2 && t2 == AD_TYPE_G2) || (t2 == AD_TYPE_CG2 && t1 == AD_TYPE_G2) ||
		(t1 == AD_TYPE_CG3 && t2 == AD_TYPE_G3) || (t2 == AD_TYPE_CG3 && t1 == AD_TYPE_G3))
		return true;
	else
		return false;
}

void model::initialize_pairs(const distance_type_matrix& mobility) {
	/* Interactions:
	- ligand_i - ligand_i : YES (1-4 only) (ligand.pairs)
	- flex_i   - flex_i   : YES (1-4 only) (other_pairs)
	- flex_i   - flex_j   : YES (other_pairs)
	- intra macrocycle closure interactions: NO (1-2, 1-3, 1-4)
	*/
	VINA_FOR_IN(i, atoms) {
		sz i_lig = find_ligand(i);
		szv bonded_atoms = bonded_to(i, 3);   // up to 1-4

		VINA_RANGE(j, i + 1, atoms.size()) {
			if (mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
                if (is_closure_clash(i, j)) continue;  // 1-2, 1-3 or 1-4 interaction around CG-CG bond
				
				sz t1 = atoms[i].get(atom_typing_used());
				sz t2 = atoms[j].get(atom_typing_used());
				sz n  = num_atom_types(atom_typing_used());

				if (t1 < n && t2 < n) { // exclude, say, Hydrogens
					sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
					interacting_pair ip(type_pair_index, i, j);
					
					if (is_glue_pair(i, j)) {
						// Add glue_i - glue_i interaction pair
						glue_pairs.push_back(ip);
					} else if (i_lig < ligands.size() && find_ligand(j) == i_lig) {
						// Add INTRAmolecular ligand_i - ligand_i
						ligands[i_lig].pairs.push_back(ip);
					} else if (!is_atom_in_ligand(i) && !is_atom_in_ligand(j)) {
						// Add INTRAmolecular flex_i  - flex_i but also flex_i - flex_j
						other_pairs.push_back(ip);
					}
				}
			}
		}
	}
}

void model::initialize(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, ligands)
		ligands[i].set_range();
	assign_bonds(mobility);
	assign_types();
	initialize_pairs(mobility);
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////


sz model::num_internal_pairs() const {
	sz tmp = 0;
	VINA_FOR_IN(i, ligands)
		tmp += ligands[i].pairs.size();
	return tmp;
}

szv model::get_movable_atom_types(atom_type::t atom_typing_used_) const {
	szv tmp;
	sz n = num_atom_types(atom_typing_used_);
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a = atoms[i];
		sz t = a.get(atom_typing_used_);
		if(t < n && !has(tmp, t))
			tmp.push_back(t);
	}
	return tmp;
}

conf_size model::get_size() const {
	conf_size tmp;
	tmp.ligands = ligands.count_torsions();
	tmp.flex    = flex   .count_torsions();
	return tmp;
}

conf model::get_initial_conf() const { // torsions = 0, orientations = identity, ligand positions = current
	conf_size cs = get_size();
	conf tmp(cs);
	tmp.set_to_null();
	VINA_FOR_IN(i, ligands)
		tmp.ligands[i].rigid.position = ligands[i].node.get_origin();
	return tmp;
}

vecv model::get_ligand_coords() const { // FIXME rm
	VINA_CHECK(ligands.size() == 1);
	vecv tmp;
	const ligand &lig = ligands.front();
	VINA_RANGE(i, lig.begin, lig.end)
		tmp.push_back(coords[i]);
	return tmp;
}

std::vector<double> model::get_ligand_coords() {
	// Way to get coordinates out of the C++ world
	VINA_CHECK(ligands.size() == 1);
	std::vector<double> tmp;
	const ligand &lig = ligands.front();
	VINA_RANGE(i, lig.begin, lig.end) {
		tmp.push_back(coords[i][0]);
		tmp.push_back(coords[i][1]);
		tmp.push_back(coords[i][2]);
	}
	return tmp;
}

vecv model::get_heavy_atom_movable_coords() const { // FIXME mv
	vecv tmp;

	VINA_FOR(i, num_movable_atoms()) {
		if (atoms[i].el != EL_TYPE_H)
			tmp.push_back(coords[i]);
	}
	return tmp;
}

sz model::find_ligand(sz a) const {
	VINA_FOR_IN(i, ligands) {
		if(a >= ligands[i].begin && a < ligands[i].end)
			return i;
	}
	return ligands.size();
}

bool model::is_atom_in_ligand(sz a) const {
	VINA_FOR_IN(i, ligands) {
		if (a >= ligands[i].begin && a < ligands[i].end)
			return true;
	}
	return false;
}

bool model::is_movable_atom(sz a) const {
	if (a < num_movable_atoms())
		return true;
	else
		return false;
}

std::vector<double> model::center() const {
	std::vector<double> center(3, 0);

	VINA_FOR(i, num_movable_atoms()) {
		center[0] += coords[i][0];
		center[1] += coords[i][1];
		center[2] += coords[i][2];
	}

	center[0] /= num_movable_atoms();
	center[1] /= num_movable_atoms();
	center[2] /= num_movable_atoms();

	return center;
}

void string_write_coord(sz i, fl x, std::string& str) {
	VINA_CHECK(i > 0);
	--i;
	std::ostringstream out;
	out.setf(std::ios::fixed, std::ios::floatfield);
	out.setf(std::ios::showpoint);
	out << std::setw(8) << std::setprecision(3) << x;
	VINA_CHECK(out.str().size() == 8); 
	VINA_CHECK(str.size() > i + 8);
	VINA_FOR(j, 8)
		str[i+j] = out.str()[j];
}

std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) {
	std::string tmp(str);
	string_write_coord(31, coords[0], tmp);
	string_write_coord(39, coords[1], tmp);
	string_write_coord(47, coords[2], tmp);
	return tmp;
}

void model::write_context(const context& c, ofile& out) const {
	verify_bond_lengths();
	VINA_FOR_IN(i, c) {
		const std::string& str = c[i].first;
		if(c[i].second) {
			out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
		}
		else
			out << str << '\n';
	}
}

void model::write_context(const context &c, std::ostringstream& out) const {
	verify_bond_lengths();

	VINA_FOR_IN(i, c) {
		const std::string &str = c[i].first;
		if (c[i].second)
			out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
		else
			out << str << '\n';
	}
}

std::string model::write_model(sz model_number, const std::string &remark) {
	std::ostringstream out;

	out << "MODEL " << model_number << '\n';
	out << remark;

	VINA_FOR_IN(i, ligands)
		write_context(ligands[i].cont, out);
	if (num_flex() > 0) // otherwise remark is written in vain
		write_context(flex_context, out);

	out << "ENDMDL\n";

	return out.str();
}

void model::set         (const conf& c) {
	ligands.set_conf(atoms, coords, c.ligands);
	flex   .set_conf(atoms, coords, c.flex);
}

fl model::gyration_radius(sz ligand_number) const {
	VINA_CHECK(ligand_number < ligands.size());
	const ligand& lig = ligands[ligand_number];
	fl acc = 0;
	unsigned counter = 0;
	VINA_RANGE(i, lig.begin, lig.end) {
		if(atoms[i].el != EL_TYPE_H) { // only heavy atoms are used
			acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
			++counter;
		}
	}
	return (counter > 0) ? std::sqrt(acc/counter) : 0;
}


fl eval_interacting_pairs(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, const bool with_max_cutoff) { // clean up
	fl e = 0;
	fl cutoff_sqr = p.cutoff_sqr();

	if (with_max_cutoff) {
		cutoff_sqr = p.max_cutoff_sqr();
	}

	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);
		if(r2 < cutoff_sqr) {
			fl tmp = p.eval_fast(ip.a, ip.b, r2);
			curl(tmp, v);
			e += tmp;
		}
	}
	return e;
}

fl eval_interacting_pairs_deriv(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces, const bool with_max_cutoff) { // adds to forces  // clean up
	fl e = 0;
	fl cutoff_sqr = p.cutoff_sqr();

	if (with_max_cutoff) {
		cutoff_sqr = p.max_cutoff_sqr();
	}

	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		vec r = coords[ip.b] - coords[ip.a]; // a -> b
		fl r2 = sqr(r);
		if(r2 < cutoff_sqr) {
			pr tmp = p.eval_deriv(ip.a, ip.b, r2);
			vec force; force = tmp.second * r;
			curl(tmp.first, force, v);
			e += tmp.first;

			// FIXME inefficient, if using hard curl
			forces[ip.a] -= force; // we could omit forces on inflex here
			forces[ip.b] += force;
		}
	}
	return e;
}

fl model::evalo(const precalculate_byatom& p, const vec& v) const { // clean up
	fl e = eval_interacting_pairs(p, v[2], other_pairs, coords);
	return e;
}

fl model::eval_inter(const precalculate_byatom& p, const vec& v) const { // clean up
	fl e = eval_interacting_pairs(p, v[2], inter_pairs, coords);
	return e;
}

fl model::evali(const precalculate_byatom& p, const vec& v) const { // clean up
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // probably might was well use coords here
	return e;
}

fl model::eval_deriv(const precalculate_byatom& p, const igrid& ig, const vec& v, change& g) { // clean up
	// INTER ligand - grid
	fl e = ig.eval_deriv(*this, v[1]); // sets minus_forces, except inflex

	// INTRA ligand_i - ligand_i
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords, minus_forces); // adds to minus_forces

	// INTER ligand_i - ligand_j and ligand_i - flex_i
	if (!inter_pairs.empty()) 
		e += eval_interacting_pairs_deriv(p, v[2], inter_pairs, coords, minus_forces); // adds to minus_forces
	// INTRA flex_i - flex_i and flex_i - flex_j
	if (!other_pairs.empty())
		e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
	// glue_i - glue_i and glue_i - glue_j
	if (!glue_pairs.empty())
		e += eval_interacting_pairs_deriv(p, v[2], glue_pairs, coords, minus_forces, true); // adds to minus_forces

	// calculate derivatives
	ligands.derivative(coords, minus_forces, g.ligands);
	flex.derivative(coords, minus_forces, g.flex); // inflex forces are ignored
	return e;
}

fl model::eval_intramolecular(const precalculate_byatom& p, const igrid& ig, const vec& v) {
	fl e = 0;
	const fl cutoff_sqr = p.cutoff_sqr();

	// internal for each ligand
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

	// flex - rigid
    e += ig.eval_intra(*this, v[1]);

	// flex_i - flex_i and flex_i - flex_j
	VINA_FOR_IN(i, other_pairs) {
		const interacting_pair& pair = other_pairs[i];
		fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
		if (r2 < cutoff_sqr) {
			fl this_e = p.eval_fast(pair.a, pair.b, r2);
			curl(this_e, v[2]);
			e += this_e;
		}
	}

	// glue_i - glue_i and glue_i - glue_j interactions (no cutoff)
	VINA_FOR_IN(i, glue_pairs) {
		const interacting_pair& pair = glue_pairs[i];
		fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
		fl this_e = p.eval_fast(pair.a, pair.b, r2);
		curl(this_e, v[2]);
		e += this_e;
	}

	return e;
}

fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
	sz n = x.m_num_movable_atoms; 
	VINA_CHECK(n == y.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, n) {
		const atom& a =   x.atoms[i];
		if(a.el != EL_TYPE_H) {
			fl r2 = max_fl;
			VINA_FOR(j, n) {
				const atom& b = y.atoms[j];
				if(a.same_element(b) && !b.is_hydrogen()) {
					fl this_r2 = vec_distance_sqr(x.coords[i], 
					                              y.coords[j]);
					if(this_r2 < r2)
						r2 = this_r2;
				}
			}
			assert(not_max(r2));
			sum += r2;
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_lower_bound(const model& m) const {
	return (std::max)(rmsd_lower_bound_asymmetric(*this, m), rmsd_lower_bound_asymmetric(m, *this));
}

fl model::rmsd_upper_bound(const model& m) const {
	VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a =   atoms[i];
		const atom& b = m.atoms[i];
		assert(a.ad == b.ad);
		assert(a.xs == b.xs);
		if(a.el != EL_TYPE_H) {
			sum += vec_distance_sqr(coords[i], m.coords[i]);
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_ligands_upper_bound(const model& m) const {
	VINA_CHECK(ligands.size() == m.ligands.size());
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR_IN(ligand_i, ligands) {
		const ligand&   lig =   ligands[ligand_i];
		const ligand& m_lig = m.ligands[ligand_i];
		VINA_CHECK(lig.begin == m_lig.begin);
		VINA_CHECK(lig.end   == m_lig.end);
		VINA_RANGE(i, lig.begin, lig.end) {
			const atom& a =   atoms[i];
			const atom& b = m.atoms[i];
			assert(a.ad == b.ad);
			assert(a.xs == b.xs);
			if(a.el != EL_TYPE_H) {
				sum += vec_distance_sqr(coords[i], m.coords[i]);
				++counter;
			}
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}


void model::verify_bond_lengths() const {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		const atom& a = get_atom(ai);
		VINA_FOR_IN(j, a.bonds) {
			const bond& b = a.bonds[j];
			fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
			bool ok = eq(d, b.length);
			if(!ok) {
				VINA_SHOW(d);
				VINA_SHOW(b.length);
			}
			VINA_CHECK(ok);
		}
	}
}

void model::check_ligand_internal_pairs() const {
	VINA_FOR_IN(i, ligands) {
		const ligand& lig = ligands[i];
		const interacting_pairs& pairs = lig.pairs;
		VINA_FOR_IN(j, pairs) {
			const interacting_pair& ip = pairs[j];
			VINA_CHECK(ip.a >= lig.begin);
			VINA_CHECK(ip.b  < lig.end);
		}
	}
}

void model::about() const {
	VINA_SHOW(atom_typing_used());
	VINA_SHOW(num_movable_atoms());
	VINA_SHOW(num_internal_pairs());
	VINA_SHOW(num_other_pairs());
	VINA_SHOW(num_ligands());
	VINA_SHOW(num_flex());
}

void model::show_pairs() const {

    std::cout << "INTER PAIRS\n";
    interacting_pairs inter_pairs = get_inter_pairs();
    VINA_FOR_IN(i, inter_pairs) {
        const interacting_pair &ip = inter_pairs[i];

		if (is_atom_in_ligand(ip.a)) {
			sz lig_i = find_ligand(ip.a);
			std::cout << "LIGAND (" << lig_i << ") : ";
		} else {
			std::cout << "  FLEX     : ";
		}

		if (is_atom_in_ligand(ip.b)) {
			sz lig_i = find_ligand(ip.b);
			std::cout << "LIGAND (" << lig_i << ") ";
		} else {
			std::cout << "  FLEX     ";
		}

		std::cout << " - " << ip.a << " : " << ip.b
				  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
				  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
				  << "\n";
	}

	std::cout << "INTRA LIG PAIRS\n";
    VINA_FOR(i, num_ligands()) {
        ligand lig = get_ligand(i);
        VINA_FOR_IN(j, lig.pairs) {
            const interacting_pair &ip = lig.pairs[j];
			sz lig_i = find_ligand(ip.a);
			std::cout << "LIGAND (" << lig_i << ") ";
			std::cout << " - " << ip.a << " : " << ip.b
					  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
					  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
					  << "\n";
		}
    }

    std::cout << "INTRA FLEX PAIRS\n";
    interacting_pairs other_pairs = get_other_pairs();
    VINA_FOR_IN(i, other_pairs) {
		const interacting_pair& ip = other_pairs[i];
		std::cout << "FLEX       ";
		std::cout << " - " << ip.a << " : " << ip.b
				  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
				  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
				  << "\n";
	}

	std::cout << "GLUE - GLUE PAIRS\n";
    interacting_pairs glue_pairs = get_glue_pairs();
    VINA_FOR_IN(i, glue_pairs) {
		const interacting_pair& ip = glue_pairs[i];
		std::cout << "FLEX       ";
		std::cout << " - " << ip.a << " : " << ip.b
				  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
				  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
				  << "\n";
	}
}

void model::show_atoms() const {
	std::cout << "ATOM INFORMATION\n";
	VINA_FOR_IN(i, atoms) {
		const atom &a = atoms[i];
		if (i < num_movable_atoms()) {
			std::cout << "     MOVABLE: ";
		} else {
			std::cout << " NOT MOVABLE: ";
		}
		std::cout << i << " - " << coords[i][0] << " " << coords[i][1] << " " << coords[i][2]
				  << " - " << a.ad << " - " << a.xs << " - " << a.charge << "\n";
	}
}

void model::show_forces() const {
	std::cout << "ATOM FORCES\n";
	VINA_FOR_IN(i, atoms) {
		std::cout << i << " " << minus_forces[i][0] << " " << minus_forces[i][1] << " " << minus_forces[i][2] << "\n";
	}
}

void model::print_stuff(bool show_coords, bool show_internal, bool show_atoms, bool show_grid, bool show_about) const {
	
	if (show_coords) {
		std::cout << "coords:\n";
		VINA_FOR_IN(i, coords)
			printnl(coords[i]);
	}

	if (show_atoms) {
		std::cout << "atoms:\n";
		VINA_FOR_IN(i, atoms) {
			const atom& a = atoms[i];
			std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
			std::cout << a.bonds.size() << "  "; printnl(a.coords);
		}
	}

	if (show_grid) {
		std::cout << "grid_atoms:\n";
		VINA_FOR_IN(i, grid_atoms) {
			const atom& a = grid_atoms[i];
			std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
			std::cout << a.bonds.size() << "  "; printnl(a.coords);
		}
	}

	if (show_about) {
		about();
	}
}

fl pairwise_clash_penalty(fl r, fl covalent_r) {
	// r = 0          -> max_penalty 
	// r = covalent_r -> 1
	// elsewhere      -> hyperbolic function
	assert(r >= 0);
	assert(covalent_r > epsilon_fl);
	const fl x = r / covalent_r;
	if(x > 2) return 0;
	return 1-x*x/4;
}

fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
		const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();
		e += pairwise_clash_penalty(r, covalent_r);
	}
	return e;
}

fl model::clash_penalty() const {
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += clash_penalty_aux(ligands[i].pairs);
	e += clash_penalty_aux(other_pairs);
	return e;
}
