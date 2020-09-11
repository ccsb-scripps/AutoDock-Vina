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

#include "cache.h"
#include "model.h"
#include "precalculate.h"


cache::cache(const std::string& scoring_function_version_, const grid_dims& gd_, fl slope_, atom_type::t atom_typing_used_) 
: scoring_function_version(scoring_function_version_), gd(gd_), slope(slope_), atu(atom_typing_used_), grids(num_atom_types(atom_typing_used_)) {}

fl cache::eval(const model& m, fl v) const { // needs m.coords
	fl e = 0;
	sz nat = num_atom_types(atu);

	VINA_FOR(i, m.num_movable_atoms()) {
		const atom& a = m.atoms[i];
		sz t = a.get(atu);
		if(t >= nat) continue;
		const grid& g = grids[t];
		assert(g.initialized());
		e += g.evaluate(m.coords[i], slope, v);
	}
	return e;
}

fl cache::eval_intra(model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atu);

	VINA_FOR(i, m.num_movable_atoms()) {
        if(m.find_ligand(i) < m.ligands.size()) continue; // we only want flex-rigid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atu);
		if(t >= nat) continue;
		const grid& g = grids[t];
		assert(g.initialized());
		e += g.evaluate(m.coords[i], slope, v);
	}
	return e;
}

fl cache::eval_deriv(model& m, fl v) const { // needs m.coords, sets m.minus_forces
	fl e = 0;
	sz nat = num_atom_types(atu);

	VINA_FOR(i, m.num_movable_atoms()) {
		const atom& a = m.atoms[i];
		sz t = a.get(atu);
		if(t >= nat) { m.minus_forces[i].assign(0); continue; }
		const grid& g = grids[t];
		assert(g.initialized());
		vec deriv;
		e += g.evaluate(m.coords[i], slope, v, deriv);
		m.minus_forces[i] = deriv;
	}
	return e;
}

std::string convert_XS_to_string(sz t) {
	switch(t) {
		case XS_TYPE_C_H   : return "C_H";
		case XS_TYPE_C_P   : return "C_P";
		case XS_TYPE_N_P   : return "N_P";
		case XS_TYPE_N_D   : return "N_D";
		case XS_TYPE_N_A   : return "N_A";
		case XS_TYPE_N_DA  : return "N_DA";
		case XS_TYPE_O_P   : return "O_P";
		case XS_TYPE_O_D   : return "O_D";
		case XS_TYPE_O_A   : return "O_A";
		case XS_TYPE_O_DA  : return "O_DA";
		case XS_TYPE_S_P   : return "S_P";
		case XS_TYPE_P_P   : return "P_P";
		case XS_TYPE_F_H   : return "F_H";
		case XS_TYPE_Cl_H  : return "Cl_H";
		case XS_TYPE_Br_H  : return "Br_H";
		case XS_TYPE_I_H   : return "I_H";
		case XS_TYPE_Met_D : return "Met_D";
		default: VINA_CHECK(false);
	}
}

void cache::write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename,
                  const std::string& fld_filename, const std::string& receptor_filename) {
	std::string atom_type;
	std::string filename;

	VINA_FOR_IN(t, atom_types) {
		atom_type = convert_XS_to_string(t);
		filename = out_prefix + "." + atom_type + ".map";
		path p(filename);
		ofile out(p);

		// write header
        out << "GRID_PARAMETER_FILE " << gpf_filename << "\n";
        out << "GRID_DATA_FILE " << fld_filename << "\n";
        out << "MACROMOLECULE " << receptor_filename << "\n";

		// m_factor_inv is spacing
		// check that it's the same in every dimension (it must be)
		// check that == operator is OK
		if ((grids[t].m_factor_inv[0] != grids[t].m_factor_inv[1]) & (grids[t].m_factor_inv[0] != grids[t].m_factor_inv[2])) {
			printf("m_factor_inv x=%f, y=%f, z=%f\n", grids[t].m_factor_inv[0], grids[t].m_factor_inv[1], grids[t].m_factor_inv[2]);
			return;
		}

		out << "SPACING " << grids[t].m_factor_inv[0] << "\n";

		// The number of elements has to be an odd number. Who said AD was not odd?
        int size_x = (grids[t].m_data.dim0() % 2 == 0) ?  grids[t].m_data.dim0() - 1 : grids[t].m_data.dim0();
        int size_y = (grids[t].m_data.dim1() % 2 == 0) ?  grids[t].m_data.dim1() - 1 : grids[t].m_data.dim1();
        int size_z = (grids[t].m_data.dim2() % 2 == 0) ?  grids[t].m_data.dim2() - 1 : grids[t].m_data.dim2();
        out << "NELEMENTS " << size_x << " " << size_y  << " " << size_z << "\n";

		// center
		fl cx = grids[t].m_init[0] + grids[t].m_range[0] * 0.5;
		fl cy = grids[t].m_init[1] + grids[t].m_range[1] * 0.5;
		fl cz = grids[t].m_init[2] + grids[t].m_range[2] * 0.5;
		out << "CENTER " << cx << " " << cy << " " << cz << "\n";

		// write data
		VINA_FOR(z, grids[t].m_data.dim2()) {
			VINA_FOR(y, grids[t].m_data.dim1()) {
				VINA_FOR(x, grids[t].m_data.dim0()) {
					out << std::setprecision(4) << grids[t].m_data(x, y, z) << "\n"; // slow?
				} // x
			} // y
		} // z
	} // map atom type
} // cache::write

void cache::populate(const model& m, const precalculate& p, const szv& atom_types_needed) {
	szv needed;
	VINA_FOR_IN(i, atom_types_needed) {
		sz t = atom_types_needed[i];
		if(!grids[t].initialized()) {
			needed.push_back(t);
			grids[t].init(gd);
		}
	}
	if(needed.empty())
		return;
	flv affinities(needed.size());

	sz nat = num_atom_types(atu);

	grid& g = grids[needed.front()];

	const fl cutoff_sqr = p.cutoff_sqr();

	grid_dims gd_reduced = szv_grid_dims(gd);
	szv_grid ig(m, gd_reduced, cutoff_sqr);

	VINA_FOR(x, g.m_data.dim0()) {
		VINA_FOR(y, g.m_data.dim1()) {
			VINA_FOR(z, g.m_data.dim2()) {
				std::fill(affinities.begin(), affinities.end(), 0);
				vec probe_coords; probe_coords = g.index_to_argument(x, y, z);
				const szv& possibilities = ig.possibilities(probe_coords);
				VINA_FOR_IN(possibilities_i, possibilities) {
					const sz i = possibilities[possibilities_i];
					const atom& a = m.grid_atoms[i];
					const sz t1 = a.get(atu);
					if(t1 >= nat) continue;
					const fl r2 = vec_distance_sqr(a.coords, probe_coords);
					if(r2 <= cutoff_sqr) {
						VINA_FOR_IN(j, needed) {
							const sz t2 = needed[j];
							assert(t2 < nat);
							const sz type_pair_index = triangular_matrix_index_permissive(num_atom_types(atu), t1, t2);
							affinities[j] += p.eval_fast(type_pair_index, r2);
						}
					}
				}
				VINA_FOR_IN(j, needed) {
					sz t = needed[j];
					assert(t < nat);
					grids[t].m_data(x, y, z) = affinities[j];
				}
			}
		}
	}
}
