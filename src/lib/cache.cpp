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


namespace fs = boost::filesystem;

std::string convert_XS_to_string(sz t) {
	switch(t) {
		case XS_TYPE_C_H     : return "C_H";
		case XS_TYPE_C_P     : return "C_P";
		case XS_TYPE_N_P     : return "N_P";
		case XS_TYPE_N_D     : return "N_D";
		case XS_TYPE_N_A     : return "N_A";
		case XS_TYPE_N_DA    : return "N_DA";
		case XS_TYPE_O_P     : return "O_P";
		case XS_TYPE_O_D     : return "O_D";
		case XS_TYPE_O_A     : return "O_A";
		case XS_TYPE_O_DA    : return "O_DA";
		case XS_TYPE_S_P     : return "S_P";
		case XS_TYPE_P_P     : return "P_P";
		case XS_TYPE_F_H     : return "F_H";
		case XS_TYPE_Cl_H    : return "Cl_H";
		case XS_TYPE_Br_H    : return "Br_H";
		case XS_TYPE_I_H     : return "I_H";
		case XS_TYPE_Si      : return "Si";
		case XS_TYPE_At      : return "At";
		case XS_TYPE_Met_D   : return "Met_D";
		case XS_TYPE_W       : return "W";
		default: VINA_CHECK(false);
	}
}

fl cache::eval(const model& m, fl v) const { // needs m.coords
	fl e = 0;
	sz nat = num_atom_types(atom_type::XS);

	VINA_FOR(i, m.num_movable_atoms()) {
		if(!m.is_atom_in_ligand(i)) continue; // we only want ligand - grid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::XS);

		if (t >= nat) { continue; }
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				t = XS_TYPE_C_H;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				t = XS_TYPE_C_P;
				break;
		}

		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v);
	}
	return e;
}

fl cache::eval_intra(model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atom_type::XS);

	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.is_atom_in_ligand(i)) continue; // we only want flex-rigid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::XS);

		if (t >= nat) { continue; }
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				t = XS_TYPE_C_H;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				t = XS_TYPE_C_P;
				break;
		}

		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v);
	}
	return e;
}

fl cache::eval_deriv(model& m, fl v) const { // needs m.coords, sets m.minus_forces
	fl e = 0;
	sz nat = num_atom_types(atom_type::XS);

	VINA_FOR(i, m.num_movable_atoms()) {
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::XS);

		if (t >= nat) { m.minus_forces[i].assign(0); continue; }
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				m.minus_forces[i].assign(0);
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				t = XS_TYPE_C_H;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				t = XS_TYPE_C_P;
				break;
		}

		vec deriv;
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v, deriv);
		m.minus_forces[i] = deriv;
	}
	return e;
}

bool cache::is_in_grid(const model& m, fl margin) const {
	VINA_FOR(i, m.num_movable_atoms()) {
		if (m.atoms[i].is_hydrogen()) continue;

		const vec& a_coords = m.coords[i];

		VINA_FOR_IN(j, m_gd) {
			if (m_gd[j].n_voxels > 0) {
				if (a_coords[j] < m_gd[j].begin - margin || a_coords[j] > m_gd[j].end + margin) {
					return false;
				}
			}
		}
	}

	return true;
}

bool cache::are_atom_types_grid_initialized(szv atom_types) const {
	sz nat = num_atom_types(atom_type::XS);

	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		if (t >= nat) { continue; }
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				t = XS_TYPE_C_H;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				t = XS_TYPE_C_P;
				break;
		}
		if (!is_atom_type_grid_initialized(t)) {
			std::cerr << "ERROR: Affinity map for atom type " <<  convert_XS_to_string(t) << " is not present.\n";
			return false;
		}
	}

	return true;
}

std::vector<std::string> vina_split(std::string str) {
	std::vector<std::string> fields;
	std::string field;
	std::istringstream iss(str);
	while (std::getline(iss, field, ' '))
	{
		fields.push_back(field);
	};
	return fields;
}

void read_vina_map(path &filename, std::vector<grid_dims> &gds, grid &g) {
	sz line_counter = 0;
	sz pt_counter = 0;
	sz x = 0;
	sz y = 0;
	sz z = 0;
	grid_dims gd;
	std::string line;
	fl spacing, center, halfspan;

	ifile in(filename);

	while (std::getline(in, line))
	{
		line_counter++;
		if (line_counter == 4)
		{
			std::vector<std::string> fields = vina_split(line);
			spacing = std::atof(fields[1].c_str());
		}
		if (line_counter == 5)
		{
			std::vector<std::string> fields = vina_split(line);
			VINA_FOR(i, 3) {
				// n_voxels must be EVEN
				// because the number of sampled points in the grid is always ODD
				// (number of sampled points == n_voxels + 1)
				gd[i].n_voxels = std::atoi(fields[i + 1].c_str());
				if (gd[i].n_voxels % 2 == 1) {
					std::cerr << "ERROR: number of voxels (NELEMENTS) must be even\n";
					exit(EXIT_FAILURE);
				}
			}
		}
		if (line_counter == 6)
		{
			std::vector<std::string> fields = vina_split(line);
			VINA_FOR(i, 3)
			{
				center = std::atof(fields[i + 1].c_str());
				halfspan = (gd[i].n_voxels) * spacing / 2.0;
				gd[i].begin = center - halfspan;
				gd[i].end = center + halfspan;
				// std::cout << center << " " << halfspan << " " << gd[i].begin << " " << gd[i].end << "\n";
			}
			gds.push_back(gd);
			g.init(gd);

		}
		if (line_counter > 6)
		{
			// std::cout << pt_counter << " " << x << " " << y << " " << z << " " << std::atof(line.c_str()) << "\n";
			g.m_data(x, y, z) = std::atof(line.c_str());
			y += sz(x == (gd[0].n_voxels + 1));
			z += sz(y == (gd[1].n_voxels + 1));
			x = x % (gd[0].n_voxels + 1);
			y = y % (gd[1].n_voxels + 1);
			pt_counter++;
			x++;
		}
	} // line loop
}

void cache::read(const std::string &map_prefix) {
	sz nat = num_atom_types(atom_type::XS);
	std::string type, filename;
	std::vector<grid_dims> gds; // to check all maps have same dims (TODO)

	bool got_C_H_already = false;
	bool got_C_P_already = false;
	bool found_at_least_1_map = false;

	VINA_FOR(atom_type, XS_TYPE_SIZE)
	{
		sz t = atom_type;

		if (t >= nat) { continue; }
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				if (got_C_H_already) continue;
				t = XS_TYPE_C_H;
				got_C_H_already = true;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				if (got_C_P_already) continue;
				t = XS_TYPE_C_P;
				got_C_P_already = true;
				break;
		}
		type = convert_XS_to_string(t);
		filename = map_prefix + "." + type + ".map";
		path p(filename);
		if (fs::exists(p)) {
			read_vina_map(p, gds, m_grids[t]);
			found_at_least_1_map = true;
		}
	} // map loop

	if (!found_at_least_1_map){
		std::cerr << "\nERROR: No *.map files with prefix \"" << map_prefix << "\"\n";
		exit(EXIT_FAILURE) ;
	}

	// Store in Cache object
	m_gd = gds[0];
}

void cache::write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename,
				      const std::string& fld_filename, const std::string& receptor_filename) {
	sz nat = num_atom_types(atom_type::XS);
	std::string atom_type;
	std::string filename;
	bool got_C_H_already = false;
	bool got_C_P_already = false;

	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		if (t >= nat) { continue; }
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				if (got_C_H_already) continue;
				t = XS_TYPE_C_H;
				got_C_H_already = true;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				if (got_C_P_already) continue;
				t = XS_TYPE_C_P;
				got_C_P_already = true;
				break;
		}
		if (m_grids[t].initialized()) {
			int nx = m_grids[t].m_data.dim0() - 1;
			int ny = m_grids[t].m_data.dim1() - 1;
			int nz = m_grids[t].m_data.dim2() - 1;
			atom_type = convert_XS_to_string(t);

			// For writing a .map, the number of points in the grid must be odd, which means that the number
			// of voxels (NELEMENTS) must be even... n_voxels = n_grid_points - 1
			if (nx % 2 == 1 || ny % 2 == 1 || nz % 2 == 1){
				std::cerr << "ERROR: Can't write maps. Number of voxels (NELEMENTS) is odd. Use --force_even_voxels.\n";
				exit(EXIT_FAILURE);
			} else {

				filename = out_prefix + "." + atom_type + ".map";
				path p(filename);
				ofile out(p);

				// write header
				out << "GRID_PARAMETER_FILE " << gpf_filename << "\n";
				out << "GRID_DATA_FILE " << fld_filename << "\n";
				out << "MACROMOLECULE " << receptor_filename << "\n";
				out << "SPACING " << m_grids[t].m_factor_inv[0] << "\n";
				
				out << "NELEMENTS " << nx << " " << ny << " " << nz << "\n";

				// center
				fl cx = m_grids[t].m_init[0] + m_grids[t].m_range[0] * 0.5;
				fl cy = m_grids[t].m_init[1] + m_grids[t].m_range[1] * 0.5;
				fl cz = m_grids[t].m_init[2] + m_grids[t].m_range[2] * 0.5;
				out << "CENTER " << cx << " " << cy << " " << cz << "\n";

				// write data
				VINA_FOR(z, m_grids[t].m_data.dim2())
				{
					VINA_FOR(y, m_grids[t].m_data.dim1())
					{
						VINA_FOR(x, m_grids[t].m_data.dim0())
						{
							out << std::setprecision(4) << m_grids[t].m_data(x, y, z) << "\n"; // slow?
						} // x
					} // y
				} // z
			} // even voxels
		} // map initialized
	} // map atom type
} // cache::write

void cache::populate(const model &m, const precalculate &p, const szv &atom_types_needed) {
	szv needed;
	bool got_C_H_already = false;
	bool got_C_P_already = false;

	VINA_FOR_IN(i, atom_types_needed) {
		sz t = atom_types_needed[i];
		switch (t)
		{
			case XS_TYPE_G0:
			case XS_TYPE_G1:
			case XS_TYPE_G2:
			case XS_TYPE_G3:
				continue;
			case XS_TYPE_C_H_CG0:
			case XS_TYPE_C_H_CG1:
			case XS_TYPE_C_H_CG2:
			case XS_TYPE_C_H_CG3:
				if (got_C_H_already) continue;
				t = XS_TYPE_C_H;
				got_C_H_already = true;
				break;
			case XS_TYPE_C_P_CG0:
			case XS_TYPE_C_P_CG1:
			case XS_TYPE_C_P_CG2:
			case XS_TYPE_C_P_CG3:
				if (got_C_P_already) continue;
				t = XS_TYPE_C_P;
				got_C_P_already = true;
				break;
		}
		if(!m_grids[t].initialized()) {
			needed.push_back(t);
			m_grids[t].init(m_gd);
		}
	}
	if(needed.empty())
		return;
	flv affinities(needed.size());

	sz nat = num_atom_types(atom_type::XS);

	grid& g = m_grids[needed.front()];

	const fl cutoff_sqr = p.cutoff_sqr();

	grid_dims gd_reduced = szv_grid_dims(m_gd);
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
					const sz t1 = a.get(atom_type::XS);
					if(t1 >= nat) continue;
					const fl r2 = vec_distance_sqr(a.coords, probe_coords);
					if(r2 <= cutoff_sqr) {
						VINA_FOR_IN(j, needed) {
							const sz t2 = needed[j];
							assert(t2 < nat);
							const sz type_pair_index = triangular_matrix_index_permissive(nat, t1, t2);
							affinities[j] += p.eval_fast(type_pair_index, r2);
						}
					}
				}
				VINA_FOR_IN(j, needed) {
					sz t = needed[j];
					assert(t < nat);
					m_grids[t].m_data(x, y, z) = affinities[j];
				}
			}
		}
	}
}
