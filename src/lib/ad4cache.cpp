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

#include "ad4cache.h"


namespace fs = boost::filesystem;

std::string get_adtype_str(sz& t) {
	switch(t) {
		case AD_TYPE_C : return "C";
		case AD_TYPE_A : return "A";
		case AD_TYPE_N : return "N";
		case AD_TYPE_O : return "O";
		case AD_TYPE_P : return "P";
		case AD_TYPE_S : return "S";
		case AD_TYPE_H : return "H";
		case AD_TYPE_F : return "F";
		case AD_TYPE_I : return "I";
		case AD_TYPE_NA: return "NA";
		case AD_TYPE_OA: return "OA";
		case AD_TYPE_SA: return "SA";
		case AD_TYPE_HD: return "HD";
		case AD_TYPE_Mg: return "Mg";
		case AD_TYPE_Mn: return "Mn";
		case AD_TYPE_Zn: return "Zn";
		case AD_TYPE_Ca: return "Ca";
		case AD_TYPE_Fe: return "Fe";
		case AD_TYPE_Cl: return "Cl";
		case AD_TYPE_Br: return "Br";
		case AD_TYPE_Si: return "Si";
		case AD_TYPE_At: return "At";
		case AD_TYPE_W : return "W";
		default: VINA_CHECK(false);
	}
}

fl ad4cache::eval(const model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

	VINA_FOR(i, m.num_movable_atoms()) {
		if(!m.is_atom_in_ligand(i)) continue; // we only want ligand interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		// HB + vdW
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v);

		// elec
		const grid& ge = m_grids[AD_TYPE_SIZE];
		e += ge.evaluate(m.coords[i], m_slope, v) * a.charge;

		// desolv
		const grid& gd = m_grids[AD_TYPE_SIZE + 1];
		e += gd.evaluate(m.coords[i], m_slope, v) * std::abs(a.charge);
	}
	return e;
}

fl ad4cache::eval_intra(model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.is_atom_in_ligand(i)) continue; // we only want flex-rigid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		// HB + vdW
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v);

		// elec
		const grid& ge = m_grids[AD_TYPE_SIZE];
		e += ge.evaluate(m.coords[i], m_slope, v) * a.charge;

		// desolv
		const grid& gd = m_grids[AD_TYPE_SIZE + 1];
		e += gd.evaluate(m.coords[i], m_slope, v) * std::abs(a.charge);
	}
	return e;
}

fl ad4cache::eval_deriv(model& m, fl v) const { // sets m.minus_forces
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

	VINA_FOR(i, m.num_movable_atoms()) {
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				m.minus_forces[i].assign(0);
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		// HB + vdW
		vec deriv;
		const grid& g = m_grids[t];
		e += g.evaluate(m.coords[i], m_slope, v, deriv);
		m.minus_forces[i] = deriv;

		// elec
		const grid& ge = m_grids[AD_TYPE_SIZE];
		e += ge.evaluate(m.coords[i], m_slope, v, deriv) * a.charge;
		deriv *= a.charge;
		m.minus_forces[i] += deriv;

		// desolv
		const grid& gd = m_grids[AD_TYPE_SIZE + 1];
		e += gd.evaluate(m.coords[i], m_slope, v, deriv) * std::abs(a.charge);
		deriv *= std::abs(a.charge);
		m.minus_forces[i] += deriv;
	}
	return e;
}

bool ad4cache::is_in_grid(const model& m, fl margin) const {
	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.atoms[i].is_hydrogen()) continue;

		const vec& a_coords = m.coords[i];
		VINA_FOR_IN(j, m_gd) {
			if(m_gd[j].n_voxels > 0)
				if(a_coords[j] < m_gd[j].begin - margin || a_coords[j] > m_gd[j].end + margin) 
					return false;
		}
	}
	return true;
}

bool ad4cache::are_atom_types_grid_initialized(szv atom_types) const {
	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		if (!is_atom_type_grid_initialized(t)) {
			std::cerr << "ERROR: Affinity map for atom type " << get_adtype_str(t) << " is not present.\n";
			return false;
		}
	}

	if (!is_atom_type_grid_initialized(AD_TYPE_SIZE)) {
		std::cerr << "ERROR: Electrostatic map is not present.\n";
		return false;
	}

	if (!is_atom_type_grid_initialized(AD_TYPE_SIZE + 1)) {
		std::cerr << "ERROR: Desolvation map is not present.\n";
		return false;
	}

	return true;
}

std::vector<std::string> split(std::string str) {
	std::vector<std::string> fields;
	std::string field;
	std::istringstream iss(str);
	while(std::getline(iss, field, ' '))
	{
		fields.push_back(field);
	};
	return fields;
}

void read_ad4_map(path& filename, std::vector<grid_dims>& gds, grid& g) {
	sz line_counter = 0;
	sz pt_counter = 0;
	sz x = 0;
	sz y = 0;
	sz z = 0;
	grid_dims gd;
	std::string line;
	fl spacing, center, halfspan;

	ifile in(filename);

	while(std::getline(in, line)) {
		line_counter++;
		if (line_counter == 4) {
			std::vector<std::string> fields = split(line);
			spacing = std::atof(fields[1].c_str());
		}
		if (line_counter == 5) {
			std::vector<std::string> fields = split(line);
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
		if (line_counter == 6) {
			std::vector<std::string> fields = split(line);
			VINA_FOR(i, 3) {
				center = std::atof(fields[i+1].c_str());
				halfspan = (gd[i].n_voxels) * spacing / 2.0;
				gd[i].begin = center - halfspan;
				gd[i].end = center + halfspan;
				// std::cout << center << " " << halfspan << " " << gd[i].begin << " " << gd[i].end << "\n";
			}
			gds.push_back(gd);
			g.init(gd);
		}
		if (line_counter > 6) {
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

void ad4cache::read(const std::string& map_prefix) {

	std::string type, filename;
	std::vector<grid_dims> gds; // to check all maps have same dims (TODO)

	bool got_C_already = false;

	VINA_FOR(atom_type, AD_TYPE_SIZE){
		sz t = atom_type;

		switch (t)
	{
		case AD_TYPE_G0:
		case AD_TYPE_G1:
		case AD_TYPE_G2:
		case AD_TYPE_G3:
			continue;
		case AD_TYPE_CG0:
		case AD_TYPE_CG1:
		case AD_TYPE_CG2:
		case AD_TYPE_CG3:
			if (got_C_already) continue;
			t = AD_TYPE_C;
			got_C_already = true;
			break;
	}

		type = get_adtype_str(t);
		filename = map_prefix + "." + type + ".map";
		path p(filename);
		if (fs::exists(p)) {
			read_ad4_map(p, gds, m_grids[t]);

		} // if file exists
	} // map loop

	//  elec map
	filename = map_prefix + ".e.map";
	path pe(filename);
	read_ad4_map(pe, gds, m_grids[AD_TYPE_SIZE]);

	//  dsolv map
	filename = map_prefix + ".d.map";
	path pd(filename);
	read_ad4_map(pd, gds, m_grids[AD_TYPE_SIZE + 1]);

	// Store in Ad4cache object
	m_gd = gds[0];
}

void ad4cache::write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename,
					      const std::string& fld_filename, const std::string& receptor_filename) {
	std::string atom_type;
	std::string filename;
	bool got_C_already = false;

	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		if (t == AD_TYPE_C && got_C_already)
			continue;
		if (t == AD_TYPE_C)
			got_C_already = true;

		if (m_grids[t].initialized()) {
			if (t < AD_TYPE_SIZE)
				atom_type = get_adtype_str(t);
			else if (t == AD_TYPE_SIZE)
				atom_type = "e";
			else if (t == AD_TYPE_SIZE + 1)
				atom_type = "d";

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
			if ((m_grids[t].m_factor_inv[0] != m_grids[t].m_factor_inv[1]) & (m_grids[t].m_factor_inv[0] != m_grids[t].m_factor_inv[2])) {
				printf("m_factor_inv x=%f, y=%f, z=%f\n", m_grids[t].m_factor_inv[0], m_grids[t].m_factor_inv[1], m_grids[t].m_factor_inv[2]);
				return;
			}

			out << "SPACING " << m_grids[t].m_factor_inv[0] << "\n";

			// The number of elements in the grid is an odd number. But NELEMENTS has to be an even number.
			int size_x = (m_grids[t].m_data.dim0() % 2 == 0) ? m_grids[t].m_data.dim0() : m_grids[t].m_data.dim0() - 1;
			int size_y = (m_grids[t].m_data.dim1() % 2 == 0) ? m_grids[t].m_data.dim1() : m_grids[t].m_data.dim1() - 1;
			int size_z = (m_grids[t].m_data.dim2() % 2 == 0) ? m_grids[t].m_data.dim2() : m_grids[t].m_data.dim2() - 1;
			out << "NELEMENTS " << size_x << " " << size_y  << " " << size_z << "\n";

			// center
			fl cx = m_grids[t].m_init[0] + m_grids[t].m_range[0] * 0.5;
			fl cy = m_grids[t].m_init[1] + m_grids[t].m_range[1] * 0.5;
			fl cz = m_grids[t].m_init[2] + m_grids[t].m_range[2] * 0.5;
			out << "CENTER " << cx << " " << cy << " " << cz << "\n";

			// write data
			VINA_FOR(z, m_grids[t].m_data.dim2()) {
				VINA_FOR(y, m_grids[t].m_data.dim1()) {
					VINA_FOR(x, m_grids[t].m_data.dim0()) {
						out << std::setprecision(4) << m_grids[t].m_data(x, y, z) << "\n"; // slow?
					} // x
				} // y
			} // z
		} // map initialized
	} // map atom type
} // cache::write
