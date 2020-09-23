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
        case AD_TYPE_W : return "W";
        default: VINA_CHECK(false);
    }
}

ad4cache::ad4cache(fl slope_): slope(slope_), grids(AD_TYPE_SIZE + 2) {}

fl ad4cache::eval(const model& m, fl v) const {
	fl e = 0;
    sz nat = num_atom_types(atom_type::AD);

    VINA_FOR(i, m.num_movable_atoms()) {
        const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);
        if (t == AD_TYPE_G0 || t == AD_TYPE_G1 || t == AD_TYPE_G2 || t == AD_TYPE_G3)
            continue;
        else if (t == AD_TYPE_CG0 || t == AD_TYPE_CG1 || t == AD_TYPE_CG2 || t == AD_TYPE_CG3)
            t = AD_TYPE_C;

        // HB + vdW
        const grid& g = grids[t];
		assert(g.initialized());
		e += g.evaluate(m.coords[i], slope, v);

        // elec
		const grid& ge = grids[AD_TYPE_SIZE];
		assert(ge.initialized());
		e += ge.evaluate(m.coords[i], slope, v) * a.charge;

        // desolv
		const grid& gd = grids[AD_TYPE_SIZE + 1];
		assert(gd.initialized());
		e += gd.evaluate(m.coords[i], slope, v) * std::abs(a.charge);

	}
	return e;
}

fl ad4cache::eval_intra(model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atom_type::AD);

    VINA_FOR(i, m.num_movable_atoms()) {
        if(m.find_ligand(i) < m.ligands.size()) continue; // we only want flex-rigid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);
        if (t == AD_TYPE_G0 || t == AD_TYPE_G1 || t == AD_TYPE_G2 || t == AD_TYPE_G3)
            continue;
        else if (t == AD_TYPE_CG0 || t == AD_TYPE_CG1 || t == AD_TYPE_CG2 || t == AD_TYPE_CG3)
            t = AD_TYPE_C;

        // HB + vdW
        const grid& g = grids[t];
		assert(g.initialized());
		e += g.evaluate(m.coords[i], slope, v);

        // elec
		const grid& ge = grids[AD_TYPE_SIZE];
		assert(ge.initialized());
		e += ge.evaluate(m.coords[i], slope, v) * a.charge;

        // desolv
		const grid& gd = grids[AD_TYPE_SIZE + 1];
		assert(gd.initialized());
		e += gd.evaluate(m.coords[i], slope, v) * std::abs(a.charge);
	}
	return e;
}

fl ad4cache::eval_deriv(model& m, fl v) const { // sets m.minus_forces
	fl e = 0;
    sz nat = num_atom_types(atom_type::AD);

    VINA_FOR(i, m.num_movable_atoms()) {
        const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);
        if (t==AD_TYPE_G0 || t==AD_TYPE_G1 || t==AD_TYPE_G2 || t==AD_TYPE_G3) {
            m.minus_forces[i].assign(0);
            continue;
        } else if (t==AD_TYPE_CG0|| t==AD_TYPE_CG1|| t==AD_TYPE_CG2|| t==AD_TYPE_CG3) {
            t = AD_TYPE_C;
        }

        // HB + vdW
        const grid& g = grids[t];
		assert(g.initialized());
		vec deriv;
		e += g.evaluate(m.coords[i], slope, v, deriv);
		m.minus_forces[i] = deriv;

        // elec
		const grid& ge = grids[AD_TYPE_SIZE];
		assert(ge.initialized());
		e += ge.evaluate(m.coords[i], slope, v, deriv) * a.charge;
        deriv *= a.charge;
		m.minus_forces[i] += deriv;

        // desolv
		const grid& gd = grids[AD_TYPE_SIZE + 1];
		assert(gd.initialized());
		e += gd.evaluate(m.coords[i], slope, v, deriv) * std::abs(a.charge);
        deriv *= std::abs(a.charge);
		m.minus_forces[i] += deriv;
	}
	return e;
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
    sz nx, ny, nz;

    ifile in(filename);

    while(std::getline(in, line)) {
        line_counter++;
        if (line_counter == 4) {
            std::vector<std::string> fields = split(line);
            spacing = std::atof(fields[1].c_str());
        }
        if (line_counter == 5) {
            std::vector<std::string> fields = split(line);
            nx = std::atoi(fields[1].c_str());
            ny = std::atoi(fields[2].c_str());
            nz = std::atoi(fields[3].c_str());
            gd[0].n = nx + sz(nx % 2 == 0) - 1; // n voxels, not sample points
            gd[1].n = ny + sz(nx % 2 == 0) - 1;
            gd[2].n = nz + sz(nx % 2 == 0) - 1;
            //std::cout << spacing << " " << gd[0].n << "\n";
        }
        if (line_counter == 6) {
            std::vector<std::string> fields = split(line);
            VINA_FOR(i, 3) {
                center = std::atof(fields[i+1].c_str());
                halfspan = (gd[i].n) * spacing / 2.0;
                gd[i].begin = center - halfspan;
                gd[i].end = center + halfspan;
            }
            gds.push_back(gd);
            g.init(gd);
        }
        if (line_counter > 6) {
            //  std::cout << line << "\n";
            g.m_data(x, y, z) = std::atof(line.c_str());
            y += sz(x == (gd[0].n + 1));
            z += sz(y == (gd[1].n + 1));
            x = x % (gd[0].n + 1);
            y = y % (gd[1].n + 1);
            pt_counter++;
            std::cout << x << " " << y << " " << z << " " << std::atof(line.c_str()) << "\n";
            x++;
        }

    } // line loop
}


grid_dims ad4cache::read(const std::string& map_prefix) {

    std::string type, filename;
    std::vector<grid_dims> gds; // to check all maps have same dims (TODO)

    bool got_C_already = false;
    sz t;

    VINA_FOR(ad_type, AD_TYPE_SIZE){
        t = ad_type;

        if ((t==AD_TYPE_G0) || (t==AD_TYPE_G1) ||
            (t==AD_TYPE_G2) || (t==AD_TYPE_G3))
            continue;
        if ((t==AD_TYPE_CG0) || (t==AD_TYPE_CG1) ||
            (t==AD_TYPE_CG2) || (t==AD_TYPE_CG3))
            t = AD_TYPE_C;
        if (t==AD_TYPE_C && got_C_already)
            continue;
        if (t==AD_TYPE_C)
            got_C_already=true;

        type = get_adtype_str(t);
        filename = map_prefix + "." + type + ".map";
        path p{filename};
        if (fs::exists(p)) {
            std::cout << "Reading " + filename + "\n";
            read_ad4_map(p, gds, grids[t]);

        } // if file exists
    } // map loop

    //  elec map
    filename = map_prefix + ".e.map";
    path pe{filename};
    std::cout << "Reading " + filename + "\n";
    read_ad4_map(pe, gds, grids[AD_TYPE_SIZE]);

    //  dsolv map
    filename = map_prefix + ".d.map";
    path pd{filename};
    std::cout << "Reading " + filename + "\n";
    read_ad4_map(pd, gds, grids[AD_TYPE_SIZE + 1]);

    // TODO verify grid_dims consistency

    return gds[0];
}

void ad4cache::write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename,
                     const std::string& fld_filename, const std::string& receptor_filename) {
    std::string atom_type;
    std::string filename;

    VINA_FOR_IN(t, atom_types) {
		atom_type = get_adtype_str(t);
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
