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

#include <boost/serialization/split_member.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/static_assert.hpp>
#include "ad4cache.h"
#include "file.h"
#include "szv_grid.h"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <boost/optional.hpp> // needed to write files?

namespace fs = boost::filesystem;

ad4cache::ad4cache(fl slope_): slope(slope_), rawgrids(AD_TYPE_SIZE + 2) {}

fl ad4cache::eval(const model& m, fl v) const {
	fl e = 0;

	VINA_FOR(i, m.num_movable_atoms()) {
        const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);
		const grid& g = rawgrids[t];
		assert(g.initialized());
		e += g.evaluate(m.coords[i], slope, v);

        // elec
		const grid& ge = rawgrids[AD_TYPE_SIZE];
		assert(ge.initialized());
		e += ge.evaluate(m.coords[i], slope, v) * a.charge;

        // desolv
		const grid& gd = rawgrids[AD_TYPE_SIZE + 1];
		assert(gd.initialized());
		e += gd.evaluate(m.coords[i], slope, v) * std::abs(a.charge);

	}
	return e;
}

fl ad4cache::eval_intra(model& m, fl v) const {
	fl e = 0;
	sz nat = num_atom_types(atu);

	VINA_FOR(i, m.num_movable_atoms()) {
        if(m.find_ligand(i) < m.ligands.size()) continue; // we only want flex-rigid interaction
		const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);
		const grid& g = rawgrids[t];
		assert(g.initialized());
		e += g.evaluate(m.coords[i], slope, v);

        // elec
		const grid& ge = rawgrids[AD_TYPE_SIZE];
		assert(ge.initialized());
		e += ge.evaluate(m.coords[i], slope, v) * a.charge;

        // desolv
		const grid& gd = rawgrids[AD_TYPE_SIZE + 1];
		assert(gd.initialized());
		e += gd.evaluate(m.coords[i], slope, v) * std::abs(a.charge);
	}
	return e;
}
fl ad4cache::eval_deriv(model& m, fl v) const { // sets m.minus_forces
	fl e = 0;

	VINA_FOR(i, m.num_movable_atoms()) {
        const atom& a = m.atoms[i];
		sz t = a.get(atom_type::AD);
		const grid& g = rawgrids[t];
		assert(g.initialized());
		vec deriv;
		e += g.evaluate(m.coords[i], slope, v, deriv);
		m.minus_forces[i] = deriv;

        // elec
		const grid& ge = rawgrids[AD_TYPE_SIZE];
		assert(ge.initialized());
		e += ge.evaluate(m.coords[i], slope, v, deriv) * a.charge;
        deriv *= a.charge;
		m.minus_forces[i] += deriv;

        // desolv
		const grid& gd = rawgrids[AD_TYPE_SIZE + 1];
		assert(gd.initialized());
		e += gd.evaluate(m.coords[i], slope, v, deriv) * std::abs(a.charge);
        deriv *= std::abs(a.charge);
		m.minus_forces[i] += deriv;
	}
	return e;
}

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
        default: VINA_CHECK(false);
    }
}

std::vector<std::string> split(std::string str) {
    std::istringstream iss(str);
    std::vector<std::string> fields{
        std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{}
    };
    return fields;
}

void readmap(path& filename, std::vector<grid_dims>& gds, grid& g) {

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
            // std::cout << x << " " << y << " " << z << std::atof(line.c_str()) << "\n";
            x++;
        }

    } // line loop
}


grid_dims ad4cache::read(const std::string& map_prefix) {

    std::string type, filename;
    std::vector<grid_dims> gds; // to check all maps have same dims (TODO)

    VINA_FOR(ad_type, AD_TYPE_SIZE){
        type = get_adtype_str(ad_type);
        filename = map_prefix + "." + type + ".map";
        path p{filename};
        if (fs::exists(p)) {
            std::cout << "Reading " + filename + "\n";
            readmap(p, gds, rawgrids[ad_type]);

        } // if file exists
    } // map loop

    //  elec map
    filename = map_prefix + ".e.map";
    path pe{filename};
    std::cout << "Reading " + filename + "\n";
    readmap(pe, gds, rawgrids[AD_TYPE_SIZE]);

    //  dsolv map
    filename = map_prefix + ".d.map";
    path pd{filename};
    std::cout << "Reading " + filename + "\n";
    readmap(pd, gds, rawgrids[AD_TYPE_SIZE + 1]);

    // TODO verify grid_dims consistency

    return gds[0];
}

void ad4cache::write(const std::string& out_prefix) {

    std::string type, filename;

    VINA_FOR(t, AD_TYPE_SIZE) {
		type = get_adtype_str(t);
        filename = out_prefix + "." + type + ".map";
        path p{filename};
        ofile out(p);

        // write header
        out << "GRID_PARAMETER_FILE NULL\n";
        out << "GRID_DATA_FILE NULL\n";
        out << "MACROMOLECULE NULL\n";

        // m_factor_inv is spacing
        // check that it's the same in every dimension (it must be)
        // check that == operator is OK
        if ((rawgrids[t].m_factor_inv[0] != rawgrids[t].m_factor_inv[1]) & (rawgrids[t].m_factor_inv[0] != rawgrids[t].m_factor_inv[2])) {
            printf("m_factor_inv x=%f, y=%f, z=%f\n", rawgrids[t].m_factor_inv[0], rawgrids[t].m_factor_inv[1], rawgrids[t].m_factor_inv[2]);
            return;
        }

        out << "SPACING " << rawgrids[t].m_factor_inv[0] << "\n";
        out << "NELEMENTS " << rawgrids[t].m_data.dim0() << " " << rawgrids[t].m_data.dim1() << " " << rawgrids[t].m_data.dim2() << "\n";

        // center
        fl cx = rawgrids[t].m_init[0] + rawgrids[t].m_range[0] * 0.5;
        fl cy = rawgrids[t].m_init[1] + rawgrids[t].m_range[1] * 0.5;
        fl cz = rawgrids[t].m_init[2] + rawgrids[t].m_range[2] * 0.5;
        printf("m_init %f\n",  rawgrids[t].m_init[0]);
        printf("m_range %f\n", rawgrids[t].m_range[0]);
        out << "CENTER " << cx << " " << cy << " " << cz << "\n";
    
        // write data
        VINA_FOR(z, rawgrids[t].m_data.dim2()) {
            VINA_FOR(y, rawgrids[t].m_data.dim1()) {
                VINA_FOR(x, rawgrids[t].m_data.dim0()) {
                    out << std::setprecision(4) << rawgrids[t].m_data(x, y, z) << "\n"; // slow?
                } // x
            } // y
        } // z
    } // map atom type
} // cache::write


//VINA_FOR_IN(i, fields){
//    std::cout << "Field: " << fields[i] << "\n";
//}
