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

#include "pdb.h"
#include "parse_error.h"
#include "file.h"
#include "convert_substring.h"

void pdb::check(fl min_distance) const {
	VINA_FOR_IN(i, atoms) {
		const pdb_atom& a = atoms[i];
		VINA_RANGE(j, i+1, atoms.size()) {
			const pdb_atom& b = atoms[j];
			fl d2 = vec_distance_sqr(a.coords, b.coords);
			if(d2 < sqr(min_distance)) {
				std::cout << "The distance between " 
					<< a.id << ":" << a.name << ":" << a.element
					<< " and " 
					<< b.id << ":" << b.name << ":" << b.element
					<< " is " << std::sqrt(d2) << '\n';
			}
		}
	}
}

pdb_atom string_to_pdb_atom(const std::string& str) {
	if(str.size() < 66) throw bad_conversion(); // b-factor is in 61-66
	pdb_atom tmp;
	tmp.id           = convert_substring<unsigned>   (str,  7, 11);
	tmp.name         = convert_substring<std::string>(str, 13, 16);
	tmp.residue_id   = convert_substring<int>        (str, 23, 26);
	tmp.residue_name = convert_substring<std::string>(str, 18, 20);
	tmp.coords[0]    = convert_substring<fl>         (str, 31, 38);
	tmp.coords[1]    = convert_substring<fl>         (str, 39, 46);
	tmp.coords[2]    = convert_substring<fl>         (str, 47, 54);
	tmp.b_factor     = convert_substring<fl>         (str, 61, 66);
	tmp.element      = convert_substring<std::string>(str, 77, 78);
	return tmp;
}

pdb parse_pdb(const path& name) {
	ifile in(name);

	pdb tmp;

	std::string str;
	unsigned count = 0;
	while(std::getline(in, str)) {
		++count;
		if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				tmp.atoms.push_back(string_to_pdb_atom(str));
			}
			catch(...) { // bad_conversion, presumably; but can lexical_cast throw its own errors?
				throw parse_error(name, count, "ATOM syntax incorrect");
			}
		}
	}
	return tmp;
}
