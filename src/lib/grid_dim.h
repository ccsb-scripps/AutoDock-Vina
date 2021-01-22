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

#ifndef VINA_GRID_DIM_H
#define VINA_GRID_DIM_H

#include <boost/array.hpp>

#include "common.h"

struct grid_dim {
	fl begin;
	fl end;
	sz n_voxels; // number of intervals == number of sample points - 1
	grid_dim() : begin(0), end(0), n_voxels(0) {}
	fl span() const { return end - begin; }
	bool enabled() const { return (n_voxels > 0); }
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & begin;
		ar & end;
		ar & n_voxels;
	}
};

inline bool eq(const grid_dim& a, const grid_dim& b) {
	return a.n_voxels == b.n_voxels && eq(a.begin, b.begin) && eq(a.end, b.end);
}

typedef boost::array<grid_dim, 3> grid_dims;

inline bool eq(const grid_dims& a, const grid_dims& b) {
	return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

inline void print(const grid_dims& gd, std::ostream& out = std::cout) {
	VINA_FOR_IN(i, gd)
		std::cout << gd[i].n_voxels << " [" << gd[i].begin << " .. " << gd[i].end << "]\n";
}

inline vec grid_dims_begin(const grid_dims& gd) {
	vec tmp;
	VINA_FOR_IN(i, gd)
		tmp[i] = gd[i].begin;
	return tmp;
}

inline vec grid_dims_end(const grid_dims& gd) {
	vec tmp;
	VINA_FOR_IN(i, gd)
		tmp[i] = gd[i].end;
	return tmp;
}

#endif
