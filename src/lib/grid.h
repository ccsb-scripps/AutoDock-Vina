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

#ifndef VINA_GRID_H
#define VINA_GRID_H

#include "array3d.h"
#include "grid_dim.h"
#include "curl.h"

class grid { // FIXME rm 'm_', consistent with my new style
    vec m_init;
    vec m_range;
    vec m_factor;
    vec m_dim_fl_minus_1;
	vec m_factor_inv;
public:
	array3d<fl> m_data; // FIXME? - make cache a friend, and convert this back to private?
	grid() : m_init(0, 0, 0), m_range(1, 1, 1), m_factor(1, 1, 1), m_dim_fl_minus_1(-1, -1, -1), m_factor_inv(1, 1, 1) {} // not private
	grid(const grid_dims& gd) { init(gd); }
    void init(const grid_dims& gd);
	vec index_to_argument(sz x, sz y, sz z) const {
		return vec(m_init[0] + m_factor_inv[0] * x,
		           m_init[1] + m_factor_inv[1] * y,
		           m_init[2] + m_factor_inv[2] * z);
	}
	bool initialized() const {
		return m_data.dim0() > 0 && m_data.dim1() > 0 && m_data.dim2() > 0;
	}
	fl evaluate(const vec& location, fl slope, fl c)             const { return evaluate_aux(location, slope, c, NULL);   }
	fl evaluate(const vec& location, fl slope, fl c, vec& deriv) const { return evaluate_aux(location, slope, c, &deriv); } // sets deriv
private:
	fl evaluate_aux(const vec& location, fl slope, fl v, vec* deriv) const; // sets *deriv if not NULL
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & m_init;
		ar & m_data;
		ar & m_range;
		ar & m_factor;
		ar & m_dim_fl_minus_1;
		ar & m_factor_inv;
	}
};

#endif
