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

#ifndef VINA_QUATERNION_H
#define VINA_QUATERNION_H

#include <boost/math/quaternion.hpp>
#include <boost/serialization/split_free.hpp>

#include "common.h"
#include "random.h"

typedef boost::math::quaternion<fl> qt;

// non-intrusive free function split serialization
namespace boost {
	namespace serialization {
		template<class Archive>
		void save(Archive& ar, const qt& q, const unsigned version) {
			fl q1 = q.R_component_1();
			fl q2 = q.R_component_2();
			fl q3 = q.R_component_3();
			fl q4 = q.R_component_4();

			ar & q1;
			ar & q2;
			ar & q3;
			ar & q4;
		}
		template<typename Archive>
		void load(Archive& ar, qt& q, const unsigned version) {
			fl a, b, c, d;
			ar & a;
			ar & b;
			ar & c;
			ar & d;
			q = qt(a, b, c, d);
		}
	}
}
BOOST_SERIALIZATION_SPLIT_FREE(qt)

bool eq(const qt& a, const qt& b); // elementwise approximate equality - may return false for equivalent rotations
const qt qt_identity(1, 0, 0, 0);
qt angle_to_quaternion(const vec& axis, fl angle); // axis is assumed to be a unit vector
qt angle_to_quaternion(const vec& rotation); // rotation == angle * axis
vec quaternion_to_angle(const qt& q);
mat quaternion_to_r3(const qt& q);

inline fl quaternion_norm_sqr(const qt& q) { // equivalent to sqr(boost::math::abs(const qt&))
	return sqr(q.R_component_1()) + sqr(q.R_component_2()) + sqr(q.R_component_3()) + sqr(q.R_component_4());
}

inline void quaternion_normalize(qt& q) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs(q))));
    const fl a = std::sqrt(s);
	assert(a > epsilon_fl);
	q *= 1/a;
	assert(quaternion_is_normalized(q));
}

inline void quaternion_normalize_approx(qt& q, const fl tolerance = 1e-6) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs(q))));
    if(std::abs(s - 1) < tolerance)
        ; // most likely scenario
    else {
        const fl a = std::sqrt(s);
        assert(a > epsilon_fl);
        q *= 1/a;
        assert(quaternion_is_normalized(q));
    }
}

qt random_orientation(rng& generator);
void quaternion_increment(qt& q, const vec& rotation);
vec quaternion_difference(const qt& b, const qt& a); // rotation that needs to be applied to convert a to b
void print(const qt& q, std::ostream& out = std::cout); // print as an angle

#endif
