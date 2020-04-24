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

#include "quaternion.h"

bool quaternion_is_normalized(const qt& q) { // not in the interface, used in assertions
	return eq(quaternion_norm_sqr(q), 1) && eq(boost::math::abs(q), 1);
}

bool eq(const qt& a, const qt& b) { // elementwise approximate equality - may return false for equivalent rotations
	return eq(a.R_component_1(), b.R_component_1()) && \
		   eq(a.R_component_2(), b.R_component_2()) && \
		   eq(a.R_component_3(), b.R_component_3()) && \
		   eq(a.R_component_4(), b.R_component_4());
}

qt angle_to_quaternion(const vec& axis, fl angle) { // axis is assumed to be a unit vector
	//assert(eq(tvmet::norm2(axis), 1));
	assert(eq(axis.norm(), 1));
	normalize_angle(angle); // this is probably only necessary if angles can be very big
	fl c = std::cos(angle/2);
	fl s = std::sin(angle/2);
	return qt(c, s*axis[0], s*axis[1], s*axis[2]);
}

qt angle_to_quaternion(const vec& rotation) {
	//fl angle = tvmet::norm2(rotation); 
	fl angle = rotation.norm(); 
	if(angle > epsilon_fl) {
		//vec axis; 
		//axis = rotation / angle;	
		vec axis = (1/angle) * rotation;
		return angle_to_quaternion(axis, angle);
	}
	return qt_identity;
}

vec quaternion_to_angle(const qt& q) {
	assert(quaternion_is_normalized(q));
	const fl c = q.R_component_1();
	if(c > -1 && c < 1) { // c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
		fl angle = 2*std::acos(c); // acos is in [0, pi]
		if(angle > pi)
			angle -= 2*pi; // now angle is in [-pi, pi]
		vec axis(q.R_component_2(), q.R_component_3(), q.R_component_4());
		fl s = std::sin(angle/2); // perhaps not very efficient to calculate sin of acos
		if(std::abs(s) < epsilon_fl)
			return zero_vec;
		axis *= (angle / s);
		return axis;
	}
	else // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
		return zero_vec;
}

mat quaternion_to_r3(const qt& q) {
	assert(quaternion_is_normalized(q));

	const fl a = q.R_component_1();
	const fl b = q.R_component_2();
	const fl c = q.R_component_3();
	const fl d = q.R_component_4();

	const fl aa = a*a;
	const fl ab = a*b;
	const fl ac = a*c;
	const fl ad = a*d;
	const fl bb = b*b;
	const fl bc = b*c;
	const fl bd = b*d;
	const fl cc = c*c;
	const fl cd = c*d;
	const fl dd = d*d;

	assert(eq(aa+bb+cc+dd, 1));

	mat tmp;

	// from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
	tmp(0, 0) = (aa+bb-cc-dd);
	tmp(0, 1) = 2*(-ad+bc);
	tmp(0, 2) = 2*(ac+bd);

	tmp(1, 0) = 2*(ad+bc);
	tmp(1, 1) = (aa-bb+cc-dd);
	tmp(1, 2) = 2*(-ab+cd);

	tmp(2, 0) = 2*(-ac+bd);
	tmp(2, 1) = 2*(ab+cd);
	tmp(2, 2) = (aa-bb-cc+dd);

	return tmp;
}

qt random_orientation(rng& generator) {
	qt q(random_normal(0, 1, generator), 
		 random_normal(0, 1, generator), 
		 random_normal(0, 1, generator), 
		 random_normal(0, 1, generator));
	fl nrm = boost::math::abs(q);
	if(nrm > epsilon_fl) {
		q /= nrm;
		assert(quaternion_is_normalized(q));
		return q;
	}
	else 
		return random_orientation(generator); // this call should almost never happen
}

void quaternion_increment(qt& q, const vec& rotation) {
	assert(quaternion_is_normalized(q));
	q = angle_to_quaternion(rotation) * q;
	quaternion_normalize_approx(q); // normalization added in 1.1.2
	//quaternion_normalize(q); // normalization added in 1.1.2
}

vec quaternion_difference(const qt& b, const qt& a) { // rotation that needs to be applied to convert a to b
	quaternion_is_normalized(a);
	quaternion_is_normalized(b);
	qt tmp = b;
	tmp /= a; // b = tmp * a    =>   b * inv(a) = tmp 
	return quaternion_to_angle(tmp); // already assert normalization
}

void print(const qt& q, std::ostream& out) { // print as an angle
	print(quaternion_to_angle(q), out);
}
