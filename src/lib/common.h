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

#ifndef VINA_COMMON_H
#define VINA_COMMON_H

#include <cassert>
#include <string>
#include <limits>
#include <utility> // pair
#include <algorithm> // too common
#include <vector> // used in typedef, and commonly used overall
#include <cmath> // commonly used
#include <iostream> // various debugging everywhere
#include <fstream> // print_coords
#include <iomanip> // to_string
#include <sstream> // to_string
#include <string> // probably included by the above anyway, common anyway

#include <boost/serialization/vector.hpp> // can't come before the above two - wart fixed in upcoming Boost versions
#include <boost/serialization/base_object.hpp> // movable_atom needs it - (derived from atom)
#include <boost/filesystem/path.hpp> // typedef'ed

#include "macros.h"

typedef double fl;

template<typename T>
T sqr(T x) {
	return x*x;
}

const fl not_a_num = std::sqrt(fl(-1)); // FIXME? check 

typedef std::size_t sz;
typedef std::pair<fl, fl> pr;

struct vec {
	fl data[3];
	vec() {
#ifndef NDEBUG
		data[0] = data[1] = data[2] = not_a_num;
#endif
	}
	vec(fl x, fl y, fl z) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	const fl& operator[](sz i) const { assert(i < 3); return data[i]; }
	      fl& operator[](sz i)       { assert(i < 3); return data[i]; }
    fl norm_sqr() const {
		return sqr(data[0]) + sqr(data[1]) + sqr(data[2]);
	}
	fl norm() const {
		return std::sqrt(norm_sqr());
	}
	fl operator*(const vec& v) const {
		return data[0] * v[0] + data[1] * v[1] + data[2] * v[2];
	}
	const vec& operator+=(const vec& v) {
		data[0] += v[0];
		data[1] += v[1];
		data[2] += v[2];
		return *this;
	}
	const vec& operator-=(const vec& v) {
		data[0] -= v[0];
		data[1] -= v[1];
		data[2] -= v[2];
		return *this;
	}
	const vec& operator+=(fl s) {
		data[0] += s;
		data[1] += s;
		data[2] += s;
		return *this;
	}
	const vec& operator-=(fl s) {
		data[0] -= s;
		data[1] -= s;
		data[2] -= s;
		return *this;
	}
	vec operator+(const vec& v) const {
		return vec(data[0] + v[0],
			       data[1] + v[1],
				   data[2] + v[2]);
	}
	vec operator-(const vec& v) const {
		return vec(data[0] - v[0],
			       data[1] - v[1],
				   data[2] - v[2]);
	}
	const vec& operator*=(fl s) {
		data[0] *= s;
		data[1] *= s;
		data[2] *= s;
		return *this;
	}
	void assign(fl s) {
		data[0] = data[1] = data[2] = s;
	}
	sz size() const { return 3; }
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & data;
	}
};

inline vec operator*(fl s, const vec& v) {
	return vec(s * v[0], s * v[1], s * v[2]);
}

inline vec cross_product(const vec& a, const vec& b) {
	return vec(a[1]*b[2] - a[2]*b[1],
		       a[2]*b[0] - a[0]*b[2],
		       a[0]*b[1] - a[1]*b[0]);
}

inline vec elementwise_product(const vec& a, const vec& b) {
	return vec(a[0] * b[0],
		       a[1] * b[1],
			   a[2] * b[2]);
}

struct mat {
	fl data[9];
	mat() {
#ifndef NDEBUG
		data[0] = data[1] = data[2] =
		data[3] = data[4] = data[5] =
		data[6] = data[7] = data[8] = not_a_num;
#endif
	}
	// column-major
	const fl& operator()(sz i, sz j) const { assert(i < 3); assert(j < 3); return data[i + 3*j]; }
	      fl& operator()(sz i, sz j)       { assert(i < 3); assert(j < 3); return data[i + 3*j]; }

	mat(fl xx, fl xy, fl xz,
		fl yx, fl yy, fl yz,
		fl zx, fl zy, fl zz) {

		data[0] = xx; data[3] = xy; data[6] = xz;
		data[1] = yx; data[4] = yy; data[7] = yz;
		data[2] = zx; data[5] = zy; data[8] = zz;
	}
	vec operator*(const vec& v) const {
		return vec(data[0]*v[0] + data[3]*v[1] + data[6]*v[2], 
			       data[1]*v[0] + data[4]*v[1] + data[7]*v[2],
				   data[2]*v[0] + data[5]*v[1] + data[8]*v[2]);
	}
	const mat& operator*=(fl s) {
		VINA_FOR(i, 9)
			data[i] *= s;
		return *this;
	}
};


typedef std::vector<vec> vecv;
typedef std::pair<vec, vec> vecp;
typedef std::vector<fl> flv;
typedef std::vector<pr> prv;
typedef std::vector<sz> szv;
typedef boost::filesystem::path path;

struct internal_error {
	std::string file;
	unsigned line;
	internal_error(const std::string& file_, unsigned line_) : file(file_), line(line_) {}
};

#ifdef NDEBUG
	#define VINA_CHECK(P) do { if(!(P)) throw internal_error(__FILE__, __LINE__); } while(false)
#else
	#define VINA_CHECK(P) assert(P)
#endif

const fl pi = fl(3.1415926535897931);

inline sz fl_to_sz(fl x, sz max_sz) { // return a value in [0, max_sz]
    if(x <= 0) return 0;
    if(x >= max_sz) return max_sz;
	sz tmp = static_cast<sz>(x);
	return (std::min)(tmp, max_sz); // sz -> fl cast loses precision. 'min' is to guard against returning values out of range
}

const fl fl_tolerance = fl(0.001);

inline bool eq(fl a, fl b) {
	return std::abs(a - b) < fl_tolerance; 
}

inline bool eq(const vec& a, const vec& b) {
	return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

template<typename T>
bool eq(const std::vector<T>& a, const std::vector<T>& b) {
	if(a.size() != b.size()) return false;
	VINA_FOR_IN(i, a)
		if(!eq(a[i], b[i]))
			return false;
	return true;
}

const fl max_fl = (std::numeric_limits<fl>::max)();
const sz max_sz = (std::numeric_limits<sz>::max)();
const unsigned max_unsigned = (std::numeric_limits<unsigned>::max)();
const fl epsilon_fl = std::numeric_limits<fl>::epsilon();

const vec zero_vec(0, 0, 0);
const vec max_vec(max_fl, max_fl, max_fl);

inline bool not_max(fl x) {
	return (x < 0.1 * max_fl);
}

inline fl vec_distance_sqr(const vec& a, const vec& b) {
	return sqr(a[0] - b[0]) + \
		   sqr(a[1] - b[1]) + \
		   sqr(a[2] - b[2]);
}

inline fl sqr(const vec& v) {
	return sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
}

template<typename T>
sz vector_append(std::vector<T>& x, const std::vector<T>& y) { // return old size
	sz old_size = x.size();
	x.insert(x.end(), y.begin(), y.end());
	return old_size;
}

template<typename T>
sz find_min(const std::vector<T>& v) { // returns v.size() i.e. 0 for empty vectors; the index of the smallest elt otherwise
	sz tmp = v.size();
	VINA_FOR_IN(i, v)
		if(i == 0 || v[i] < v[tmp])
			tmp = i;
	return tmp;
}

inline void normalize_angle(fl& x) { // subtract or add enough 2*pi's to make x be in [-pi, pi]
	if     (x >  3*pi) { // very large
		fl n = ( x - pi) / (2*pi); // how many 2*pi's do you want to subtract?
		x -= 2*pi*std::ceil(n); // ceil can be very slow, but this should not be called often
		normalize_angle(x);
	}
	else if(x < -3*pi) { // very small
		fl n = (-x - pi) / (2*pi); // how many 2*pi's do you want to add?
		x += 2*pi*std::ceil(n); // ceil can be very slow, but this should not be called often
		normalize_angle(x);
	}
	else if(x >    pi) { // in (   pi, 3*pi]
		x -= 2*pi;
	}
	else if(x <   -pi) { // in [-3*pi,  -pi)
		x += 2*pi;
	}
	assert(x >= -pi && x <= pi);
	// in [-pi, pi]
}

inline fl normalized_angle(fl x) {
	normalize_angle(x);
	return x;
}

template<typename T>
std::string to_string(const T& x, std::streamsize width = 0, char fill = ' ') { // default 0 width means no restrictions on width
	std::ostringstream out;
	out.fill(fill);
	if(width > 0)
		out << std::setw(width);
	out << x;
	return out.str();
}

template<typename T>
T sum(const std::vector<T>& v) {
	T acc = 0;
	VINA_FOR_IN(i, v)
		acc += v[i];
	return acc;
}

// multiply pK by this to get free energy in kcal/mol:
	// K = exp(E/RT)  -- lower K and lower E == better binder
	// pK = -log10(K)   => K = 10^(-pK)
	// E = RT ln(K) = RT ln (10^(-pK)) = - RT * ln(10) * pK
const fl pK_to_energy_factor = -8.31 /* RT in J/K/mol */ * 0.001 /* kilo */ * 300 /* K */ / 4.184 /* J/cal */ * std::log(10.0); //  -0.6 kcal/mol * log(10) = -1.38

inline fl pK_to_energy(fl pK) { return pK_to_energy_factor * pK; }

inline void print(fl x, std::ostream& out = std::cout) {
	out << x;
}

inline void print(sz x, std::ostream& out = std::cout) {
	out << x;
}

inline void print(const vec& v, std::ostream& out = std::cout) {
	out << "(";
	VINA_FOR_IN(i, v) {
		if(i != 0) 
			out << ", ";
		print(v[i], out);
	}
	out << ")";
}

template<typename T>
void print(const std::vector<T>& v, std::ostream& out = std::cout) {
	out << "[";
	VINA_FOR_IN(i, v) {
		if(i != 0) 
			out << " ";
		print(v[i], out);
	}
	out << "]";
}

template<typename T>
void printnl(const T& x, std::ostream& out = std::cout) {
	print(x, out);
	out << '\n';
}

inline bool starts_with(const std::string& str, const std::string& start) {
	return str.size() >= start.size() && str.substr(0, start.size()) == start;
}

template<typename T>
bool has(const std::vector<T>& v, const T& element) {
	return std::find(v.begin(), v.end(), element) != v.end();
}

#endif
