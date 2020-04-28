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

#ifndef VINA_ARRAY3D_H
#define VINA_ARRAY3D_H

#include <exception> // std::bad_alloc
#include "common.h"

inline sz checked_multiply(sz i, sz j) {
	if(i == 0 || j == 0) return 0;
	const sz tmp = i * j;
	if(tmp < i || tmp < j || tmp / i != j)
		throw std::bad_alloc(); // can't alloc if the size makes sz wrap around
	return tmp;
}

inline sz checked_multiply(sz i, sz j, sz k) {
	return checked_multiply(checked_multiply(i, j), k);
}

template<typename T>
class array3d {
	sz m_i, m_j, m_k;
	std::vector<T> m_data;
	friend class boost::serialization::access;
	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & m_i;
		ar & m_j;
		ar & m_k;
		ar & m_data;
	}
public:
	array3d() : m_i(0), m_j(0), m_k(0) {}
	array3d(sz i, sz j, sz k) : m_i(i), m_j(j), m_k(k), m_data(checked_multiply(i, j, k)) {}
	sz dim0() const { return m_i; }
	sz dim1() const { return m_j; }
	sz dim2() const { return m_k; }
	sz dim(sz i) const {
		switch(i) {
			case 0: return m_i;
			case 1: return m_j;
			case 2: return m_k;
			default: assert(false); return 0; // to get rid of the warning
		}
	}
	void resize(sz i, sz j, sz k) { // data is essentially garbled
		m_i = i;
		m_j = j;
		m_k = k;
		m_data.resize(checked_multiply(i, j, k));
	}
	T&       operator()(sz i, sz j, sz k)       { return m_data[i + m_i*(j + m_j*k)]; }
	const T& operator()(sz i, sz j, sz k) const { return m_data[i + m_i*(j + m_j*k)]; }
};

#endif
