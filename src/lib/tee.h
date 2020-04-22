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

#ifndef VINA_TEE_H
#define VINA_TEE_H

#include <iostream>
#include "file.h"

struct tee {
	ofile* of;
	tee() : of(NULL) {}
	void init(const path& name) {
		of = new ofile(name);
	}
	virtual ~tee() { delete of; }
	void flush() {
		std::cout << std::flush;
		if(of)
			(*of) << std::flush;
	}
	void endl() {
		std::cout << std::endl;
		if(of)
			(*of) << std::endl;
	}
	void setf(std::ios::fmtflags a) {
		std::cout.setf(a);
		if(of)
			of->setf(a);
	}
	void setf(std::ios::fmtflags a, std::ios::fmtflags b) {
		std::cout.setf(a, b);
		if(of)
			of->setf(a, b);
	}
};

template<typename T>
tee& operator<<(tee& out, const T& x) {
	std::cout << x;
	if(out.of)
		(*out.of) << x;
	return out;
}

#endif
