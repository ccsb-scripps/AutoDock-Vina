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

#ifndef VINA_CURL_H
#define VINA_CURL_H

#include "common.h"

#if 1 // prefer this to "hard curl"? 
template<typename T> // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
	if(e > 0 && not_max(v)) { // FIXME authentic_v can be gotten rid of everywhere now
		fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
		e *= tmp;
		deriv *= sqr(tmp);
	}
}

inline void curl(fl& e, fl v) {
	if(e > 0 && not_max(v)) {
		fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
		e *= tmp;
	}
}

#else

template<typename T> // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
	if(e > v) {
		e = v;
		deriv = 0;
	}
}

inline void curl(fl& e, fl v) {
	if(e > v) {
		e = v;
	}
}
#endif

#endif
