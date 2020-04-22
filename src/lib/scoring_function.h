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

#ifndef VINA_SCORING_FUNCTION_H
#define VINA_SCORING_FUNCTION_H

#include "atom_type.h"

struct model; // forward declaration

struct scoring_function {
	virtual atom_type::t atom_typing_used() const = 0;
	virtual fl cutoff() const = 0;
	virtual fl eval(sz t1, sz t2, fl r) const = 0;
	virtual fl conf_independent(const model& m, fl e) const = 0;
};

#endif
