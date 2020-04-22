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

#ifndef VINA_TRIANGULAR_MATRIX_INDEX_H
#define VINA_TRIANGULAR_MATRIX_INDEX_H

#include "common.h"

inline sz triangular_matrix_index(sz n, sz i, sz j) {
	assert(j < n);
	assert(i <= j); 

	return i + j*(j+1)/2; 
}

inline sz triangular_matrix_index_permissive(sz n, sz i, sz j) {
	return (i <= j) ? triangular_matrix_index(n, i, j)
		            : triangular_matrix_index(n, j, i);
}

#endif
