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

#include "terms.h"
#include "int_pow.h"

enum scoring_function_choice {SF_VINA, SF_AD42, SF_VINARDO};

struct everything : public terms {
	everything(scoring_function_choice sfchoice);
};

inline bool is_glued(sz xs_t1, sz xs_t2) {
    return (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_H_CG0) ||
       (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_P_CG0) ||
       (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_H_CG0) ||
       (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_P_CG0) ||

       (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_H_CG1) ||
       (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_P_CG1) ||
       (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_H_CG1) ||
       (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_P_CG1) ||

       (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_H_CG2) ||
       (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_P_CG2) ||
       (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_H_CG2) ||
       (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_P_CG2) ||

       (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_H_CG3) ||
       (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_P_CG3) ||
       (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_H_CG3) ||
       (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_P_CG3);
}
