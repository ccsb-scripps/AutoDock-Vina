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

#include "current_weights.h"

flv current_weights(const terms& t) { // verifies size compatibility with t
	const fl a[] = {-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 1.923}; // design.out227
	//const fl a[] =   {-0.035579, -0.005156, 0.840245, -0.035069, 0, -0.587439, 0.3, 0}; // FROM design.out227 (FIXME 0.37)
	//const fl a[] = {-0.035579, -0.005156, 0.840245, -0.035069, 0, -0.587439, 1.923}; // from design.out227
	//const fl a[] = {-0.035579, -0.005156, 0.840245, -0.035069, 0, -0.587439, 0, 0}; // from design.out227
	//const fl a[] = {    -0.03729492,     -0.00447574,      0.49765303,     -0.01407780,      0.00604362,     -0.37693501,     -0.04244564,      0.41294104}; // some optimization
	//const fl a[]   = {    -0.06420507,     -0.00321022,      0.53200870,     -0.00858916,      0.01888147,     -0.43265646,     -0.15275852,      0.40187766};
	//const fl a[] = {    -0.26023247,      0.00995843,      3.22632287,     -1.03265955,     -2.17710460,     21.60339910}; // SOME DEBUGGING FIXME

	//              {-0.059, -0.002, 0.732, -0.045, -0.574, 1.022}
	//{    -0.04001575,     -0.00305351,      0.69244502,     -0.05569879,     -0.59632554,      1.30150720}
	//{    -0.05536152,     -0.00219298,      0.84076170,     -0.05453476,     -0.61922654,      1.00723685}
	//{    -0.04174947,     -0.00277397,      0.67397970,     -0.05836878,     -0.57969306,      1.13869379}
	//
	//{    -0.02008723,     -0.00412019,      0.60584644,     -0.05019072,     -0.58100053,      1.30362853}
	sz n = sizeof(a) / sizeof(const fl);
	flv tmp;
	VINA_FOR(i, n)
		tmp.push_back(a[i]);
	// FIXME rm
	sz tmp_size = tmp.size();
	sz conf_indep_size  = t.size_conf_independent(true);
	sz t_size = t.size();
	sz t_names_size = t.get_names(true).size();

	VINA_CHECK(tmp.size() == t.size_conf_independent(true) + t.get_names(true).size());
	return tmp;
}
