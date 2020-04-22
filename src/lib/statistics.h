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

#ifndef VINA_STATISTICS_H
#define VINA_STATISTICS_H

#include <algorithm> // sort
#include <cmath> // sqrt
#include "common.h" // for flv

// simple but inefficient implementations

inline fl mean(const flv& v) {
	fl acc = 0;
	VINA_FOR_IN(i, v) 
		acc += v[i];
	return v.empty() ? 0 : (acc/v.size());
}

inline fl deviation(const flv& v) {
	fl m = mean(v);
	fl acc = 0;
	VINA_FOR_IN(i, v)
		acc += sqr(v[i] - m);
	return v.empty() ? 0 : std::sqrt(acc/v.size());
}

inline fl rmsd(const flv& a, const flv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a)
		acc += sqr(a[i] - b[i]);
	return a.empty() ? 0 : std::sqrt(acc/a.size());
}

inline fl average_difference(const flv& b, const flv& a) { // b - a
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a)
		acc += b[i] - a[i];
	return a.empty() ? 0 : (acc/a.size());
}

inline fl pearson(const flv& x, const flv& y) {
	sz n = x.size();
	VINA_CHECK(n == y.size());
	if(n == 0) return 0; // ?
	fl sum_x = 0;
	fl sum_y = 0;
	fl sum_x_sq = 0;
	fl sum_y_sq = 0;
	fl sum_prod = 0;
	
	VINA_FOR(i, n) {
		sum_x += x[i];
		sum_y += y[i];
		sum_x_sq += sqr(x[i]);
		sum_y_sq += sqr(y[i]);
		sum_prod += x[i] * y[i];
	}
	fl sd_x = std::sqrt(sum_x_sq/n - sqr(sum_x/n)); // FIXME the argument is supposed to be >= 0, but ...
	fl sd_y = std::sqrt(sum_y_sq/n - sqr(sum_y/n));
	fl cov = sum_prod/n - (sum_x/n) * (sum_y/n);
	fl tmp = sd_x * sd_y;
	if(std::abs(tmp) < epsilon_fl) return 0; // ?
	return cov / tmp;
}

struct spearman_aux {
	fl x;
	sz i;
	spearman_aux(fl x, sz i) : x(x), i(i) {}
};

inline bool operator<(const spearman_aux& a, const spearman_aux& b) {
	return a.x < b.x;
}

inline flv get_rankings(const flv& x) {
	std::vector<spearman_aux> to_sort;
	VINA_FOR_IN(i, x)
		to_sort.push_back(spearman_aux(x[i], i));
	std::sort(to_sort.begin(), to_sort.end());

	flv tmp(x.size(), 0);
	VINA_FOR_IN(rank, to_sort)
		tmp[to_sort[rank].i] = rank;
	return tmp;
}

inline fl spearman(const flv& x, const flv& y) {
	return pearson(get_rankings(x), get_rankings(y));
}

#endif
