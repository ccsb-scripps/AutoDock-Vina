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

#ifndef VINA_BFGS_H
#define VINA_BFGS_H

#include "matrix.h"

typedef triangular_matrix<fl> flmat;

inline void minus_mat_vec_product(const flmat& m, const fl* in, fl* out, sz n) {
	VINA_FOR(i, n) {
		fl sum = 0;
		VINA_FOR(j, n)
			sum += m(m.index_permissive(i, j)) * in[j];
		out[i] = -sum;
	}
}

inline fl scalar_product(const fl* a, const fl* b, sz n) {
	fl tmp = 0;
	VINA_FOR(i, n)
		tmp += a[i] * b[i];
	return tmp;
}

inline bool bfgs_update(flmat& h, const fl* p, const fl* y, fl* minus_hy, const fl alpha, sz n) {
	const fl yp  = scalar_product(y, p, n);
	if(alpha * yp < epsilon_fl) return false; // FIXME?
	minus_mat_vec_product(h, y, minus_hy, n);
	const fl yhy = - scalar_product(y, minus_hy, n);
	const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
	VINA_FOR(i, n)
		VINA_RANGE(j, i, n) // includes i
			h(i, j) +=   alpha * r * (minus_hy[i] * p[j]
	                                + minus_hy[j] * p[i]) +
			           + alpha * alpha * (r*r * yhy  + r) * p[i] * p[j]; // s * s == alpha * alpha * p * p
	return true;
}

template<typename F, typename Conf, typename Change>
fl line_search(F& f, sz n, const Conf& x, const fl* g_flat, const fl f0, const Change& p, const fl* p_flat, Conf& x_new, Change& g_new, fl& f1, int& evalcount) { // returns alpha
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	const fl pg = scalar_product(p_flat, g_flat, n);

	VINA_U_FOR(trial, max_trials) {
		x_new = x; x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
		evalcount++;
		if(f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
			break;
		alpha *= multiplier;
	}
	return alpha;
}

inline void set_diagonal(flmat& m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

template<typename F, typename Conf, typename Change>
fl bfgs(F& f, Conf& x, Change& g, const unsigned max_steps, const fl average_required_improvement, const sz over,
		int& evalcount) { // x is I/O, final value is returned
	sz n = g.num_floats();
	flmat h(n, 0);
	set_diagonal(h, 1);

	Change g_new(g);
	Conf x_new(x);
	fl f0 = f(x, g);
	evalcount++;

	fl f_orig = f0;
	Change g_orig(g);
	Conf x_orig(x);

	Change p(g);

	flv f_values; f_values.reserve(max_steps+1);
	f_values.push_back(f0);

	// Flat copies of g, p, y and -H*y, allocated once per call
	flv flat(4 * n);
	fl* g_flat = n ? &flat[0] : NULL;
	fl* p_flat = g_flat + n;
	fl* y_flat = p_flat + n;
	fl* minus_hy_flat = y_flat + n;
	g.read_flat(g_flat);

	VINA_U_FOR(step, max_steps) {
		minus_mat_vec_product(h, g_flat, p_flat, n);
		p.write_flat(p_flat);
		fl f1 = 0;
		const fl alpha = line_search(f, n, x, g_flat, f0, p, p_flat, x_new, g_new, f1, evalcount);
		g_new.read_flat(y_flat);
		VINA_FOR(i, n)
			y_flat[i] -= g_flat[i]; // y = g_new - g

		f_values.push_back(f1);
		f0 = f1;
		x = x_new;
		if(!(std::sqrt(scalar_product(g_flat, g_flat, n)) >= 1e-5)) break; // breaks for nans too // FIXME !!??
		g = g_new; // ?
		g.read_flat(g_flat);

		if(step == 0) {
			const fl yy = scalar_product(y_flat, y_flat, n);
			if(std::abs(yy) > epsilon_fl)
				set_diagonal(h, alpha * scalar_product(y_flat, p_flat, n) / yy);
		}

		bool h_updated = bfgs_update(h, p_flat, y_flat, minus_hy_flat, alpha, n);
	}
	if(!(f0 <= f_orig)) { // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
	return f0;
}

#endif
