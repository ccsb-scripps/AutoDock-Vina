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

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include "scoring_function.h"
#include "matrix.h"

struct precalculate_element {
	precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) {}
	fl eval_fast(fl r2) const {
		assert(r2 * factor < fast.size());
		sz i = sz(factor * r2);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i < fast.size()); 
		return fast[i];
	}
	pr eval_deriv(fl r2) const {
		fl r2_factored = factor * r2;
		assert(r2_factored + 1 < smooth.size());
		sz i1 = sz(r2_factored); 
		sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i1 < smooth.size());
		assert(i2 < smooth.size());
		fl rem = r2_factored - i1;
		assert(rem >= -epsilon_fl);
		assert(rem < 1 + epsilon_fl);
		const pr& p1 = smooth[i1];
		const pr& p2 = smooth[i2];
		fl e   = p1.first  + rem * (p2.first  - p1.first);
		fl dor = p1.second + rem * (p2.second - p1.second); 
		return pr(e, dor);
	}
	void init_from_smooth_fst(const flv& rs) {
		sz n = smooth.size();
		VINA_CHECK(rs.size() == n);
		VINA_CHECK(fast.size() == n);
		VINA_FOR(i, n) {
			// calculate dor's
			fl& dor = smooth[i].second;
			if(i == 0 || i == n-1)
				dor = 0;
			else {
				fl delta = rs[i+1] - rs[i-1];
				fl r = rs[i];
				dor = (smooth[i+1].first - smooth[i-1].first) / (delta * r);
			}
			// calculate fast's from smooth.first's
			fl f1 = smooth[i].first;
			fl f2 = (i+1 >= n) ? 0 : smooth[i+1].first;
			fast[i] = (f2 + f1) / 2;
		}
	}
	sz min_smooth_fst() const {
		sz tmp = 0; // returned if smooth.empty()
		VINA_FOR_IN(i_inv, smooth) {
			sz i = smooth.size() - i_inv - 1; // i_inv < smooth.size()  => i_inv + 1 <= smooth.size()
			if(i_inv == 0 || smooth[i].first < smooth[tmp].first)
				tmp = i;
		}
		return tmp;
	}
	void widen_smooth_fst(const flv& rs, fl left, fl right) {
		flv tmp(smooth.size(), 0); // the new smooth[].first's
		sz min_index = min_smooth_fst();
		VINA_CHECK(min_index < rs.size()); // won't hold for n == 0
		VINA_CHECK(rs.size() == smooth.size());
		fl optimal_r   = rs[min_index];
		VINA_FOR_IN(i, smooth) {
			fl r = rs[i];
			if     (r < optimal_r - left ) r += left;
			else if(r > optimal_r + right) r -= right;
			else                           r = optimal_r;

			if(r < 0) r = 0;
			if(r > rs.back()) r = rs.back();

			tmp[i] = eval_deriv(sqr(r)).first;
		}
		VINA_FOR_IN(i, smooth)
			smooth[i].first = tmp[i];
	}
	void widen(const flv& rs, fl left, fl right) {
		widen_smooth_fst(rs, left, right);
		init_from_smooth_fst(rs);
	}
	flv fast;
	prv smooth; // [(e, dor)]
	fl factor;
};

struct precalculate {
	precalculate(const scoring_function& sf, fl v = max_fl, fl factor_ = 32) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
		m_cutoff_sqr(sqr(sf.cutoff())),
		n(sz(factor_ * m_cutoff_sqr) + 3),  // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
		factor(factor_),

		data(num_atom_types(sf.atom_typing_used()), precalculate_element(n, factor_)),
		m_atom_typing_used(sf.atom_typing_used()) {

		VINA_CHECK(factor > epsilon_fl);
		VINA_CHECK(sz(m_cutoff_sqr*factor) + 1 < n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
		VINA_CHECK(m_cutoff_sqr*factor + 1 < n);

		flv rs = calculate_rs();

		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim()) {
				precalculate_element& p = data(t1, t2);
				// init smooth[].first
				VINA_FOR_IN(i, p.smooth)
					p.smooth[i].first = (std::min)(v, sf.eval(t1, t2, rs[i]));

				// init the rest
				p.init_from_smooth_fst(rs);
			}
	}
	fl eval_fast(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_fast(r2);
	}
	pr eval_deriv(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_deriv(r2);
	}
	sz index_permissive(sz t1, sz t2) const { return data.index_permissive(t1, t2); }
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }
	fl cutoff_sqr() const { return m_cutoff_sqr; }
	void widen(fl left, fl right) {
		flv rs = calculate_rs();
		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim())
				data(t1, t2).widen(rs, left, right);
	}
private:
	flv calculate_rs() const {
		flv tmp(n, 0);
		VINA_FOR(i, n)
			tmp[i] = std::sqrt(i / factor);
		return tmp;
	}
	fl m_cutoff_sqr;
	sz n;
	fl factor;
	atom_type::t m_atom_typing_used;

	triangular_matrix<precalculate_element> data;
};

#endif
