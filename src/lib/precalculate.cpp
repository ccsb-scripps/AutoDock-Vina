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

#include "precalculate.h"
#include "model.h"


fl precalculate_element::eval_fast(fl r2) const
{
    assert(r2 * factor < fast.size());
    sz i = sz(factor * r2); // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
    assert(i < fast.size());
    return fast[i];
}

pr precalculate_element::eval_deriv(fl r2) const
{
    fl r2_factored = factor * r2;
    assert(r2_factored + 1 < smooth.size());
    sz i1 = sz(r2_factored);
    sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
    assert(i1 < smooth.size());
    assert(i2 < smooth.size());
    fl rem = r2_factored - i1;
    assert(rem >= -epsilon_fl);
    assert(rem < 1 + epsilon_fl);
    const pr &p1 = smooth[i1];
    const pr &p2 = smooth[i2];
    fl e = p1.first + rem * (p2.first - p1.first);
    fl dor = p1.second + rem * (p2.second - p1.second);
    return pr(e, dor);
}

void precalculate_element::init_from_smooth_fst(const flv &rs)
{
    sz n = smooth.size();
    VINA_CHECK(rs.size() == n);
    VINA_CHECK(fast.size() == n);
    VINA_FOR(i, n) {
        // calculate dor's
        fl &dor = smooth[i].second;
        if (i == 0 || i == n - 1)
            dor = 0;
        else {
            fl delta = rs[i + 1] - rs[i - 1];
            fl r = rs[i];
            dor = (smooth[i + 1].first - smooth[i - 1].first) / (delta * r);
        }
        // calculate fast's from smooth.first's
        fl f1 = smooth[i].first;
        fl f2 = (i + 1 >= n) ? 0 : smooth[i + 1].first;
        fast[i] = (f2 + f1) / 2;
    }
}

sz precalculate_element:: min_smooth_fst() const
{
    sz tmp = 0; // returned if smooth.empty()
    VINA_FOR_IN(i_inv, smooth) {
        sz i = smooth.size() - i_inv - 1; // i_inv < smooth.size()  => i_inv + 1 <= smooth.size()
        if (i_inv == 0 || smooth[i].first < smooth[tmp].first)
            tmp = i;
    }
    return tmp;
}

void precalculate_element::widen_smooth_fst(const flv &rs, fl left, fl right)
{
    flv tmp(smooth.size(), 0); // the new smooth[].first's
    sz min_index = min_smooth_fst();
    VINA_CHECK(min_index < rs.size()); // won't hold for n == 0
    VINA_CHECK(rs.size() == smooth.size());
    fl optimal_r = rs[min_index];
    VINA_FOR_IN(i, smooth)
    {
        fl r = rs[i];
        if (r < optimal_r - left)
            r += left;
        else if (r > optimal_r + right)
            r -= right;
        else
            r = optimal_r;

        if (r < 0)
            r = 0;
        if (r > rs.back())
            r = rs.back();

        tmp[i] = eval_deriv(sqr(r)).first;
    }
    VINA_FOR_IN(i, smooth)
    smooth[i].first = tmp[i];
}

void precalculate_element::widen(const flv &rs, fl left, fl right)
{
    widen_smooth_fst(rs, left, right);
    init_from_smooth_fst(rs);
}

precalculate::precalculate(const ScoringFunction &sf, fl v, fl factor)
{   // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
    m_factor = factor;
    m_cutoff_sqr = sqr(sf.get_max_cutoff());
    m_n = sz(m_factor * m_cutoff_sqr) + 3; // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
    //std::cout << "-- DEBUG -- sf.cutoff^2 in precalculate = " << m_cutoff_sqr << "\n";
    triangular_matrix<precalculate_element> data(num_atom_types(sf.get_atom_typing()), precalculate_element(m_n, m_factor));

    VINA_CHECK(m_factor > epsilon_fl);
    VINA_CHECK(sz(m_cutoff_sqr * m_factor) + 1 < m_n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
    VINA_CHECK(m_cutoff_sqr * m_factor + 1 < m_n);

    flv rs = calculate_rs();

    VINA_FOR(t1, data.dim())
    {
        VINA_RANGE(t2, t1, data.dim())
        {
            precalculate_element &p = data(t1, t2);
            // init smooth[].first
            VINA_FOR_IN(i, p.smooth)
            {
                p.smooth[i].first = (std::min)(v, sf.eval(t1, t2, rs[i]));
            }

            // init the rest
            p.init_from_smooth_fst(rs);
        }
    }

    m_data = data;
}

fl precalculate::eval_fast(sz type_pair_index, fl r2) const
{
    assert(r2 <= m_cutoff_sqr);
    return m_data(type_pair_index).eval_fast(r2);
}

pr precalculate::eval_deriv(sz type_pair_index, fl r2) const
{
    assert(r2 <= m_cutoff_sqr);
    return m_data(type_pair_index).eval_deriv(r2);
}

void precalculate::widen(fl left, fl right)
{
    flv rs = calculate_rs();
    VINA_FOR(t1, m_data.dim())
    VINA_RANGE(t2, t1, m_data.dim())
    m_data(t1, t2).widen(rs, left, right);
}

flv precalculate::calculate_rs() const
{
    flv tmp(m_n, 0);
    VINA_FOR(i, m_n)
    tmp[i] = std::sqrt(i / m_factor);
    return tmp;
}

precalculate_byatom::precalculate_byatom(const ScoringFunction &sf, const model &model, fl v, fl factor)
{   // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
    m_factor = factor;
    m_cutoff_sqr = sqr(sf.get_max_cutoff());
    m_n = sz(m_factor * m_cutoff_sqr) + 3; // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
    //std::cout << "-- DEBUG -- sf.cutoff^2 in precalculate = " << m_cutoff_sqr << "\n";
    sz n_atoms = model.num_movable_atoms();
    atomv atoms = model.get_atoms();
    triangular_matrix<precalculate_element> data(n_atoms, precalculate_element(m_n, m_factor));

    VINA_CHECK(m_factor > epsilon_fl);
    VINA_CHECK(sz(m_cutoff_sqr * m_factor) + 1 < m_n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
    VINA_CHECK(m_cutoff_sqr * m_factor + 1 < m_n);

    flv rs = calculate_rs();

    VINA_FOR(i, data.dim())
    {
        VINA_RANGE(j, i, data.dim())
        {
            precalculate_element &p = data(i, j);
            // init smooth[].first
            VINA_FOR_IN(k, p.smooth)
            {
                p.smooth[k].first = (std::min)(v, sf.eval(atoms[i], atoms[j], rs[k]));
            }

            // init the rest
            p.init_from_smooth_fst(rs);
        }
    }

    m_data = data;
}

fl precalculate_byatom::eval_fast(sz i, sz j, fl r2) const
{
    assert(r2 <= m_cutoff_sqr);
    return m_data(i, j).eval_fast(r2);
}

pr precalculate_byatom::eval_deriv(sz i, sz j, fl r2) const
{
    assert(r2 <= m_cutoff_sqr);
    return m_data(i, j).eval_deriv(r2);
}

void precalculate_byatom::widen(fl left, fl right)
{
    flv rs = calculate_rs();
    VINA_FOR(t1, m_data.dim())
    VINA_RANGE(t2, t1, m_data.dim())
    m_data(t1, t2).widen(rs, left, right);
}

flv precalculate_byatom::calculate_rs() const
{
    flv tmp(m_n, 0);
    VINA_FOR(i, m_n)
    tmp[i] = std::sqrt(i / m_factor);
    return tmp;
}
