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


//Forward declaration
struct model;

class precalculate_element
{
public:
    precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) { }
    fl eval_fast(fl r2) const;
    pr eval_deriv(fl r2) const;
    void init_from_smooth_fst(const flv& rs);
    sz min_smooth_fst() const;
    void widen_smooth_fst(const flv& rs, fl left, fl right);
    void widen(const flv& rs, fl left, fl right);

    prv smooth; // [(e, dor)]

private:
    flv fast;
    fl factor;
};

class precalculate 
{
public:
    precalculate() { }
    precalculate(const ScoringFunction& sf, fl v=max_fl, fl factor=32);
    fl eval_fast(sz type_pair_index, fl r2) const;
    pr eval_deriv(sz type_pair_index, fl r2) const;
    sz index_permissive(sz t1, sz t2) const { return m_data.index_permissive(t1, t2); }
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }
    void widen(fl left, fl right);
private:
    flv calculate_rs() const;

    fl m_cutoff_sqr;
    fl m_max_cutoff_sqr;
    sz m_n;
    fl m_factor;

    triangular_matrix<precalculate_element> m_data;
};

class precalculate_byatom 
{
public:
    precalculate_byatom() { }
    precalculate_byatom(const ScoringFunction &sf, const model &model, fl v=max_fl, fl factor=32);
    fl eval_fast(sz i, sz j, fl r2) const;
    pr eval_deriv(sz i, sz j, fl r2) const;
    sz index_permissive(sz t1, sz t2) const { return m_data.index_permissive(t1, t2); }
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }
    sz get_factor() const { return m_factor; }
    void widen(fl left, fl right);
private:
    flv calculate_rs() const;

    fl m_cutoff_sqr;
    fl m_max_cutoff_sqr;
    sz m_n;
    fl m_factor;

    triangular_matrix<precalculate_element> m_data;
};

#endif
