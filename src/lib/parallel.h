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

#ifndef VINA_PARALLEL_H
#define VINA_PARALLEL_H

#include <vector>

#include "common.h"

#include <boost/optional.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>


template<typename F, bool Sync = false>
struct parallel_for : private boost::thread_group {
	parallel_for(const F* f, sz num_threads) : m_f(f), destructing(false), size(0), thread_finished(num_threads, true), count_finished(0), num_threads(num_threads) {
        VINA_FOR(i, num_threads)
            create_thread(aux(i, this));
    }
	void run(sz size_) {
		boost::mutex::scoped_lock self_lk(self);
        size = size_;
        count_finished = 0;
		VINA_FOR_IN(i, thread_finished) 
			thread_finished[i] = false;
        cond.notify_all(); // many things modified
        while(count_finished < num_threads) // wait until processing of all elements is thread_finished
            busy.wait(self_lk);
    }
    virtual ~parallel_for() {
        {
			boost::mutex::scoped_lock self_lk(self);
            destructing = true;
            cond.notify_all(); // destructing modified
        }
        join_all(); 
    }
private:
	void loop(sz offset) {
		while(boost::optional<sz> sz_option = get_size(offset)) {
			sz s = sz_option.get();
			for(sz i = offset; i < s; i += num_threads)
				(*m_f)(i);
			{
				boost::mutex::scoped_lock self_lk(self);
				thread_finished[offset] = true;
				++count_finished;
				busy.notify_one();
			}
		}
	}
    struct aux {
		sz offset;
        parallel_for* par;
		aux(sz offset, parallel_for* par) : offset(offset), par(par) {}
		void operator()() const { par->loop(offset); }
    };
    const F* m_f; // does not keep a local copy!
    boost::condition cond;
    boost::condition busy;
    bool destructing; // dtor called
    sz size; // size of the vector given to run()
	std::vector<bool> thread_finished;
    sz count_finished; 
	sz num_threads;
	boost::mutex self; // any modification or reading of mutables should lock this first
	boost::optional<sz> get_size(sz offset) {
		boost::mutex::scoped_lock self_lk(self);
        while(!destructing && thread_finished[offset])
            cond.wait(self_lk);
		if(destructing) return boost::optional<sz>(); // wrap it up!
        return size;
    }
};

template<typename F>
struct parallel_for<F, true> : private boost::thread_group {
	parallel_for(const F* f, sz num_threads) : m_f(f), destructing(false), size(0), started(0), finished(0) {
		a.par = this; // VC8 warning workaround
        VINA_FOR(i, num_threads)
            create_thread(boost::ref(a));
    }
	void run(sz size_) {
		boost::mutex::scoped_lock self_lk(self);
        size = size_;
        finished = 0;
        started = 0;
        cond.notify_all(); // many things modified
        while(finished < size) // wait until processing of all elements is finished
            busy.wait(self_lk);
    }
    virtual ~parallel_for() {
        {
			boost::mutex::scoped_lock self_lk(self);
            destructing = true;
            cond.notify_all(); // destructing modified
        }
        join_all(); 
    }
private:
	void loop() {
		while(boost::optional<sz> i = get_next()) {
			(*m_f)(i.get());
			{
				boost::mutex::scoped_lock self_lk(self);
				++finished;
				busy.notify_one();
			}
		}
	}
    struct aux {
        parallel_for* par;
        aux() : par(NULL) {}
		void operator()() const { par->loop(); }
    };
    aux a;
    const F* m_f; // does not keep a local copy!
    boost::condition cond;
    boost::condition busy;
    bool destructing; // dtor called
    sz size; // size of the vector given to run() // FIXME?
    sz started; // the number of jobs given to run() the work started on
    sz finished; // the number of jobs given to run() the work finished on
	boost::mutex self; // any modification or reading of mutables should lock this first
	boost::optional<sz> get_next() {
		boost::mutex::scoped_lock self_lk(self);
        while(!destructing && started >= size)
            cond.wait(self_lk);
		if(destructing) return boost::optional<sz>(); // NOTHING
		sz tmp = started;
        ++started;
        return tmp;
    }
};


template<typename F, typename Container, typename Input, bool Sync = false>
struct parallel_iter { 
	parallel_iter(const F* f, sz num_threads) : a(f), pf(&a, num_threads) {}
	void run(Container& v) {
		a.v = &v;
		pf.run(v.size());
	}
private:
	struct aux {
		const F* f;
		Container* v;
		aux(const F* f) : f(f), v(NULL) {}
		void operator()(sz i) const { 
			VINA_CHECK(v);
			(*f)((*v)[i]); 
		}
	};
	aux a;
	parallel_for<aux, Sync> pf;
};

#endif
