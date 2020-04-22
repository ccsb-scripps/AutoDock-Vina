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

#ifndef VINA_PARALLEL_PROGRESS_H
#define VINA_PARALLEL_PROGRESS_H

#include <boost/progress.hpp>
#include <boost/thread/mutex.hpp>

#include "incrementable.h"

struct parallel_progress : public incrementable {
	parallel_progress() : p(NULL) {}
	void init(unsigned long n) { p = new boost::progress_display(n); }
	void operator++() {
		if(p) {
			boost::mutex::scoped_lock self_lk(self);
			++(*p);
		}
	}
	virtual ~parallel_progress() { delete p; }
private:
	boost::mutex self;
	boost::progress_display* p;
};

#endif

