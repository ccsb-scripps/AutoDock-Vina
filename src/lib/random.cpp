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

#include <ctime> // for time (for seeding)

#include "random.h"

fl random_fl(fl a, fl b, rng& generator) { // expects a < b, returns rand in [a, b]
	assert(a < b); // BOOST also asserts a < b
	typedef boost::uniform_real<fl> distr;
	boost::variate_generator<rng&, distr> r(generator, distr(a, b));
	fl tmp = r();
	assert(tmp >= a);
	assert(tmp <= b);
	return tmp;
}

fl random_normal(fl mean, fl sigma, rng& generator) { // expects sigma >= 0
	assert(sigma >= 0); // BOOST asserts this as well
	typedef boost::normal_distribution<fl> distr;
	boost::variate_generator<rng&, distr> r(generator, distr(mean, sigma));
	return r();
}

int random_int(int a, int b, rng& generator) { // expects a <= b, returns rand in [a, b]
	assert(a <= b); // BOOST asserts this as well
	typedef boost::uniform_int<int> distr;
	boost::variate_generator<rng&, distr> r(generator, distr(a, b));
	int tmp = r();
	assert(tmp >= a);
	assert(tmp <= b);
	return tmp;
}

sz random_sz(sz a, sz b, rng& generator) { // expects a <= b, returns rand in [a, b]
	assert(a <= b);
	assert(int(a) >= 0);
	assert(int(b) >= 0);
	int i = random_int(int(a), int(b), generator);
	assert(i >= 0);
	assert(i >= int(a));
	assert(i <= int(b));
	return static_cast<sz>(i);
}

vec random_inside_sphere(rng& generator) {
	while(true) { // on average, this will have to be run about twice
		fl r1 = random_fl(-1, 1, generator);
		fl r2 = random_fl(-1, 1, generator);
		fl r3 = random_fl(-1, 1, generator);

		vec tmp(r1, r2, r3);
		if(sqr(tmp) < 1)
			return tmp;
	}
}

vec random_in_box(const vec& corner1, const vec& corner2, rng& generator) { // expects corner1[i] < corner2[i]
	vec tmp;
	VINA_FOR_IN(i, tmp)
		tmp[i] = random_fl(corner1[i], corner2[i], generator);
	return tmp;
}

int auto_seed() {
    // Seed generator, fix previous seed generator based on PID and time
    // Source: https://stackoverflow.com/questions/22883840/c-get-random-number-from-0-to-max-long-long-integer
    std::random_device rd; // Get a random seed from the OS entropy device, or whatever
    std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937 generator
                               // and seed it with entropy.

    // Define the distribution, by default it goes from 0 to MAX(unsigned int)
    // or what have you.
    std::uniform_int_distribution<unsigned int> distr;
    return distr(eng);
}
