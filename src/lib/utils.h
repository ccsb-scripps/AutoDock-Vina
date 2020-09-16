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

#ifndef VINA_UTILS_H
#define VINA_UTILS_H

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "tee.h"

inline char separator();
path make_path(const std::string& str);
void doing(int verbosity, const std::string& str, tee& log);
void done(int verbosity, tee& log);
std::string default_output(const std::string& input_name);
std::string default_output(const std::string& input_name, const std::string& directory_pathname);
bool is_directory(const std::string& directory_pathname);
std::string get_filename(const std::string& s);

#endif