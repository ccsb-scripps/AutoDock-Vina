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

#include "utils.h"

path make_path(const std::string& str) {
    boost::filesystem::path p(str);
  return p;
}


void doing(int verbosity, const std::string& str, tee& log) {
  if(verbosity > 1) {
    log << str << std::string(" ... ");
    log.flush();
  }
}


void done(int verbosity, tee& log) {
  if(verbosity > 1) {
    log << "done.";
    log.endl();
  }
}


std::string default_output(const std::string& input_name) {
  std::string tmp = input_name;
  if(tmp.size() >= 6 && tmp.substr(tmp.size()-6, 6) == ".pdbqt")
    tmp.resize(tmp.size() - 6); // FIXME?
  return tmp + "_out.pdbqt";
}
