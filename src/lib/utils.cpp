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


inline char separator() {
    // Source: https://stackoverflow.com/questions/12971499/how-to-get-the-file-separator-symbol-in-standard-c-c-or
    #ifdef _WIN32
        return '\\';
    #else
        return '/';
    #endif
}


path make_path(const std::string& str) {
    boost::filesystem::path p(str);
  return p;
}


void doing(const std::string& str, int verbosity, int level) {
    if(verbosity > level) {
      std::cout << str << std::string(" ... ") << std::flush;
    }
}


void done(int verbosity, int level) {
    if(verbosity > level) {
        std::cout << "done.\n" << std::flush;
    }
}


std::string default_output(const std::string& input_name) {
    std::string tmp = input_name;
    if (tmp.size() >= 6 && tmp.substr(tmp.size()-6, 6) == ".pdbqt")
        tmp.resize(tmp.size() - 6); // FIXME?
    return tmp + "_out.pdbqt";
}


std::string default_output(const std::string& input_name, const std::string& directory_pathname) {
    return directory_pathname + separator() + default_output(input_name);
}


bool is_directory(const std::string& directory_pathname) {
  //Source: https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
  struct stat info;

    if (stat(directory_pathname.c_str(), &info) != 0) 
        return false;
    else if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows 
        return true;
    else
        return false;
}


std::string get_filename(const std::string& s) {
    size_t i = s.rfind(separator(), s.length());

    if (i != std::string::npos) {
        return(s.substr(i + 1, s.length() - i));
    }

    return(s);
}

std::string get_file_contents(const std::string& filename) {   
    ifile in(make_path(filename));

    if (in) {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return (contents);
    }
    throw file_error(filename, true);
}
