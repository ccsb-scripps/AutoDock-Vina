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

#ifndef VINA_PARSE_ERROR_H
#define VINA_PARSE_ERROR_H

#include "common.h"


class pdbqt_parse_error : public std::exception {
    public:    
        explicit pdbqt_parse_error(const std::string & message)
            : m_message("\n\nPDBQT parsing error: " + message + "\n") {}
        explicit pdbqt_parse_error(const std::string & message, const std::string & pdbqt_line)
            : m_message("\n\nPDBQT parsing error: " + message + "\n > " + pdbqt_line + "\n") {}

        virtual const char* what() const throw () {
            return m_message.c_str();
        }

    private:
        const std::string m_message;
};

#endif
