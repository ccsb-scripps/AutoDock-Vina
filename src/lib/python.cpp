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

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/to_python_converter.hpp>
#include "vina.h"
#include "parse_pdbqt.h"


void(Vina::*set_receptor_1)(const std::string&) = &Vina::set_receptor;
void(Vina::*set_receptor_2)(const std::string&, const std::string&) = &Vina::set_receptor;

void(Vina::*set_ligand_1)(const std::string&) = &Vina::set_ligand;
void(Vina::*set_ligand_2)(const std::vector<std::string>&) = &Vina::set_ligand;

BOOST_PYTHON_MODULE(vina)
{   
    using namespace boost::python;

    class_<Vina>("Vina")
        .def("global_search", &Vina::global_search)
        .def("set_receptor", set_receptor_1)
        .def("set_receptor", set_receptor_2)
        .def("set_ligand", set_ligand_1)
        .def("set_ligand", set_ligand_2)
        .def("set_box", &Vina::set_box)
        .def("set_weights", &Vina::set_weights)
        .def("randomize", &Vina::randomize)
        .def("optimize", &Vina::optimize)
        .def("score", &Vina::score)
        .def("compute_grid", &Vina::compute_grid)
        .def("set_forcefield", &Vina::set_forcefield)
        .def("write_results", &Vina::write_results)
        .def("write_pose", &Vina::write_pose)
    ;
}