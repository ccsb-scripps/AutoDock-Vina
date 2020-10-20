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

#ifndef VINA_FORCEFIELD_H
#define VINA_FORCEFIELD_H


#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <map>
#include <regex>


class ForceField {
public:
    ForceField(const std::string& forcefield_filename=std::string(), const std::string& forcefield_type="vina") {
        if (!forcefield_filename.empty() && forcefield_type == "vina") {
            read_vina_forcefield(forcefield_filename);
            m_is_vina = true;
        } else if (!forcefield_filename.empty() && forcefield_type == "AD4") {
            read_autodock_forcefield(forcefield_filename);
            m_is_autodock = true;
        } else if (!forcefield_filename.empty()) {
            std::cout << "Error: forcefield " << forcefield_type << " type not recognized. \n";
            std::cout << "ForceField available: vina or AD4\n";
            exit (EXIT_FAILURE);
        }

        m_covalent_radius = {
            {"C", 0.77}, {"N", 0.75}, {"O", 0.73}, {"P", 1.06},
            {"S", 1.02}, {"H", 0.37}, {"F", 0.71}, {"I", 1.33},
            {"Mg", 1.30}, {"Mn", 1.39}, {"Zn", 1.31}, {"Ca", 1.74},
            {"Fe", 1.25}, {"Cl", 0.99}, {"Br", 1.14}
        };
    }

    // Destructor
    virtual ~ForceField() { };

    bool is_vina();
    bool is_autodock4();

    int atom_type(const std::string& atom_type, const int hbond=0, const int bonded=0);
    double weight(const std::string& weight_name);
    double rii_vdw(const int atom_type);
    double epsii_vdw(const int atom_type);
    double sol_par(const int atom_type);
    double vol(const int atom_type);
    double rij_hb(const int atom_type);
    double epsij_hb(const int atom_type);
    double covalent_radius(const int atom_type);
    const std::string element(const int atom_type);
    bool is_atom_type_defined(const std::string& atom_type);
    bool is_atom_type_defined(const int atom_type);
    bool is_hydrogen(const int atom_type);
    bool is_hydrophobic(const int atom_type);
    bool is_heteroatom(const int atom_type);
    bool is_metal(const int atom_type);
    bool is_donor(const int atom_type);
    bool is_acceptor(const int atom_type);
    bool is_donor_acceptor(const int atom_type);
    bool is_hbond_possible(const int atom_type1, const int atom_type2);
    double optimal_covalent_bond_length(const int atom_type1, const int atom_type2);

private:
    bool m_is_vina = false;
    bool m_is_autodock = false;
    std::map<std::string, int> m_atom_types_to_int;
    std::map<int, std::string> m_int_to_atom_types;
    std::map<std::string, std::vector<std::string>> m_conversion_ad_vina;

    std::map<std::string, double> m_weights;
    std::map<std::string, double> m_covalent_radius;
    std::map<int, double> m_rii_vdw;
    std::map<int, double> m_epsii_vdw;
    std::map<int, double> m_sol_par;
    std::map<int, double> m_vol;
    std::map<int, double> m_rij_hb;
    std::map<int, double> m_epsij_hb;
    std::map<int, int> m_hbond;
    std::map<int, int> m_bond_index;
    std::map<int, int> m_bonded;
    std::map<int, std::string> m_element;

    void read_autodock_forcefield(const std::string& forcefield_filename);
    void read_vina_forcefield(const std::string& forcefield_filename);
};

#endif
