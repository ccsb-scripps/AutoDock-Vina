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

#include "forcefield.h"


std::vector<std::string> split(const std::string &s) {
    // Source: https://www.quora.com/How-do-I-split-a-string-by-space-into-an-array-in-c++
    std::vector<std::string> result;
    std::istringstream iss(s);

    for(std::string s; iss >> s; )
        result.push_back(s);

    return result;
}

bool ForceField::is_vina() {
    return m_is_vina;
}

bool ForceField::is_autodock4() {
    return m_is_autodock;
}

int ForceField::atom_type(const std::string& atom_type, const int hbond, const int bonded) {
    try {
        if (m_is_autodock) {
            // We just return the AD atom type int when it is AD
            return m_atom_types_to_int.at(atom_type);
        } else {
            // Vina atom types are based on the connectivity
            // We can have multiple vina atom types for one AD atom type
            std::vector<std::string> vina_atom_types = m_conversion_ad_vina.at(atom_type);

            for (auto& vina_atom_type : vina_atom_types) {
                // Get the atom type int of the vina type
                const int i = m_atom_types_to_int.at(vina_atom_type);
                // Return the first atom type int that match hbond and bonded parameters
                if ((m_hbond.at(i) == hbond) && (m_bonded.at(i) == bonded)) {
                    return i;
                }
            }
        }
    } catch (std::out_of_range& e) {
        std::cerr << "Error: Atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::weight(const std::string& weight_name) {
    try {
        return m_weights.at(weight_name);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: Weight parameter " << weight_name << " is not defined." << std::endl;
    }
}

double ForceField::rii_vdw(const int atom_type) {
    try {
        return m_rii_vdw.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: rii vdW for atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::epsii_vdw(const int atom_type) {
    try {
        return m_epsii_vdw.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: epsii vdW for atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::sol_par(const int atom_type) {
    try {
        return m_sol_par.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: solvation parameter for atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::vol(const int atom_type) {
    try {
        return m_vol.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: volume for atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::rij_hb(const int atom_type) {
    try {
        return m_rij_hb.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: rij hb for atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::epsij_hb(const int atom_type) {
    try {
        return m_epsij_hb.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: epsij hb for atom type " << atom_type << " is not defined." << std::endl;
    }
}

double ForceField::covalent_radius(const int atom_type) {
    const std::string current_element = element(atom_type);
    
    try {
        return m_covalent_radius.at(current_element);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: covalent radius for element " << current_element << " is not defined." << std::endl;
    }
}

const std::string ForceField::element(const int atom_type) {
    try {
        return m_element.at(atom_type);
    } catch (std::out_of_range& e) {
        std::cerr << "Error: element for atom type " << atom_type << " is not defined." << std::endl;
    }
}

bool ForceField::is_atom_type_defined(const std::string& atom_type) {
    if (m_atom_types_to_int.count(atom_type) > 0) {
        return true;
    } else {
        return false;
    }
}

bool ForceField::is_atom_type_defined(const int atom_type) {
    if (m_int_to_atom_types.count(atom_type) > 0) {
        return true;
    } else {
        return false;
    }
}

bool ForceField::is_hydrophobic(const int atom_type) {
    try {
        if (m_hbond.at(atom_type) == 0) {
            return true;
        } else {
            return false;
        }
    } catch (std::out_of_range& e) {
        std::cerr << "Error: hydrophobicity (m_hbond) for atom type " << atom_type << " is not defined." << std::endl;
    }

    return false;
}

bool ForceField::is_donor(const int atom_type) {
    try {
        if (m_hbond.at(atom_type) == 1 || m_hbond.at(atom_type) == 2) {
            return true;
        } else {
            return false;
        }
    } catch (std::out_of_range& e) {
        std::cerr << "Error: m_hbond for atom type " << atom_type << " is not defined." << std::endl;
    }

    return false;
}

bool ForceField::is_acceptor(const int atom_type) {
    try {
        if (m_hbond.at(atom_type) == 3 || m_hbond.at(atom_type) == 4 || m_hbond.at(atom_type) == 5) {
            return true;
        } else {
            return false;
        }
    } catch (std::out_of_range& e) {
        std::cerr << "Error: m_hbond for atom type " << atom_type << " is not defined." << std::endl;
    }

    return false;
}

bool ForceField::is_donor_acceptor(const int atom_type) {
    if (is_donor(atom_type) && is_acceptor(atom_type)) {
        return true;
    } else {
        return false;
    }
}

bool ForceField::is_hbond_possible(const int atom_type1, const int atom_type2) {
    if (is_hydrophobic(atom_type1) || is_hydrophobic(atom_type2)) {
        return false;
    }

    if ((is_donor(atom_type1) && is_acceptor(atom_type2)) || (is_acceptor(atom_type1) && is_donor(atom_type2))) {
        return true;
    } else {
        return false;
    }
}

bool ForceField::is_heteroatom(const int atom_type) {
    if (element(atom_type) != "H") {
        return true;
    } else {
        return false;
    }
}

bool ForceField::is_metal(const int atom_type) {
    std::vector<std::string> metals = {"Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"};
    std::string current_element = element(atom_type);

    if (std::count(metals.begin(), metals.end(), current_element)) {
        return true;
    } else {
        return false;
    }
}

double ForceField::optimal_covalent_bond_length(const int atom_type1, const int atom_type2) {
    return covalent_radius(atom_type1) + covalent_radius(atom_type2);
}

void ForceField::read_autodock_forcefield(const std::string& forcefield_filename) {
    int i = 0;
    std::string line;
    std::vector<std::string> results;
    std::ifstream myfile(forcefield_filename);

    if (myfile.is_open()) {
        while (getline (myfile, line)) {   
            if (line.rfind("FE_coeff_vdW", 0) == 0) {
                results = split(line);
                m_weights["vdw"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_hbond", 0) == 0) {
                results = split(line);
                m_weights["hbond"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_estat", 0) == 0) {
                results = split(line);
                m_weights["estat"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_desolv", 0) == 0) {
                results = split(line);
                m_weights["desolv"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_tors", 0) == 0) {
                results = split(line);
                m_weights["tors"] = std::stod(results[1]);
            } else if (line.rfind("atom_par", 0) == 0) {
                results = split(line);
                // For the atom_types from/to int conversion
                m_atom_types_to_int[results[1]] = i;
                m_int_to_atom_types[i] = results[1];
                // Collect atom type informations
                m_rii_vdw[i] = std::stod(results[2]);
                m_epsii_vdw[i] = std::stod(results[3]);
                m_vol[i] = std::stod(results[4]);
                m_sol_par[i] = std::stod(results[5]);
                m_rij_hb[i] = std::stod(results[6]);
                m_epsij_hb[i] = std::stod(results[7]);
                m_hbond[i] = std::stoi(results[8]);
                m_bond_index[i] = std::stoi(results[11]);
                m_element[i] = results[12];
                
                std::cout << i << " " << results[1] << " " << results[2] << " " << results[3] << "\n";

                i++;
            }
        }
        myfile.close();
    }
}

void ForceField::read_vina_forcefield(const std::string& forcefield_filename) {
    int i = 0;
    std::string line;
    std::vector<std::string> results;
    std::ifstream myfile(forcefield_filename);

    if (myfile.is_open()) {
        while (std::getline (myfile, line)) {   
            if (line.rfind("FE_coeff_gauss1", 0) == 0) {
                results = split(line);
                m_weights["gauss1"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_gauss2", 0) == 0) {
                results = split(line);
                m_weights["gauss2"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_repulsion", 0) == 0) {
                results = split(line);
                m_weights["repulsion"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_hydrophobic", 0) == 0) {
                results = split(line);
                m_weights["hydrophobic"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_hydrogen", 0) == 0) {
                results = split(line);
                m_weights["hydrogen"] = std::stod(results[1]);
            } else if (line.rfind("FE_coeff_rot", 0) == 0) {
                results = split(line);
                m_weights["rot"] = std::stod(results[1]);
            } else if (line.rfind("atom_par", 0) == 0) {
                results = split(line);
                // For the atom_types from/to int conversion
                m_atom_types_to_int[results[3]] = i;
                m_int_to_atom_types[i] = results[3];
                // We have more than one Vina atom types for some AD atom types
                m_conversion_ad_vina[results[1]].push_back(results[3]);
                // Collect atom type informations
                m_rii_vdw[i] = std::stod(results[3]);
                m_hbond[i] = std::stoi(results[4]);
                m_bonded[i] = std::stoi(results[5]);
                m_element[i] = results[6];

                std::cout << i << " " << results[1] << " " << results[3] << " " << results[4] << " " << results[5] << " " << results[6] << "\n";

                i++;
            }
        }
        myfile.close();
    }
}
