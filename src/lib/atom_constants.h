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

#ifndef VINA_ATOM_CONSTANTS_H
#define VINA_ATOM_CONSTANTS_H

#include "common.h"

// based on SY_TYPE_* but includes H
const sz EL_TYPE_H    =  0;
const sz EL_TYPE_C    =  1;
const sz EL_TYPE_N    =  2;
const sz EL_TYPE_O    =  3;
const sz EL_TYPE_S    =  4;
const sz EL_TYPE_P    =  5;
const sz EL_TYPE_F    =  6;
const sz EL_TYPE_Cl   =  7;
const sz EL_TYPE_Br   =  8;
const sz EL_TYPE_I    =  9;
const sz EL_TYPE_Met  = 10;
const sz EL_TYPE_SIZE = 11;

// AutoDock4
const sz AD_TYPE_C    =  0;
const sz AD_TYPE_A    =  1;
const sz AD_TYPE_N    =  2;
const sz AD_TYPE_O    =  3;
const sz AD_TYPE_P    =  4;
const sz AD_TYPE_S    =  5;
const sz AD_TYPE_H    =  6; // non-polar hydrogen
const sz AD_TYPE_F    =  7;
const sz AD_TYPE_I    =  8;
const sz AD_TYPE_NA   =  9;
const sz AD_TYPE_OA   = 10;
const sz AD_TYPE_SA   = 11;
const sz AD_TYPE_HD   = 12;
const sz AD_TYPE_Mg   = 13;
const sz AD_TYPE_Mn   = 14;
const sz AD_TYPE_Zn   = 15;
const sz AD_TYPE_Ca   = 16;
const sz AD_TYPE_Fe   = 17;
const sz AD_TYPE_Cl   = 18;
const sz AD_TYPE_Br   = 19;
const sz AD_TYPE_SIZE = 20;

// X-Score
const sz XS_TYPE_C_H   =  0;
const sz XS_TYPE_C_P   =  1;
const sz XS_TYPE_N_P   =  2;
const sz XS_TYPE_N_D   =  3;
const sz XS_TYPE_N_A   =  4;
const sz XS_TYPE_N_DA  =  5;
const sz XS_TYPE_O_P   =  6;
const sz XS_TYPE_O_D   =  7;
const sz XS_TYPE_O_A   =  8;
const sz XS_TYPE_O_DA  =  9;
const sz XS_TYPE_S_P   = 10;
const sz XS_TYPE_P_P   = 11;
const sz XS_TYPE_F_H   = 12;
const sz XS_TYPE_Cl_H  = 13;
const sz XS_TYPE_Br_H  = 14;
const sz XS_TYPE_I_H   = 15;
const sz XS_TYPE_Met_D = 16;
const sz XS_TYPE_SIZE  = 17;

// DrugScore-CSD
const sz SY_TYPE_C_3   =  0;
const sz SY_TYPE_C_2   =  1;
const sz SY_TYPE_C_ar  =  2;
const sz SY_TYPE_C_cat =  3;
const sz SY_TYPE_N_3   =  4;
const sz SY_TYPE_N_ar  =  5;
const sz SY_TYPE_N_am  =  6;
const sz SY_TYPE_N_pl3 =  7;
const sz SY_TYPE_O_3   =  8;
const sz SY_TYPE_O_2   =  9;
const sz SY_TYPE_O_co2 = 10;
const sz SY_TYPE_S     = 11;
const sz SY_TYPE_P     = 12;
const sz SY_TYPE_F     = 13;
const sz SY_TYPE_Cl    = 14;
const sz SY_TYPE_Br    = 15;
const sz SY_TYPE_I     = 16;
const sz SY_TYPE_Met   = 17;
const sz SY_TYPE_SIZE  = 18;

struct atom_kind {
    std::string name;
    fl radius;
	fl depth;
	fl solvation;
	fl volume;
	fl covalent_radius; // from http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
};

// generated from edited AD4_parameters.data using a script, 
// then covalent radius added from en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
const atom_kind atom_kind_data[] = { // name, radius, depth, solvation parameter, volume, covalent radius
	{ "C",    2.00000,    0.15000,   -0.00143,   33.51030,   0.77}, //  0
	{ "A",    2.00000,    0.15000,   -0.00052,   33.51030,   0.77}, //  1
	{ "N",    1.75000,    0.16000,   -0.00162,   22.44930,   0.75}, //  2
	{ "O",    1.60000,    0.20000,   -0.00251,   17.15730,   0.73}, //  3
	{ "P",    2.10000,    0.20000,   -0.00110,   38.79240,   1.06}, //  4
	{ "S",    2.00000,    0.20000,   -0.00214,   33.51030,   1.02}, //  5
	{ "H",    1.00000,    0.02000,    0.00051,    0.00000,   0.37}, //  6
	{ "F",    1.54500,    0.08000,   -0.00110,   15.44800,   0.71}, //  7
	{ "I",    2.36000,    0.55000,   -0.00110,   55.05850,   1.33}, //  8
	{"NA",    1.75000,    0.16000,   -0.00162,   22.44930,   0.75}, //  9
	{"OA",    1.60000,    0.20000,   -0.00251,   17.15730,   0.73}, // 10
	{"SA",    2.00000,    0.20000,   -0.00214,   33.51030,   1.02}, // 11
	{"HD",    1.00000,    0.02000,    0.00051,    0.00000,   0.37}, // 12
	{"Mg",    0.65000,    0.87500,   -0.00110,    1.56000,   1.30}, // 13
	{"Mn",    0.65000,    0.87500,   -0.00110,    2.14000,   1.39}, // 14
	{"Zn",    0.74000,    0.55000,   -0.00110,    1.70000,   1.31}, // 15
	{"Ca",    0.99000,    0.55000,   -0.00110,    2.77000,   1.74}, // 16
	{"Fe",    0.65000,    0.01000,   -0.00110,    1.84000,   1.25}, // 17
	{"Cl",    2.04500,    0.27600,   -0.00110,   35.82350,   0.99}, // 18
	{"Br",    2.16500,    0.38900,   -0.00110,   42.56610,   1.14}  // 19
};

const fl metal_solvation_parameter = -0.00110;

const fl metal_covalent_radius = 1.75; // for metals not on the list // FIXME this info should be moved to non_ad_metals

const sz atom_kinds_size =  sizeof(atom_kind_data) / sizeof(const atom_kind);

struct atom_equivalence {
	std::string name;
	std::string to;
};

const atom_equivalence atom_equivalence_data[] = {
	{"Se",  "S"}
};

const sz atom_equivalences_size = sizeof(atom_equivalence_data) / sizeof(const atom_equivalence);

struct acceptor_kind {
	sz ad_type;
	fl radius;
	fl depth;
};

const acceptor_kind acceptor_kind_data[] = { // ad_type, optimal length, depth
	{AD_TYPE_NA, 1.9, 5.0},
	{AD_TYPE_OA, 1.9, 5.0},
	{AD_TYPE_SA, 2.5, 1.0}
};

const sz acceptor_kinds_size = sizeof(acceptor_kind_data) / sizeof(acceptor_kind);

inline bool ad_is_hydrogen(sz ad) {
	return ad == AD_TYPE_H || ad == AD_TYPE_HD;
}

inline bool ad_is_heteroatom(sz ad) { // returns false for ad >= AD_TYPE_SIZE
	return ad != AD_TYPE_A && ad != AD_TYPE_C  && 
		   ad != AD_TYPE_H && ad != AD_TYPE_HD && 
		   ad < AD_TYPE_SIZE;
}

inline sz ad_type_to_el_type(sz t) {
	switch(t) {
		case AD_TYPE_C    : return EL_TYPE_C;
		case AD_TYPE_A    : return EL_TYPE_C;
		case AD_TYPE_N    : return EL_TYPE_N;
		case AD_TYPE_O    : return EL_TYPE_O;
		case AD_TYPE_P    : return EL_TYPE_P;
		case AD_TYPE_S    : return EL_TYPE_S;
		case AD_TYPE_H    : return EL_TYPE_H;
		case AD_TYPE_F    : return EL_TYPE_F;
		case AD_TYPE_I    : return EL_TYPE_I;
		case AD_TYPE_NA   : return EL_TYPE_N;
		case AD_TYPE_OA   : return EL_TYPE_O;
		case AD_TYPE_SA   : return EL_TYPE_S;
		case AD_TYPE_HD   : return EL_TYPE_H;
		case AD_TYPE_Mg   : return EL_TYPE_Met;
		case AD_TYPE_Mn   : return EL_TYPE_Met;
		case AD_TYPE_Zn   : return EL_TYPE_Met;
		case AD_TYPE_Ca   : return EL_TYPE_Met;
		case AD_TYPE_Fe   : return EL_TYPE_Met;
		case AD_TYPE_Cl   : return EL_TYPE_Cl;
		case AD_TYPE_Br   : return EL_TYPE_Br;
		case AD_TYPE_SIZE : return EL_TYPE_SIZE;
		default: VINA_CHECK(false);
	}
	return EL_TYPE_SIZE; // to placate the compiler in case of warnings - it should never get here though
}

const fl xs_vdw_radii[] = {
	1.9, // C_H
	1.9, // C_P
	1.8, // N_P
	1.8, // N_D
	1.8, // N_A
	1.8, // N_DA
	1.7, // O_P
	1.7, // O_D
	1.7, // O_A
	1.7, // O_DA
	2.0, // S_P
	2.1, // P_P
	1.5, // F_H
	1.8, // Cl_H
	2.0, // Br_H
	2.2, // I_H
	1.2  // Met_D
};

inline fl xs_radius(sz t) {
	const sz n = sizeof(xs_vdw_radii) / sizeof(const fl);
	assert(n == XS_TYPE_SIZE);
	assert(t < n);
	return xs_vdw_radii[t];
}

const std::string non_ad_metal_names[] = { // expand as necessary
	"Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};

inline bool is_non_ad_metal_name(const std::string& name) {
	const sz s = sizeof(non_ad_metal_names) / sizeof(const std::string);
	VINA_FOR(i, s)
		if(non_ad_metal_names[i] == name)
			return true;
	return false;
}

inline bool xs_is_hydrophobic(sz xs) {
	return xs == XS_TYPE_C_H || 
		   xs == XS_TYPE_F_H ||
		   xs == XS_TYPE_Cl_H ||
		   xs == XS_TYPE_Br_H || 
		   xs == XS_TYPE_I_H;
}

inline bool xs_is_acceptor(sz xs) {
	return xs == XS_TYPE_N_A ||
		   xs == XS_TYPE_N_DA ||
		   xs == XS_TYPE_O_A ||
		   xs == XS_TYPE_O_DA;
}

inline bool xs_is_donor(sz xs) {
	return xs == XS_TYPE_N_D ||
		   xs == XS_TYPE_N_DA ||
		   xs == XS_TYPE_O_D ||
		   xs == XS_TYPE_O_DA ||
		   xs == XS_TYPE_Met_D;
}

inline bool xs_donor_acceptor(sz t1, sz t2) {
	return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(sz t1, sz t2) {
	return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

inline const atom_kind& ad_type_property(sz i) {
	assert(AD_TYPE_SIZE == atom_kinds_size);
    assert(i < atom_kinds_size);
    return atom_kind_data[i];
}

inline sz string_to_ad_type(const std::string& name) { // returns AD_TYPE_SIZE if not found (no exceptions thrown, because metals unknown to AD4 are not exceptional)
    VINA_FOR(i, atom_kinds_size)
		if(atom_kind_data[i].name == name)
			return i;
	VINA_FOR(i, atom_equivalences_size)
		if(atom_equivalence_data[i].name == name)
			return string_to_ad_type(atom_equivalence_data[i].to);
    return AD_TYPE_SIZE;
}

inline fl max_covalent_radius() {
	fl tmp = 0;
	VINA_FOR(i, atom_kinds_size)
		if(atom_kind_data[i].covalent_radius > tmp)
			tmp = atom_kind_data[i].covalent_radius;
	return tmp;
}

#endif
