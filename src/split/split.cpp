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

#include <iostream>
#include <string>
#include <exception>
#include <vector> // ligand paths
#include <list>
#include <fstream> //  getline?
#include <cmath> // for ceila
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename

#include "file.h"
#include "parse_error.h"

path make_path(const std::string& str) {
    boost::filesystem::path p(str);
	return p;
}

std::string default_prefix(const std::string& input_name, const std::string& add) {
	std::string tmp = input_name;
	if(tmp.size() >= 6 && tmp.substr(tmp.size()-6, 6) == ".pdbqt")
		tmp.resize(tmp.size() - 6); // FIXME?
	return tmp + add;
}

struct usage_error : public std::runtime_error {
	usage_error(const std::string& message) : std::runtime_error(message) {}
};

typedef std::list<std::string> strl;

struct model {
	strl ligand;
	strl flex;
};

typedef std::list<model> models;

models parse_multimodel_pdbqt(const std::string& input) {
	const path p = make_path(input);
	ifile in(p);
	models tmp;
	unsigned count = 0;
	std::string str;
	bool parsing_model = false;
	bool parsing_ligand = true;

	while(std::getline(in, str)) {
		++count;
		if(starts_with(str, "MODEL")) {
			if(parsing_model == true || parsing_ligand == false) 
				throw pdbqt_parse_error("Misplaced MODEL tag at line " + std::to_string(count) + ".");
			tmp.push_back(model());
			parsing_model = true;
		}
		else if(starts_with(str, "ENDMDL")) {
			if(parsing_model == false || parsing_ligand == false)
				throw pdbqt_parse_error("Misplaced ENDMDL tag at line " + std::to_string(count) + ".");
			parsing_model = false;
		}
		else if(starts_with(str, "BEGIN_RES")) {
			if(parsing_model == false || parsing_ligand == false)
				throw pdbqt_parse_error("Misplaced BEGIN_RES tag at line " + std::to_string(count) + ".");
			parsing_ligand = false;
			tmp.back().flex.push_back(str);
		}
		else if(starts_with(str, "END_RES")) {
			if(parsing_model == false || parsing_ligand == true)
				throw pdbqt_parse_error("Misplaced END_RES tag at line " + std::to_string(count) + ".");
			parsing_ligand = true;
			tmp.back().flex.push_back(str);
		}
		else {
			if(parsing_model == false)
				throw pdbqt_parse_error("Input occurs outside MODEL at line " + std::to_string(count) + ".");
			if(parsing_ligand) {
				tmp.back().ligand.push_back(str);
			}
			else {
				tmp.back().flex.push_back(str);
			}
		}
	}
	if(parsing_model == true)
		throw pdbqt_parse_error("Missing ENDMDL tag at line " + std::to_string(count + 1) + ".");
	return tmp;
}

void write_pdbqt(const strl& lines, const std::string& name) {
	if(lines.size() > 0) {
		ofile out(make_path(name));
		for(strl::const_iterator it = lines.begin(); it != lines.end(); ++it)
			out << *it << '\n';
	}
}

void write_multimodel_pdbqt(const models& m, const std::string& ligand_prefix, const std::string& flex_prefix) {
	sz how_many = m.size();
	std::streamsize w = static_cast<std::streamsize>(to_string(how_many).size());
	sz counter = 0;
	for(models::const_iterator it = m.begin(); it != m.end(); ++it) {
		++counter;
		const std::string add = to_string(counter, w, '0') + ".pdbqt";
		write_pdbqt(it->ligand, ligand_prefix + add);
		write_pdbqt(it->flex,   flex_prefix   + add);
	}
}

int main(int argc, char* argv[]) {
	using namespace boost::program_options;
	const bool advanced = false;
	const std::string git_version = VERSION;
	const std::string version_string = "AutoDock Vina PDBQT Split " + git_version;

	const std::string error_message = "\n\n\
Please contact the author, Dr. Oleg Trott <ot14@columbia.edu>, so\n\
that this problem can be resolved. The reproducibility of the\n\
error may be vital, so please remember to include the following in\n\
your problem report:\n\
* the EXACT error message,\n\
* your version of the program,\n\
* the computer system you are running it on,\n\
* command line and configuration file options,\n\
* input (if possible),\n\
\n\
Thank you!\n";

	try {
		std::string input_name, ligand_prefix = "ligand_", flex_prefix = "flex_";
		bool help = false, version = false;
		options_description inputs("Input");
		inputs.add_options()
			("input", value<std::string>(&input_name), "input to split (PDBQT)")
		;
		options_description outputs("Output (optional) - defaults are chosen based on the input file name");
		outputs.add_options()
			("ligand", value<std::string>(&ligand_prefix), "prefix for ligands")
			("flex", value<std::string>(&flex_prefix), "prefix for side chains")
		;
		options_description info("Information (optional)");
		info.add_options()
			("help", bool_switch(&help), "print this message")
			("version", bool_switch(&version), "print program version")
		;
		options_description desc;
		desc.add(inputs).add(outputs).add(info);

		std::cout << version_string << '\n';
		positional_options_description positional; // remains empty
		variables_map vm;
		try {
			store(command_line_parser(argc, argv)
				.options(desc)
				.style(command_line_style::default_style ^ command_line_style::allow_guessing)
				.positional(positional)
				.run(), 
				vm);
			notify(vm); 
		}
		catch(boost::program_options::error& e) {
			std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc << '\n';
			return 1;
		}
		if(help) {
			std::cout << desc << '\n';
			return 0;
		}
		if(version) {
			return 0;
		}

		if(vm.count("input") <= 0) {
			std::cerr << "Missing input.\n" << "\nCorrect usage:\n" << desc << '\n';
			return 1;
		}
		if(vm.count("ligand") <= 0) {
			ligand_prefix = default_prefix(input_name, "_ligand_");
			std::cout << "Prefix for ligands will be " << ligand_prefix << '\n';
		}
		if(vm.count("flex") <= 0) {
			flex_prefix = default_prefix(input_name, "_flex_");
			std::cout << "Prefix for flexible side chains will be " << flex_prefix << '\n';
		}
		const models tmp = parse_multimodel_pdbqt(input_name);
		write_multimodel_pdbqt(tmp, ligand_prefix, flex_prefix);
	}
	catch(file_error& e) {
		std::cerr << "\n\nError: could not open \"" << e.name.string() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
		return 1;
	}
	catch(boost::filesystem::filesystem_error& e) {
		std::cerr << "\n\nFile system error: " << e.what() << '\n';
		return 1;
	}
	catch(usage_error& e) {
		std::cerr << "\n\nUsage error: " << e.what() << ".\n";
		return 1;
	}
	catch(pdbqt_parse_error& e) {
		std::cerr << e.what();
        return 1;
	}
	catch(std::bad_alloc&) {
		std::cerr << "\n\nError: insufficient memory!\n";
		return 1;
	}

	// Errors that shouldn't happen:

	catch(std::exception& e) { 
		std::cerr << "\n\nAn error occurred: " << e.what() << ". " << error_message;
		return 1; 
	}
	catch(internal_error& e) {
		std::cerr << "\n\nAn internal error occurred in " << e.file << "(" << e.line << "). " << error_message;
		return 1;
	}
	catch(...) {
		std::cerr << "\n\nAn unknown error occurred. " << error_message;
		return 1;
	}
}
