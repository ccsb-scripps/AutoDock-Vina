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
#include <boost/program_options.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include "parse_error.h"
#include "tee.h"
#include "vina.h"
#include "utils.h"
#include "scoring_function.h"

struct usage_error : public std::runtime_error {
	usage_error(const std::string& message) : std::runtime_error(message) {}
};

struct options_occurrence {
	bool some;
	bool all;
	options_occurrence() : some(false), all(true) {} // convenience
	options_occurrence& operator+=(const options_occurrence& x) {
		some = some || x.some;
		all  = all  && x.all;
		return *this;
	}
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
	options_occurrence tmp;
	VINA_FOR_IN(i, d.options()) 
		if(vm.count((*d.options()[i]).long_name())) 
			tmp.some = true;
		else 
			tmp.all = false;
	return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm, boost::program_options::options_description& d) {
	VINA_FOR_IN(i, d.options()) {
		const std::string& str = (*d.options()[i]).long_name();
		if(!vm.count(str))
			std::cerr << "Required parameter --" << str << " is missing!\n";
	}
}

int main(int argc, char* argv[]) {
	using namespace boost::program_options;
	const std::string version_string = "AutoDock Vina 1.1.2 (May 11, 2011)";
	const std::string error_message = "\n\n\
Please contact the author, Dr. Oleg Trott <ot14@columbia.edu>, so\n\
that this problem can be resolved. The reproducibility of the\n\
error may be vital, so please remember to include the following in\n\
your problem report:\n\
* the EXACT error message,\n\
* your version of the program,\n\
* the type of computer system you are running it on,\n\
* all command line options,\n\
* configuration file (if used),\n\
* ligand file as PDBQT,\n\
* receptor file as PDBQT,\n\
* flexible side chains file as PDBQT (if used),\n\
* output file as PDBQT (if any),\n\
* input (if possible),\n\
* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";

	const std::string cite_message = "\
#################################################################\n\
# If you used AutoDock Vina in your work, please cite:          #\n\
#                                                               #\n\
# O. Trott, A. J. Olson,                                        #\n\
# AutoDock Vina: improving the speed and accuracy of docking    #\n\
# with a new scoring function, efficient optimization and       #\n\
# multithreading, Journal of Computational Chemistry 31 (2010)  #\n\
# 455-461                                                       #\n\
#                                                               #\n\
# DOI 10.1002/jcc.21334                                         #\n\
#                                                               #\n\
# Please see http://vina.scripps.edu for more information.      #\n\
#################################################################\n";

	try {
		tee log;
		std::string rigid_name, ligand_name, flex_name, config_name, out_name, log_name;
        std::string ad4_maps;
        std::string sf_name = "vina";
		double center_x, center_y, center_z; 
		double size_x, size_y, size_z;
		int cpu = 0, seed, exhaustiveness = 8, verbosity = 2, num_modes = 9;
		double energy_range = 2.0;
		double granularity = 0.375;
		bool no_cache = false;

		// autodock4.2 weights
		double weight_ad4_vdw   = 0.1662;
		double weight_ad4_hb    = 0.1209;
		double weight_ad4_elec  = 0.1406;
		double weight_ad4_dsolv = 0.1322;
		double weight_ad4_rot   = 0.2983;

        // vina weights
		double weight_gauss1      = -0.035579;
		double weight_gauss2      = -0.005156;
		double weight_repulsion   =  0.840245;
		double weight_hydrophobic = -0.035069;
		double weight_hydrogen    = -0.587439;
		double weight_rot         =  0.05846;
		bool score_only = false, local_only = false, randomize_only = false, help = false, help_advanced = false, version = false; // FIXME

		positional_options_description positional; // remains empty

		options_description inputs("Input");
		inputs.add_options()
			("receptor", value<std::string>(&rigid_name), "rigid part of the receptor (PDBQT)")
			("flex", value<std::string>(&flex_name), "flexible side chains, if any (PDBQT)")
			("ligand", value<std::string>(&ligand_name), "ligand (PDBQT)")
			("ad4_maps", value<std::string>(&ad4_maps), "maps for the autodock4.2 (ad4) scoring function")
			//("vina_maps", value<std::string>(&vina_maps), "maps for vina scoring function")
			("scoring_function", value<std::string>(&sf_name)->default_value(sf_name), "vina or ad4")
		;
		//options_description search_area("Search area (required, except with --score_only)");
		options_description search_area("Search space (required)");
		search_area.add_options()
			("center_x", value<double>(&center_x), "X coordinate of the center")
			("center_y", value<double>(&center_y), "Y coordinate of the center")
			("center_z", value<double>(&center_z), "Z coordinate of the center")
			("size_x", value<double>(&size_x), "size in the X dimension (Angstroms)")
			("size_y", value<double>(&size_y), "size in the Y dimension (Angstroms)")
			("size_z", value<double>(&size_z), "size in the Z dimension (Angstroms)")
		;
		//options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
		options_description outputs("Output (optional)");
		outputs.add_options()
			("out", value<std::string>(&out_name), "output models (PDBQT), the default is chosen based on the ligand file name")
			("log", value<std::string>(&log_name), "optionally, write log file")
		;
		options_description advanced("Advanced options (see the manual)");
		advanced.add_options()
			("score_only",     bool_switch(&score_only),     "score only - search space can be omitted")
			("local_only",     bool_switch(&local_only),     "do local search only")
			("randomize_only", bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")

			("weight_gauss1", value<double>(&weight_gauss1)->default_value(weight_gauss1),                "gauss_1 weight")
			("weight_gauss2", value<double>(&weight_gauss2)->default_value(weight_gauss2),                "gauss_2 weight")
			("weight_repulsion", value<double>(&weight_repulsion)->default_value(weight_repulsion),       "repulsion weight")
			("weight_hydrophobic", value<double>(&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
			("weight_hydrogen", value<double>(&weight_hydrogen)->default_value(weight_hydrogen),          "Hydrogen bond weight")
			("weight_rot", value<double>(&weight_rot)->default_value(weight_rot),                         "N_rot weight")

			("weight_ad4_vdw",   value<double>(&weight_ad4_vdw)  ->default_value(weight_ad4_vdw),   "ad4_vdw weight")
			("weight_ad4_hb",    value<double>(&weight_ad4_hb)   ->default_value(weight_ad4_hb),    "ad4_hb weight")
			("weight_ad4_elec",  value<double>(&weight_ad4_elec) ->default_value(weight_ad4_elec),  "ad4_elec weight")
			("weight_ad4_dsolv", value<double>(&weight_ad4_dsolv)->default_value(weight_ad4_dsolv), "ad4_dsolv weight")
			("weight_ad4_rot",   value<double>(&weight_ad4_rot)  ->default_value(weight_ad4_rot),   "ad4_rot weight")
		;
		options_description misc("Misc (optional)");
		misc.add_options()
			("cpu", value<int>(&cpu), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
			("seed", value<int>(&seed), "explicit random seed")
			("exhaustiveness", value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
			("num_modes", value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
			("energy_range", value<double>(&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
		;
		options_description config("Configuration file (optional)");
		config.add_options()
			("config", value<std::string>(&config_name), "the above options can be put here")
		;
		options_description info("Information (optional)");
		info.add_options()
			("help",          bool_switch(&help), "display usage summary")
			("help_advanced", bool_switch(&help_advanced), "display usage summary with advanced options")
			("version",       bool_switch(&version), "display program version")
		;
		options_description desc, desc_config, desc_simple;
		desc       .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
		desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
		desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);

		variables_map vm;
		try {
			//store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
			store(command_line_parser(argc, argv)
				.options(desc)
				.style(command_line_style::default_style ^ command_line_style::allow_guessing)
				.positional(positional)
				.run(), 
				vm);
			notify(vm); 
		}
		catch(boost::program_options::error& e) {
			std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
			return 1;
		}
		if(vm.count("config")) {
			try {
				path name = make_path(config_name);
				ifile config_stream(name);
				store(parse_config_file(config_stream, desc_config), vm);
				notify(vm);
			}
			catch(boost::program_options::error& e) {
				std::cerr << "Configuration file parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
				return 1;
			}
		}
		if(help) {
			std::cout << desc_simple << '\n';
			return 0;
		}
		if(help_advanced) {
			std::cout << desc << '\n';
			return 0;
		}
		if(version) {
			std::cout << version_string << '\n';
			return 0;
		}

		bool output_produced   = !score_only; 
		bool receptor_needed   = !(sf_name.compare("ad4") == 0);

		if(receptor_needed) {
			if((vm.count("receptor") <= 0)){// & (vm.count("vina_maps") <= 0)){
				std::cerr << desc_simple << "\n\nMissing either receptor or vina_maps.\n";
				return 1;
			}
		}

        if(sf_name.compare("ad4") == 0) {
            //if((vm.count("receptor") > 0) | (vm.count("vina_maps") > 0) | (vm.count("ad4_maps") <= 0)) {
            //if(                               (vm.count("vina_maps") > 0) | (vm.count("ad4_maps") <= 0)) {
            if((vm.count("ad4_maps") <= 0)) {
				std::cerr << "When using AutoDock4.2 scoring function (--scoring_function ad4):\n";
				std::cerr << "  --ad4_maps  is required\n";
				//std::cerr << "  --vina_maps not accepted\n";
				//std::cerr << "  --receptor  not accepted\n";
				return 1;
            }
        } else if (sf_name.compare("vina") == 0) {
            if(vm.count("ad4_maps") > 0) {
				std::cerr << desc_simple << "\n\nNo ad4_maps with vina scoring function\n";
				return 1;
            }
        }
        //if((vm.count("vina_maps") > 0) & (vm.count("receptor") > 0)){
		//	std::cerr << desc_simple << "\n\nMust provide either --receptor or --vina_maps, not both.\n";
		//	return 1;
		//}
		if(vm.count("ligand") <= 0) {
			std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n" << desc_simple << '\n';
			return 1;
		}
		if(vm.count("flex") && !vm.count("receptor"))
			throw usage_error("Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

		if(vm.count("log") > 0)
			log.init(log_name);

		log << cite_message << '\n';

		if(output_produced) { // FIXME
			if(!vm.count("out")) {
				out_name = default_output(ligand_name);
				log << "Output will be " << out_name << '\n';
			}
		}


		doing(verbosity, "Reading input", log);
		done(verbosity, log);

		log << "Exhaustiveness: " << exhaustiveness << "\n";
		log << "CPU: " << cpu << "\n";
		log << "Seed: " << seed << "\n";
		log << "Verbosity: " << verbosity << "\n";
		log << "Rigid receptor: " << rigid_name << "\n";
		log << "Flex receptor: " << flex_name << "\n";
		log << "Ligand: " << ligand_name << "\n";
		log << "Size: X " << size_x << " Y " << size_y << " Z " << size_z << "\n";
		log << "Center: X " << center_x << " Y " << center_y << " Z " << center_z << "\n";
		log << "Grid space: " << granularity << "\n";
		log << "Scoring function : " << sf_name << "\n";
		log.endl();

		Vina v(sf_name, exhaustiveness, cpu, seed, no_cache, verbosity);

        // Can be ignored for AD4
        v.set_receptor(rigid_name, flex_name);
        v.set_ligand(ligand_name);

        // Technically we don't have to initialize weights, 
        // because they are initialized during the Vina object creation with the default weights
        // but we still do it in case the user decided to change them
		if (sf_name.compare("vina") == 0) {
            v.set_vina_weights(weight_gauss1, weight_gauss2, weight_repulsion,
                               weight_hydrophobic, weight_hydrogen, weight_rot);
			v.compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z, granularity);
		} else {
            v.set_ad4_weights(weight_ad4_vdw, weight_ad4_hb, weight_ad4_elec,
                              weight_ad4_dsolv, weight_ad4_rot);
            v.load_ad4_maps(ad4_maps);
        }

		if (score_only) {
			v.score();
            // TODO write contributions
		} else if (local_only) {
			v.optimize();
			v.write_pose(out_name);
		} else {
			v.global_search();
			v.write_results(out_name, num_modes, energy_range);
		}

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
	catch(parse_error& e) {
		std::cerr << "\n\nParse error on line " << e.line << " in file \"" << e.file.string() << "\": " << e.reason << '\n';
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
