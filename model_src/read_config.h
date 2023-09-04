/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2021 Daniel Kitzmann, Joachim Stock
*
* FastChem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* FastChem is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* FastChem directory under <license.md>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _read_config_h
#define _read_config_h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>


struct Config{
  std::string atmosphere_file = "";
  std::string chem_output_file = "";
  std::string cond_output_file = "";
  std::string monitor_output_file = "";
  
  bool calc_condensation = false;
  bool rainout_condensation = false;

  unsigned int verbose_level = 1; 
  bool output_mixing_ratios = true;

  std::string element_abundance_file = "";
  std::string species_data_file = "";
  std::string cond_species_data_file = "";
  double chemistry_accuracy = 0;
  double element_conservation_accuracy = 0;
  double newton_error = 0;

  unsigned int nb_chemistry_iterations = 0;
  unsigned int nb_newton_iterations = 0;
  unsigned int nb_nelder_mead_iterations = 0;
  unsigned int nb_bisection_iterations = 0;
};



bool readConfigFile(std::string &file_path, Config &config)
{
  std::fstream file;

  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to read config file: " << file_path << "\n";
    return false;
  }

  std::string file_name, line;
  int nb_value;

  //FastChem parameter file
  std::getline(file, line); std::getline(file, line);
  std::istringstream input(line);

  input >> config.atmosphere_file;

  if (config.atmosphere_file == "")
  {
    std::cout << "Unable to read p-T file location from: " << file_path << "\n";

    return false;
  }

  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  std::string calc_type;

  input.str(line); input.clear();
  input >> calc_type;

  if (calc_type != "g" && calc_type != "ce" && calc_type !="cr")
  {
    std::cout << "Unknown calculation type \"" << calc_type << "\" found in: " << file_path << "\n";

    return false;
  }

  if (calc_type == "ce")
  {
    config.calc_condensation = true;
    config.rainout_condensation = false;
  }

  if (calc_type == "cr")
  {
    config.calc_condensation = true;
    config.rainout_condensation = true;
  }

  //chemistry output file
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> config.chem_output_file;

  if (config.chem_output_file == "")
  {
    std::cout << "Unable to read chemistry output file location from: " << file_path << "\n";

    return false;
  }

  input >> config.cond_output_file;
  
  if (config.calc_condensation && config.cond_output_file == "")
  {
    std::cout << "Unable to read condensate output file location from: " << file_path << "\n";

    return false;
  }


  //monitor output file
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> config.monitor_output_file;

  if (config.monitor_output_file == "")
  {
    std::cout << "Unable to read monitor output file location from: " << file_path.c_str() << "\n";

    return false;
  }


  //monitor output file
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> nb_value;

  if (line == "" || !(nb_value >= 0))
  {
    std::cout << "Unable to read verbose option from: " << file_path.c_str() << "\n";

    return false;
  }

  config.verbose_level = nb_value;

  
  //output option
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> line;

  if (line == "MR")
    config.output_mixing_ratios = true;
  else
    config.output_mixing_ratios = false;

  
  //element abundance file
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> config.element_abundance_file;

  if (config.element_abundance_file == "")
  {
    std::cout << "Unable to read element abundance file from: " << file_path.c_str() << "\n";

    return false;
  }


  //species data file
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> config.species_data_file;

  if (config.species_data_file == "")
  {
    std::cout << "Unable to read species data file from: " << file_path.c_str() << "\n";

    return false;
  }

  input >> config.cond_species_data_file;

  if (config.calc_condensation 
      && (config.cond_species_data_file == "" || config.cond_species_data_file == "none"))
  {
    std::cout << "Unable to read condensate species data file from: " << file_path.c_str() << "\n";

    return false;
  }

  if (config.cond_species_data_file == "") 
    config.cond_species_data_file = "none";


  //chemistry accuracy
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> config.chemistry_accuracy;

  if (line == "" || !(config.chemistry_accuracy > 0.0))
  {
    std::cout << "Unable to chemistry accuracy parameter from: " << file_path.c_str() << "\n";

    return false;
  }

  config.newton_error = config.chemistry_accuracy;


  //element conservation accuracy
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> config.element_conservation_accuracy;

  if (line == "" || !(config.element_conservation_accuracy > 0.0))
  {
    std::cout << "Unable to element conservation accuracy parameter from: " << file_path.c_str() << "\n";

    return false;
  }

  config.newton_error = config.chemistry_accuracy;


  //chemistry iteration number
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> nb_value;

  if (!(nb_value > 0))
  {
    std::cout << "Unable to read max chemistry iteration number from: " << file_path.c_str() << "\n";

    return false;
  }
  
  config.nb_chemistry_iterations = nb_value;


  //solver iteration number
  std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> nb_value;

  if (!(nb_value > 0))
  {
    std::cout << "Unable to read max solver iteration number from: " << file_path.c_str() << "\n";

    return false;
  }
  
  config.nb_newton_iterations = nb_value;
  config.nb_nelder_mead_iterations = nb_value;
  config.nb_bisection_iterations = nb_value;


  if (config.verbose_level >= 2)
  std::cout << "\nRead from config file " << file_path << " : \n"
            << "  Temperature-pressure file: " << config.atmosphere_file << "\n"
            << "  Chemistry output file: " << config.chem_output_file << "\n"
            << "  Monitor file: " << config.monitor_output_file << "\n"
            << "  Verbose level: " << config.verbose_level << "\n"
            << "  Output mixing ratios: " << config.output_mixing_ratios << "\n"
            << "  Element abundance file: " << config.element_abundance_file << "\n"
            << "  Species data file: " << config.species_data_file << "\n"
            << "  Chemistry accuracy: " << config.chemistry_accuracy << "\n"
            << "  Element conservation accuracy: " << config.element_conservation_accuracy << "\n"
            << "  Number of chemistry iterations: " << config.nb_chemistry_iterations << "\n"
            << "  Number of solver iterations: " << config.nb_newton_iterations << "\n"
            << "\n";

  return true;
}


#endif