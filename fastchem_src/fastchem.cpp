/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2022 Daniel Kitzmann, Joachim Stock
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


#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "fastchem.h"


namespace fastchem {


//Constructor for the FastChem class
//Requires: parameter file path and the initial verbose level that is used while reading the input files
template <class double_type>
FastChem<double_type>::FastChem(
  const std::string& model_parameter_file, const unsigned int verbose_level_start)
    : options(model_parameter_file, verbose_level_start)
    , element_data(options.element_abundances_file, options.chemical_element_file)
    , gas_phase(options, element_data)
    , condensed_phase(options, element_data)
{
  if (!options.parameter_file_loaded)
  {
    std::cout << "Error reading parameters\n";
    is_initialised = false;

    return;
  }
  
  if (element_data.is_initialised == true 
      && gas_phase.is_initialised == true 
      && (condensed_phase.is_initialised == true || options.condensates_data_file == "none"))
    is_initialised = true;
  else
  {
    std::cout << "Error initialising FastChem!\n\n";
    is_initialised = false;

    return;
  }

  if (options.verbose_level >= 1)
    std::cout << "\nFastChem initialisation summary:\n"
              << "  number of species: " << gas_phase.nb_species + condensed_phase.nb_condensates
              << "    elements: " << element_data.nb_elements
              << "    molecules: " << gas_phase.nb_molecules
              << "    condensates: " << condensed_phase.nb_condensates
              << "\n\n";

  init();
}



//Constructor for the FastChem class
//Requires: file paths for element abundance and species data files 
//          and the initial verbose level that is used while reading the input files
template <class double_type>
FastChem<double_type>::FastChem(
  const std::string& element_abundances_file,
  const std::string& species_data_file,
  const std::string& cond_species_data_file,
  const unsigned int verbose_level_start) 
    : options(element_abundances_file, species_data_file, cond_species_data_file, verbose_level_start)
    , element_data(element_abundances_file, options.chemical_element_file)
    , gas_phase(options, element_data)
    , condensed_phase(options, element_data)
{
  if (element_data.is_initialised == true 
      && gas_phase.is_initialised == true 
      && (condensed_phase.is_initialised == true || cond_species_data_file == "none"))
    is_initialised = true;
  else
  {
    std::cout << "Error initialising FastChem!\n\n";
    is_initialised = false;

    return;
  }

  if (options.verbose_level >= 1)
    std::cout << "\nFastChem initialisation summary:\n"
              << "  number of species: " << gas_phase.nb_species + condensed_phase.nb_condensates
              << "    elements: " << element_data.nb_elements
              << "    molecules: " << gas_phase.nb_molecules
              << "    condensates: " << condensed_phase.nb_condensates
              << "\n\n";

  init();
}



//Constructor for the FastChem class
//This version only initiales the gas phase
//Requires: file paths for element abundance and species data files 
//          and the initial verbose level that is used while reading the input files
template <class double_type>
FastChem<double_type>::FastChem(
  const std::string& element_abundances_file,
  const std::string& species_data_file,
  const unsigned int verbose_level_start) 
    : options(element_abundances_file, species_data_file, std::string("none"), verbose_level_start)
    , element_data(element_abundances_file, options.chemical_element_file)
    , gas_phase(options, element_data)
    , condensed_phase(options, element_data)
{
  if (element_data.is_initialised == true && gas_phase.is_initialised == true)
    is_initialised = true;
  else
  {
    std::cout << "Error initialising FastChem!\n\n";
    is_initialised = false;

    return;
  }

  if (options.verbose_level >= 1)
    std::cout << "\nFastChem initialisation summary:\n"
              << "  number of species: " << gas_phase.nb_species + condensed_phase.nb_condensates
              << "    elements: " << element_data.nb_elements
              << "    molecules: " << gas_phase.nb_molecules
              << "    condensates: " << condensed_phase.nb_condensates
              << "\n\n";

  init();
}



template <class double_type>
void FastChem<double_type>::init()
{

}



//the copy constructor
template <class double_type>
FastChem<double_type>::FastChem(const FastChem &obj)
  : options(obj.options)
  , element_data(obj.element_data)
  , gas_phase(obj.gas_phase, options, element_data)
  , condensed_phase(obj.condensed_phase, options, element_data)
{

}


template class FastChem<double>;
template class FastChem<long double>;
}
