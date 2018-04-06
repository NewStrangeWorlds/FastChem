/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
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


#include "fastchem.h"

#include <string>
#include <vector>



namespace fastchem {


//Constructor for the FastChem class
//Requires: parameter file path and the initial verbose level that is used while reading the input files
template <class double_type>
FastChem<double_type>::FastChem(const std::string& model_parameter_file, const unsigned int verbose_level_init)
{
  verbose_level = verbose_level_init;

  bool parameter_file_loaded = false;

  if (model_parameter_file != "")
    parameter_file_loaded = readParameterFile(model_parameter_file);


  if (!parameter_file_loaded)
  {
    std::cout << "Error reading parameters\n";
    is_initialized = false;
  }


  if (parameter_file_loaded) init();
}



//Copy constructor
//Could be made more pretty, but this one does the job as well...
template <class double_type>
FastChem<double_type>::FastChem(const FastChem &obj)
{
  nb_chemical_elements = obj.nb_chemical_elements;
  nb_species = obj.nb_species;
  nb_molecules = obj.nb_molecules;
  nb_elements = obj.nb_elements;

  nb_max_fastchem_iter = obj.nb_max_fastchem_iter;
  nb_max_pressure_iter = obj.nb_max_pressure_iter;
  nb_max_bisection_iter = obj.nb_max_bisection_iter;
  nb_max_neldermead_iter = obj.nb_max_neldermead_iter;
  nb_max_newton_iter = obj.nb_max_newton_iter;

  element_density_minlimit = obj.element_density_minlimit;
  molecule_density_minlimit = obj.molecule_density_minlimit;

  accuracy = obj.accuracy;
  accuracy_delta = obj.accuracy_delta;
  newton_err = obj.newton_err;

  verbose_level = obj.verbose_level;
  use_scaling_factor = obj.use_scaling_factor;
  is_initialized = obj.is_initialized;


  chemical_element_file = obj.chemical_element_file;
  species_data_file = obj.species_data_file;
  element_abundances_file = obj.element_abundances_file;


  chemical_elements = obj.chemical_elements;
  elements = obj.elements;
  molecules = obj.molecules;

  e_ = obj.e_;

  element_calculation_order = obj.element_calculation_order;

  for (size_t i=0; i<nb_elements; ++i)
    species.push_back(&elements[i]);

  for (size_t i=0; i<nb_molecules; ++i)
    species.push_back(&molecules[i]);
}



template class FastChem<double>;
template class FastChem<long double>;


}
