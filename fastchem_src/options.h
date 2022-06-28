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


#ifndef _options_h
#define _options_h

#include <vector>
#include <iostream>
#include <string>


namespace fastchem {


template <class double_type>
struct FastChemOptions{
  FastChemOptions(
    const std::string& parameter_file,
    const unsigned int verbose_level_start)
      : verbose_level(verbose_level_start) {
          parameter_file_loaded = readParameterFile(parameter_file);
      };
  FastChemOptions(
    const std::string& element_abundances_file_,
    const std::string& species_data_file_,
    const unsigned int verbose_level_start)
      : verbose_level(verbose_level_start) {
          element_abundances_file = element_abundances_file_;
          species_data_file = species_data_file_;
      };

  unsigned int nb_max_fastchem_iter   = 1000;
  unsigned int nb_max_bisection_iter  = 3000;
  unsigned int nb_max_newton_iter     = 3000;
  unsigned int nb_max_neldermead_iter = 3000;

  double_type accuracy = 1e-4;
  double_type newton_err = 1e-4;

  double_type element_density_minlimit = 1e-300; //smallest allowed particle number densities
  double_type molecule_density_minlimit = 1e-300;

  unsigned int verbose_level = 1;

  bool use_scaling_factor = false;
  double additional_scaling_factor = 0.0;

  std::string chemical_element_file = "";
  std::string species_data_file = "";
  std::string element_abundances_file = "";

  bool parameter_file_loaded = false;
  bool readParameterFile(const std::string& model_parameter_file);
};



}
#endif
