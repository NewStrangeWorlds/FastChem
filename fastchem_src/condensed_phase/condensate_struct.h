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

#ifndef _condensate_struct_h
#define _condensate_struct_h

#include <vector>
#include <string>

#include "../species_struct.h"
#include "../fastchem_constants.h"


namespace fastchem {


enum phase_states { liquid, solid };

//Parent class of all condensates
template <class double_type>
struct Condensate
{
  std::string symbol;
  std::string name;

  std::vector<unsigned int> element_indices;
  std::vector<int> stoichiometric_vector;

  phase_states phase;
  double_type density = 0.0;
  std::vector<double_type> fit_coeff;
  std::vector<double> fit_temp_limits;
  std::vector<double> phase_temp_limits;

  double_type mass_action_constant = 0.0;

  double_type log_activity = 0;
  double_type activity_correction = 0;
  double_type tau = 0;
  bool linear_system_remove = false;

  double_type number_density = 0;         //fictitious number density (not a real one because a condensate "molecule" normally doesn't exist)
  double_type max_number_density = 0;
  
  unsigned int reference_element = FASTCHEM_UNKNOWN_SPECIES;     //the element, the degree of condensation is defined for
  double_type degree_of_condensation = 0;

  void calcMassActionConstant(const double temperature);
  void calcActivity(
    const double temperature, const std::vector<Element<double_type>>& elements);
  void findReferenceElement(
    const std::vector<Element<double_type>>& elements);
  void degreeOfCondensation(
    const std::vector<Element<double_type>>& elements, const double_type total_element_density);
  void maxDensity(
    const std::vector< Element<double_type> >& elements, double_type total_number_density);
};


}

#endif
