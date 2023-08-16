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


#include <cmath>
#include <iostream>
#include <algorithm>

#include "../species_struct.h"

#include "../fastchem_constants.h"


namespace fastchem {


template <class double_type>
void Condensate<double_type>::calcMassActionConstant(const double temperature)
{
  size_t set_index = 0;

  for (size_t i=0; i<fit_coeff_limits.size(); ++i)
    if (temperature <= fit_coeff_limits[i])
    {
      set_index = i;
      break;
    }

  if (fit_coeff_limits.back() < temperature) 
    set_index = fit_coeff_limits.size()-1;

  double_type log_K = fit_coeff[set_index][0]/temperature
                    + fit_coeff[set_index][1]*std::log(temperature)
                    + fit_coeff[set_index][2]
                    + fit_coeff[set_index][3]*temperature
                    + fit_coeff[set_index][4]*temperature * temperature;

  double_type sigma = 0;

  for (auto & i : stoichiometric_vector)
    sigma += i;

  const double_type pressure_scaling = 1.0e6 / (CONST_K * temperature);
  mass_action_constant = log_K - sigma * std::log(pressure_scaling);
}



template <class double_type>
void Condensate<double_type>::calcActivity(
  const double temperature,
  const std::vector<Element<double_type>>& elements,
  const bool use_data_validity_limits)
{
  if (use_data_validity_limits && temperature > fit_coeff_limits.back())
  {
    log_activity = -10.0;
    return;
  }

  log_activity = mass_action_constant;

  for (auto & i : element_indices)
    log_activity += std::log(elements[i].number_density) * stoichiometric_vector[elements[i].index];

  if (log_activity < -10.0) log_activity = -10.0;
}



template <class double_type>
double_type Condensate<double_type>::calcActivity(
  const double temperature, 
  const std::vector<Element<double_type>>& elements,
  const std::vector<double_type> elem_number_densities,
  const bool use_data_validity_limits)
{
  if (use_data_validity_limits && temperature > fit_coeff_limits.back())
    return -10.0;

  double_type log_activity = mass_action_constant;

  for (auto & i : element_indices)
    log_activity += std::log(elem_number_densities[elements[i].index]) * stoichiometric_vector[elements[i].index];

  if (log_activity < -10.0) log_activity = -10.0;

  return log_activity;
}



//the reference element is the one with the smallest abundance in the condensate
template <class double_type>
void Condensate<double_type>::findReferenceElement(
  const std::vector<Element<double_type>>& elements)
{
  reference_element = element_indices[0];
  
  double_type smallest_abundance = 
    elements[element_indices[0]].abundance/stoichiometric_vector[element_indices[0]];

  for (auto & i : element_indices)
    if (elements[i].abundance/stoichiometric_vector[i] < smallest_abundance)
    {
      reference_element = elements[i].index;
      smallest_abundance = elements[i].abundance/stoichiometric_vector[i];
    }
}



//degree of condensation for the condensate
template <class double_type>
void Condensate<double_type>::degreeOfCondensation(
  const std::vector<Element<double_type>>& elements, const double_type total_element_density)
{
  if (reference_element == FASTCHEM_UNKNOWN_SPECIES)
    findReferenceElement(elements);

  degree_of_condensation = stoichiometric_vector[reference_element] 
    * this->number_density /(total_element_density * elements[reference_element].epsilon);
}


//maximum condensate density
//see Eq. 13 in Paper III
template <class double_type>
void Condensate<double_type>::maxDensity(
  const std::vector< Element<double_type> >& elements, double_type total_number_density)
{
  max_number_density = elements[element_indices[0]].epsilon 
    * total_number_density/stoichiometric_vector[element_indices[0]];

  for (auto & i : element_indices)
  { 
    const double_type max_density = elements[i].epsilon * total_number_density/stoichiometric_vector[i];

    if (max_density < max_number_density) max_number_density = max_density;
  }
}


template struct Condensate<double>;
template struct Condensate<long double>;

}
