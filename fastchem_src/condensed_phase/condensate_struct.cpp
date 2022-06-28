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

#include "condensate_struct.h"

#include "../fastchem_constants.h"
#include "../species_struct.h"


namespace fastchem {



template <class double_type>
void Condensate<double_type>::calcActivity(
  const double temperature, const std::vector<Element<double_type>>& elements)
{
  if (!(temperature > phase_temp_limits[0] && temperature <= phase_temp_limits[1]))
  {
    log_activity = -10;
    return;
  }

  log_activity = mass_action_constant;

  for (auto & i : element_indices)
    log_activity += std::log(elements[i].number_density) * stoichiometric_vector[elements[i].index];

  //if (ln_activity < -9) ln_activity = -9;
  
  //if (std::abs(ln_activity) < 1e-10) ln_activity = 0;
  //if (ln_activity < -11) ln_activity = -11;
}



template <class double_type>
void Condensate<double_type>::calcMassActionConstant(const double temperature)
{
  const double_type thermal_energy = 1.0e6 / (CONST_K * temperature);

  double_type log_K = fit_coeff[0]/temperature
                    + fit_coeff[1]*std::log(temperature)
                    + fit_coeff[2]
                    + fit_coeff[3]*temperature
                    + fit_coeff[4]*temperature * temperature;

  double_type sigma = 0;

  for (auto & i : stoichiometric_vector)
    sigma += i;

  mass_action_constant = log_K - (sigma) * std::log(thermal_energy);
}



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



template <class double_type>
void Condensate<double_type>::degreeOfCondensation(
  const std::vector<Element<double_type>>& elements, const double_type total_element_density)
{
  if (reference_element == FASTCHEM_UNKNOWN_SPECIES)
    findReferenceElement(elements);

  degree_of_condensation = stoichiometric_vector[reference_element] 
    * number_density /(total_element_density * elements[reference_element].epsilon);
}



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
