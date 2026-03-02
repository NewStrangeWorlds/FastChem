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


#include "../fastchem_constants.h"
#include "../species_struct.h"

#include <iostream>
#include <cmath>
#include <limits>


namespace fastchem {


//Calculates the mass action constant, see Paper I, Eq. (2.9)
//Change this function if you want to implement your own parametrisation
void Molecule::calcMassActionConstant(const double temperature)
{
  double log_K = mass_action_coeff[0]/temperature
                    + mass_action_coeff[1]*std::log(temperature)
                    + mass_action_coeff[2]
                    + mass_action_coeff[3]*temperature
                    + mass_action_coeff[4]*temperature * temperature;
  
  //adjusting log_K from its standard pressure (1 bar = 1e-6 dyn cm-2) to the actual pressure 
  const double pressure_scaling = 1.0e-6 * CONST_K * temperature;
  mass_action_constant = log_K - sigma * std::log(pressure_scaling);
}



//Check for the number density of molecules
void Molecule::checkN(
  const double& min_limit, const double& gas_density)
{
  if (this->log_number_density < static_cast<double>(LOG_DENSITY_FLOOR))
    this->log_number_density = static_cast<double>(LOG_DENSITY_FLOOR);

  const double log_gas = std::log(gas_density);

  if (this->log_number_density > log_gas)
    this->log_number_density = log_gas;

  this->number_density = safeExp(this->log_number_density);
}


//Compute log(n_i) = mac + sum_l nu_il * y_l (pure log-space arithmetic)
void Molecule::calcLogNumberDensity(const std::vector< Element >& elements)
{
  this->log_number_density = mass_action_constant;

  for (auto i : element_indices)
    this->log_number_density += stoichiometric_vector[i] * elements[i].log_number_density;
}


void Molecule::calcNumberDensity(const std::vector< Element >& elements)
{
  calcLogNumberDensity(elements);
  this->number_density = safeExp(this->log_number_density);
}



}
