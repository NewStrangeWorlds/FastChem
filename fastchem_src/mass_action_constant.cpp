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

#include <cmath>



namespace fastchem {


//Calculates the mass action constant, see Eq. (2.9)
//Change this function if you want to implement your own parametrisation
template <class double_type>
void Molecule<double_type>::calcMassActionConstant(const double temperature, const unsigned int grid_index)
{
  double_type thermal_energy = 1.0e-6 * CONST_K * temperature;

  double_type log_K = mass_action_coeff[0]/temperature
                    + mass_action_coeff[1]*std::log(temperature)
                    + mass_action_coeff[2]
                    + mass_action_coeff[3]*temperature
                    + mass_action_coeff[4]*temperature * temperature;

  mass_action_constant[grid_index] = log_K - sigma * std::log(thermal_energy);
}


template struct Molecule<double>;
template struct Molecule<long double>;

}
