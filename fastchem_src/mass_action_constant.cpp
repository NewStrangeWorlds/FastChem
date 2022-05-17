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


#include "fastchem_constants.h"
#include "species_struct.h"
#include "../_deps/boost_math-src/include/boost/math/special_functions/trunc.hpp"

#include <cmath>



namespace fastchem {


//Calculates the mass action constant, see Paper I, Eq. (2.9)
//Change this function if you want to implement your own parametrisation
template <class double_type>
void Molecule<double_type>::calcMassActionConstant(const double temperature)
{
  double_type log_K = 0;

  if (this->mass_action_const_tab.size() > 1)
    log_K = this->interpolateMassActionConstant(temperature);
  else
    log_K = mass_action_coeff[0]/temperature
                    + mass_action_coeff[1]*std::log(temperature)
                    + mass_action_coeff[2]
                    + mass_action_coeff[3]*temperature
                    + mass_action_coeff[4]*temperature * temperature;


  //adjusting log_K from its standard pressure (1 bar = 1e-6 dyn cm-2) to the actual pressure 
  double_type pressure_scaling = 1.0e-6 * CONST_K * temperature;
  mass_action_constant = log_K - sigma * std::log(pressure_scaling);
}


//Interpolates the mass action constant from tabulated data
//Uses a cubic b spline from the Boost library
template <class double_type>
double_type Molecule<double_type>::interpolateMassActionConstant(const double temperature)
{ 
  if (this->lnk_spline == nullptr)
    this->lnk_spline = new boost::math::interpolators::cardinal_cubic_b_spline<double_type>(this->mass_action_const_tab.begin(), 
                                                                                            this->mass_action_const_tab.end(), 
                                                                                            this->tab_temp_start, this->tab_temp_step);


  double_type log_K = this->lnk_spline->operator()(temperature);

  return log_K;
}


template struct Molecule<double>;
template struct Molecule<long double>;
}
