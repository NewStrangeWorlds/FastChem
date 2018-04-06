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

#include <vector>
#include <cmath>



namespace fastchem {


//Iteration for the n_<H> density
//See Sect. 2.3 for details
template <class double_type>
bool FastChem<double_type>::calcTotalHydrogenDensity(const double temperature_gas, const double pressure, const unsigned int grid_index,
                                                     double_type& h_density, double_type& density_iteration_lambda, double_type& density_iteration_error)
{
  //this value will be fixed
  double_type total_density = pressure/(CONST_K * temperature_gas);

  double_type total_density_calc = 0;

  for (size_t i=0; i<nb_species; ++i)
    total_density_calc += species[i]->number_density[grid_index];


  double_type current_total_density_error = (total_density - total_density_calc)/total_density;

  if (density_iteration_error * current_total_density_error < 0)
    density_iteration_lambda = 0.1 * density_iteration_lambda + 0.9;


  bool is_converged = false;

  if (std::fabs(current_total_density_error) < accuracy_delta)
     is_converged = true;
  else
  {
    if (total_density_calc > total_density)
      h_density = density_iteration_lambda * h_density;
    else
      h_density = 1./density_iteration_lambda * h_density;

    is_converged = false;
  }


  density_iteration_error = current_total_density_error;


  if (std::isnan(total_density_calc)) h_density = total_density_calc;


  return is_converged;
}



//Alternative version for the calculation of the n_<H> density
//Uses the approach of Woitke et al. (2018), can be faster than the standard function
template <class double_type>
bool FastChem<double_type>::calcTotalHydrogenDensityAlt(const double temperature_gas, const double pressure, const unsigned int grid_index,
                                                        double_type& h_density, double_type& muH, double_type& density_iteration_error)
{
  //this value will be fixed
  double_type total_density = pressure/(CONST_K * temperature_gas);

  double_type total_density_calc = 0;

  for (size_t i=0; i<nb_species; ++i)
    total_density_calc += species[i]->number_density[grid_index];


  double_type current_total_density_error = (total_density - total_density_calc)/total_density;


  bool is_converged = false;

  if (std::fabs(current_total_density_error) < accuracy_delta)
     is_converged = true;
  else
  {
    double_type pressure_calc = total_density_calc * CONST_K * temperature_gas;

    double_type mu = h_density/pressure_calc * (CONST_K*temperature_gas)*muH;
    h_density = pressure * mu/(CONST_K * temperature_gas)/muH;


    is_converged = false;
  }



  density_iteration_error = current_total_density_error;


  if (std::isnan(total_density_calc)) h_density = total_density_calc;


  return is_converged;
}


template class FastChem<double>;
template class FastChem<long double>;


}
