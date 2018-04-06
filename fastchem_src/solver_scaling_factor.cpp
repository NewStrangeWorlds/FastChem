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

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//Computes the scaling factor psi to avoid numerical overflow
//See Appendix A for details
template <class double_type>
double_type FastChem<double_type>::solverScalingFactor(Element<double_type>& species, const double_type number_density_min, const double_type h_density, const unsigned int grid_index)
{
  double_type scaling_factor = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[species.index] < 1 || molecules[i].stoichometric_vector[species.index] > static_cast<int>(species.solver_order) )
      continue;


    if (molecules[i].abundance == species.abundance)
    {
      molecules[i].sum[grid_index] = 0.0;


      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];

        if (l != species.index)
          molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density[grid_index]);

      }


      molecules[i].sum[grid_index] += molecules[i].mass_action_constant[grid_index];
    }


    if (molecules[i].sum[grid_index] > scaling_factor)
      scaling_factor = molecules[i].sum[grid_index];

  }



  double_type xi = (number_density_min + species.number_density[grid_index]) * std::exp(-scaling_factor);

  for (size_t i=0; i<species.molecule_list.size(); ++i)
  {
    unsigned int j = species.molecule_list[i];

    if (species.abundance == molecules[j].abundance)
      xi += molecules[j].stoichometric_vector[species.index] * std::pow(molecules[j].number_density[grid_index], molecules[j].stoichometric_vector[species.index]);
  }


  if (xi == 0)
    xi = std::numeric_limits<double_type>::max_exponent / 6.0;
  else
    xi = std::numeric_limits<double_type>::max_exponent - std::log(xi);

  xi = std::sqrt(xi);


  scaling_factor -= xi;


  return scaling_factor;
}




template class FastChem<double>;
template class FastChem<long double>;

}



