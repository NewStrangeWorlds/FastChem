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


#include "species_struct.h"

#include <cmath>
#include <iostream>


namespace fastchem {


//Determination of the correction factor n_j_min, that contains the number density of an element j contained in molecules
//with elements that are less abundant than j
//See Eq. (2.23)
template <class double_type>
void Element<double_type>::calcMinorSpeciesDensities(
  const std::vector< Molecule<double_type> > &molecules)
{
  number_density_min = 0.0;

  for (auto & i : minor_molecules)
    number_density_min += (molecules[i].stoichiometric_vector[index] 
                       + epsilon * molecules[i].sigma) * molecules[i].number_density;
}



template <class double_type>
void Element<double_type>::calcEpsilon(const std::vector< Element<double_type> > &elements)
{ 
  double_type element_sum = 0.0;

  for (auto & i : elements)
    element_sum += i.abundance;

  epsilon = this->abundance/element_sum;
}



template <class double_type>
void Element<double_type>::calcSolverScalingFactor(
  const std::vector<Element<double_type>> &elements, 
  const std::vector<Molecule<double_type>> &molecules,
  const double additional_scaling_factor)
{
  solver_scaling_factor = 0.0;


  for (auto & i : molecule_list)
  {
    if (molecules[i].stoichiometric_vector[index] < 1 || molecules[i].stoichiometric_vector[index] > static_cast<int>(solver_order) )
      continue;

    double_type sum = 0.0;

    if (molecules[i].abundance == this->abundance)
    {
      for (auto & l : molecules[i].element_indices)
      {
        if (l != index)
          sum += molecules[i].stoichiometric_vector[l] * std::log(elements[l].number_density);

      }

      sum += molecules[i].mass_action_constant;
    }


    if (sum > solver_scaling_factor)
      solver_scaling_factor = sum;
  }

  //scale the factor by an additional, optional factor supplied by the user
  solver_scaling_factor -= additional_scaling_factor;
}


template struct Element<double>;
template struct Element<long double>;
template struct Molecule<double>;
template struct Molecule<long double>;
}
