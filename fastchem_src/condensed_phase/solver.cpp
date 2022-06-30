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


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

#include "solver.h"

#include "../species_struct.h"
#include "../../_ext/Eigen/Dense"


namespace fastchem {


template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::assembleJacobian(
  const std::vector<Condensate<double_type>*>& condensates,
  const size_t nb_cond_rem,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules)
{
  const size_t nb_elements = elements.size();
  const size_t nb_condensates_system = condensates.size() - nb_cond_rem;

  Eigen::MatrixXdt<double_type> jacobian;
  jacobian.setZero(nb_elements + nb_condensates_system, nb_elements + nb_condensates_system);


  /*for (size_t i=0; i<nb_condensates_system; ++i)
    jacobian(i, i) = - condensates[i]->activity_correction;


  for (size_t i=0; i<nb_elements; ++i)
  {
    double_type sum = elements[i]->number_density;
    
    for (auto j : elements[i]->molecule_list)
      sum += molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;

    for (auto c : condensates)
      sum += c->stoichiometric_vector[elements[i]->index] * c->number_density;

    double_type element_sum = sum;


    jacobian(i+nb_condensates_system,i+nb_condensates_system) = elements[i]->number_density / element_sum;

    for (size_t j=0; j<nb_elements; ++j)
    {
      for (auto l : elements[i]->molecule_list)
         jacobian(i+nb_condensates_system,j+nb_condensates_system) += molecules[l].stoichiometric_vector[elements[i]->index] * molecules[l].stoichiometric_vector[elements[j]->index] * molecules[l].number_density / element_sum;

      for (auto l : condensates)
        if (l->linear_system_remove)
          jacobian(i+nb_condensates_system,j+nb_condensates_system) += l->stoichiometric_vector[elements[i]->index] * l->stoichiometric_vector[elements[j]->index] * l->number_density / l->activity_correction / element_sum;
    }
  }


  for (size_t i=0; i<nb_condensates_system; ++i)
    for (size_t j=0; j<nb_elements; ++j)
    { 
      double_type sum = elements[j]->number_density;
    
      for (auto k : elements[j]->molecule_list)
        sum += molecules[k].stoichiometric_vector[elements[j]->index] * molecules[k].number_density;

      for (auto c : condensates)
        sum += c->stoichiometric_vector[elements[j]->index] * c->number_density;

      double_type element_sum = sum;

      jacobian(i, j+nb_condensates_system) = condensates[i]->stoichiometric_vector[elements[j]->index];
      jacobian(j+nb_condensates_system, i) = condensates[i]->stoichiometric_vector[elements[j]->index] * condensates[i]->number_density / element_sum;
    }*/


  return jacobian;
}


template class CondPhaseSolver<double>;
template class CondPhaseSolver<long double>;
}



