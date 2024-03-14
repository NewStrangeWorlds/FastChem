/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2024 Daniel Kitzmann, Joachim Stock
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


#include "solver.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

#include "../../_ext/Eigen/Dense"


namespace fastchem {

template <class double_type>
void GasPhaseSolver<double_type>::selectNewtonElements(
  std::vector<Element<double_type>>& elements,
  const std::vector<Molecule<double_type>>& molecules,
  const std::vector<double_type>& number_density_old,
  const double gas_density,
  std::vector<Element<double_type>*>& newton_elements)
{
  std::vector<int> select_elements(elements.size(), 0);

  for (size_t i=0; i<elements.size(); ++i)
  {
    if (std::fabs((elements[i].number_density - number_density_old[i])) > options.chem_accuracy*number_density_old[i])
    {
      if (elements[i].symbol != "e-" && elements[i].fixed_by_condensation == false)
        select_elements[i] = 1;
    }
  }


  for (size_t i=0; i<molecules.size(); ++i)
  {
    size_t j = i + elements.size();

    if (std::fabs((molecules[i].number_density - number_density_old[j])) > options.chem_accuracy*number_density_old[j]
            && molecules[i].number_density/gas_density > 1.e-155)
    {
      double_type max_abundance = 0;
      size_t max_index = 0;

      for (auto & e : molecules[i].element_indices)
      {
        if (elements[e].symbol == "e-" || elements[e].fixed_by_condensation) continue;

        if (elements[e].abundance > max_abundance)
        {
          max_index = e;
          max_abundance = elements[e].abundance;
        }
      }

      select_elements[max_index] = 1;
    }
  }


  for (size_t i=0; i<elements.size(); ++i)
  {
    if (select_elements[i] == 1) newton_elements.push_back(&elements[i]);
  }

}


template <class double_type>
Eigen::VectorXdt<double_type> GasPhaseSolver<double_type>::assembleJacobian(
  const std::vector<Element<double_type>*>& species,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules,
  Eigen::MatrixXdt<double_type>& jacobian)
{
  const size_t nb_newton_species = species.size();

  jacobian.setZero(nb_newton_species, nb_newton_species);

  for (size_t i=0; i<nb_newton_species; ++i)
  {
    jacobian(i,i) = species[i]->number_density;

    for (size_t j=0; j<nb_newton_species; ++j)
    {
      for (auto l : species[i]->molecule_list)
         jacobian(i,j) += 
           molecules[l].stoichiometric_vector[species[i]->index] 
           * molecules[l].stoichiometric_vector[species[j]->index] 
           * molecules[l].number_density;
    }
  }


  Eigen::VectorXdt<double_type> scaling_factors = jacobian.rowwise().maxCoeff();

  // for (auto i=0; i<jacobian.rows(); ++i)
  // {
  //   for (auto j=0; j<jacobian.rows(); ++j)
  //     jacobian(i,j) /= scaling_factors(i);
  // }

  return scaling_factors;
}


template <class double_type>
void GasPhaseSolver<double_type>::assembleRightHandSide(
  const std::vector<Element<double_type>*>& species,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules,
  const double total_element_density,
  const Eigen::VectorXdt<double_type>& scaling_factors,
  Eigen::VectorXdt<double_type>& rhs)
{
  const size_t nb_newton_species = species.size();

  rhs.setZero(nb_newton_species);

  for (size_t i=0; i<nb_newton_species; ++i)
  { 
    rhs(i) = total_element_density * species[i]->phi - species[i]->number_density;
    
    for (auto j : species[i]->molecule_list)
      rhs(i) -= molecules[j].stoichiometric_vector[species[i]->index] * molecules[j].number_density;
  }

  // for (auto i=0; i<rhs.rows(); ++i)
  //   rhs(i) /= scaling_factors(i);
}



template <class double_type>
void GasPhaseSolver<double_type>::newtonSolMult(
  std::vector<Element<double_type>*>& species,
  std::vector<Element<double_type>>& elements,
  const std::vector<Molecule<double_type>>& molecules, 
  const double_type gas_density)
{
  Eigen::MatrixXdt<double_type> jacobian;
  Eigen::VectorXdt<double_type> rhs;

  Eigen::VectorXdt<double_type> scaling_factors = assembleJacobian(
    species,
    elements,
    molecules,
    jacobian);

  assembleRightHandSide(
    species,
    elements,
    molecules,
    gas_density,
    scaling_factors,
    rhs);


  Eigen::PartialPivLU<Eigen::Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>> solver;
  solver.compute(jacobian);
  Eigen::VectorXdt<double_type> result = solver.solve(rhs);

  Eigen::VectorXdt<double_type> result_scaled = result;

  double_type max_value = result_scaled.cwiseAbs().maxCoeff();

  if (max_value > 2.0)
    result_scaled *= 2.0/max_value;

  for (size_t i=0; i<species.size(); ++i)
    species[i]->number_density = species[i]->number_density * std::exp(result_scaled[i]);
}


template class GasPhaseSolver<double>;
template class GasPhaseSolver<long double>;
}
