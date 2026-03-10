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

void GasPhaseSolver::selectNewtonElements(
  std::vector<Element>& elements,
  const std::vector<Molecule>& molecules,
  const std::vector<double>& number_density_old,
  const double total_element_density,
  std::vector<Element*>& newton_elements)
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
          && molecules[i].log_number_density > static_cast<double>(LOG_DENSITY_FLOOR) + 1)
    {
      double max_abundance = 0;
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
  { //if (!elements[i].fixed_by_condensation) select_elements[i] = 1;
    if (select_elements[i] == 1) newton_elements.push_back(&elements[i]);
  }
}


Eigen::VectorXdt GasPhaseSolver::assembleJacobian(
  const std::vector<Element*>& species,
  const std::vector< Element >& elements,
  const std::vector< Molecule >& molecules,
  Eigen::MatrixXdt& jacobian)
{
  const size_t nb_newton_species = species.size();

  jacobian.setZero(nb_newton_species, nb_newton_species);

  //Build global-element-index -> local-Newton-species-index map
  const size_t total_nb_elements = molecules.empty()
    ? 0 : molecules[0].stoichiometric_vector.size();
  std::vector<int> global_to_local(total_nb_elements, -1);
  for (size_t j=0; j<nb_newton_species; ++j)
    global_to_local[species[j]->index] = static_cast<int>(j);

  for (size_t i=0; i<nb_newton_species; ++i)
  {
    jacobian(i,i) = species[i]->number_density;

    for (auto l : species[i]->molecule_list)
    {
      const double n_l = molecules[l].number_density;
      if (n_l == 0.0) continue;
      const double scaled = static_cast<double>(
        molecules[l].stoichiometric_vector[species[i]->index]) * n_l;

      for (auto k_global : molecules[l].element_indices)
      {
        const int j = global_to_local[k_global];
        if (j < 0) continue;
        jacobian(i, j) += scaled * molecules[l].stoichiometric_vector[k_global];
      }
    }
  }


  Eigen::VectorXdt scaling_factors = jacobian.cwiseAbs().rowwise().maxCoeff();

  for (auto i = 0; i < scaling_factors.rows(); ++i)
    if (scaling_factors(i) == 0.0) scaling_factors(i) = 1.0;

  for (int j = 0; j < jacobian.cols(); ++j)
    jacobian.col(j).array() /= scaling_factors.array();

  return scaling_factors;
}


void GasPhaseSolver::assembleRightHandSide(
  const std::vector<Element*>& species,
  const std::vector< Element >& elements,
  const std::vector< Molecule >& molecules,
  const double total_element_density,
  const Eigen::VectorXdt& scaling_factors,
  Eigen::VectorXdt& rhs)
{
  const size_t nb_newton_species = species.size();

  rhs.setZero(nb_newton_species);

  for (size_t i=0; i<nb_newton_species; ++i)
  { 
    rhs(i) = total_element_density * species[i]->phi - species[i]->number_density;
    
    for (auto j : species[i]->molecule_list)
      rhs(i) -= molecules[j].stoichiometric_vector[species[i]->index] * molecules[j].number_density;
  }

  rhs.array() /= scaling_factors.array();
}



void GasPhaseSolver::newtonSolMult(
  std::vector<Element*>& species,
  std::vector<Element>& elements,
  const std::vector<Molecule>& molecules, 
  const double total_element_density)
{
  Eigen::MatrixXdt jacobian;
  Eigen::VectorXdt rhs;

  Eigen::VectorXdt scaling_factors = assembleJacobian(
    species,
    elements,
    molecules,
    jacobian);

  assembleRightHandSide(
    species,
    elements,
    molecules,
    total_element_density,
    scaling_factors,
    rhs);


  Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> solver;
  solver.compute(jacobian);
  Eigen::VectorXdt result = solver.solve(rhs);

  Eigen::VectorXdt result_scaled = result;

  double max_value = result_scaled.cwiseAbs().maxCoeff();

  if (max_value > 2.0)
    result_scaled *= 2.0/max_value;

  for (size_t i=0; i<species.size(); ++i)
  { 
    species[i]->log_number_density += result_scaled[i];
    species[i]->number_density = safeExp(species[i]->log_number_density);
  }
}


}
