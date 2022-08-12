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
#include <functional>

#include "solver.h"

#include "../species_struct.h"
#include "../../_ext/Eigen/Dense"


namespace fastchem {


/*template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::assembleJacobian(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr,
  const size_t nb_cond_rem,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules)
{
  const size_t nb_elements = elements.size();
  const size_t nb_condensates_system = condensates.size() - nb_cond_rem;

  Eigen::MatrixXdt<double_type> jacobian;
  jacobian.setZero(nb_elements + nb_condensates_system, nb_elements + nb_condensates_system);


  for (size_t i=0; i<nb_condensates_system; ++i)
    jacobian(i, i) = - condensates[i]->activity_corr;


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
    }


  return jacobian;
}*/



template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::assembleJacobian(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr,
  const std::vector<double_type>& number_densities,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules)
{
  const size_t nb_condensates = condensates_jac.size();
  const size_t nb_elements = elements.size();

  Eigen::MatrixXdt<double_type> jacobian;
  jacobian.setZero(nb_elements + nb_condensates, nb_elements + nb_condensates);


  for (size_t i=0; i<nb_condensates; ++i)
    jacobian(i, i) = - activity_corr[condensates_jac[i]];


  for (size_t i=0; i<nb_elements; ++i)
  {
    jacobian(i+nb_condensates,i+nb_condensates) = elements[i]->number_density;

    for (size_t j=0; j<nb_elements; ++j)
    {
      for (auto l : elements[i]->molecule_list)
         jacobian(i+nb_condensates,j+nb_condensates) += molecules[l].stoichiometric_vector[elements[i]->index] * molecules[l].stoichiometric_vector[elements[j]->index] * molecules[l].number_density;

      for (auto l : condensates_rem)
        if (number_densities[l] > 0)
          jacobian(i+nb_condensates,j+nb_condensates) += condensates[l]->stoichiometric_vector[elements[i]->index] * condensates[l]->stoichiometric_vector[elements[j]->index] * number_densities[l] / activity_corr[l];
    }
  }


  for (size_t i=0; i<nb_condensates; ++i)
    for (size_t j=0; j<nb_elements; ++j)
    {
      jacobian(i, j+nb_condensates) = condensates[condensates_jac[i]]->stoichiometric_vector[elements[j]->index];
      jacobian(j+nb_condensates, i) = condensates[condensates_jac[i]]->stoichiometric_vector[elements[j]->index] * number_densities[condensates_jac[i]];
    }


  return jacobian;
}



template <class double_type>
Eigen::VectorXdt<double_type> CondPhaseSolver<double_type>::assembleRightHandSide(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_densities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const double_type total_element_density, 
      const double_type log_tau)
{
  const size_t nb_elements = elements.size();
  const size_t nb_cond_rem = condensates_rem.size();
  const size_t nb_cond_jac = condensates_jac.size();

  Eigen::VectorXdt<double_type> rhs_vector;
  rhs_vector.setZero(nb_elements + nb_cond_jac);


  for (size_t i=0; i<nb_cond_jac; ++i)
  { 
    const int index = condensates_jac[i];

    rhs_vector(i) = - condensates[index]->log_activity - activity_corr[index]
                    * (1.0 + log_tau - std::log(number_densities[index]) 
                    - std::log(activity_corr[index]));
  }


  for (size_t i=0; i<nb_elements; ++i)
  { 
    rhs_vector(i+nb_cond_jac) = total_element_density * elements[i]->epsilon - elements[i]->number_density;
    
    for (auto j : elements[i]->molecule_list)
      rhs_vector(i+nb_cond_jac) -= molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;
    
    for (size_t j=0; j<condensates.size(); ++j)
      rhs_vector(i+nb_cond_jac) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j];

    for (auto j : condensates_rem)
        rhs_vector(i+nb_cond_jac) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j]
                                  * (condensates[j]->log_activity/activity_corr[j] + log_tau - std::log(number_densities[j]) 
                                                   - std::log(activity_corr[j]) + 1.0);
  }

  return rhs_vector;
}



template <class double_type> std::vector<double_type> CondPhaseSolver<double_type>::solveSystem(
      Eigen::MatrixXdt<double_type>& coeff_matrix,
      Eigen::VectorXdt<double_type>& rhs)
{
  Eigen::PartialPivLU<Eigen::Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>> solver;
  //Eigen::FullPivHouseholderQR<Eigen::Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>> dense_solver;

  //dense_solver.setThreshold(1e-100);
  solver.compute(coeff_matrix);

  //std::cout << "treshhold: " << dense_solver.threshold() << "\n";
  Eigen::VectorXdt<double_type> result = solver.solve(rhs);

  //std::cout << "Matrix is invertible: " << dense_solver.isInvertible() << "\n";

  //if (!dense_solver.isInvertible())
    //result = pseudoNewtonSolve(jacobian, rhs_vector);

  bool a_solution_exists = (coeff_matrix*result).isApprox(rhs); 

  auto test = coeff_matrix*result;


  std::vector<double_type> result_vec(result.rows(), 0.0);
  //std::vector<double_type> result_vec(result);

  for (size_t i=0; i<result_vec.size(); ++i)
    result_vec[i] = result(i);

  return result_vec;
}



template class CondPhaseSolver<double>;
template class CondPhaseSolver<long double>;
}



