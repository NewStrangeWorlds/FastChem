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


template <class double_type>
bool CondPhaseSolver<double_type>::newtonStep(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules,
  const double_type total_element_density,
  const std::vector<double_type>& cond_densities,
  const std::vector<double_type>& element_densities,
  const std::vector<double_type>& activity_corr,
  Eigen::VectorXdt<double_type>& result,
  Eigen::VectorXdt<double_type>& scaling_factors,
  double_type& objective_function)
{
  Eigen::MatrixXdt<double_type> jacobian = assembleJacobian(
    condensates,
    activity_corr,
    cond_densities,
    condensates_jac,
    condensates_rem,
    elements,
    molecules,
    total_element_density);


  scaling_factors = jacobian.rowwise().maxCoeff();

  for (auto i=0; i<jacobian.rows(); ++i)
  {
    for (auto j=0; j<jacobian.rows(); ++j)
      jacobian(i,j) /= scaling_factors(i);
  }


  Eigen::VectorXdt<double_type> rhs = assembleRightHandSideLog(
    condensates,
    condensates_jac,
    condensates_rem,
    activity_corr,
    cond_densities,
    elements,
    molecules,
    total_element_density,
    0,
    scaling_factors,
    objective_function);


  const bool jacobian_is_invertible = solveSystem(jacobian, rhs, result);

  return jacobian_is_invertible; 
}



template <class double_type>
std::vector<double_type> CondPhaseSolver<double_type>::newtonStepFull(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules,
  const double_type total_element_density,
  const std::vector<double_type>& cond_densities,
  const std::vector<double_type>& element_densities,
  const std::vector<double_type>& activity_corr,
  Eigen::VectorXdt<double_type>& scaling_factors,
  double_type& objective_function)
{

  Eigen::MatrixXdt<double_type> jacobian = assembleJacobianFullLog(
    condensates,
    activity_corr,
    cond_densities,
    elements,
    molecules,
    total_element_density);

  scaling_factors = jacobian.rowwise().maxCoeff();

  Eigen::VectorXdt<double_type> rhs = assembleRightHandSideFullLog(
    condensates,
    activity_corr,
    cond_densities,
    elements,
    molecules,
    total_element_density,
    scaling_factors,
    objective_function);


  //for (auto i=condensates_jac.size(); i<rhs.rows(); ++i)
  for (auto i=0; i<rhs.rows(); ++i)
  {
    for (auto j=0; j<rhs.rows(); ++j)
      jacobian(i,j) /= scaling_factors(i);
  }

  //std::cout << jacobian << "\n";
  //std::cout << rhs << "\n";

  Eigen::VectorXdt<double_type> result;
  const bool jacobian_is_invertible = solveSystem(jacobian, rhs, result);

  std::vector<double_type> res;
  return res;
  //return jacobian_is_invertible; 
}


template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::assembleJacobian(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr,
  const std::vector<double_type>& number_densities,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules,
  const double_type total_element_density)
{
  const size_t nb_condensates = condensates_jac.size();
  const size_t nb_elements = elements.size();

  Eigen::MatrixXdt<double_type> jacobian;
  jacobian.setZero(nb_elements + nb_condensates, nb_elements + nb_condensates);


  for (size_t i=0; i<nb_condensates; ++i)
  {
    jacobian(i, i) = - activity_corr[condensates_jac[i]];
  }
   


  for (size_t i=0; i<nb_elements; ++i)
  {
    jacobian(i+nb_condensates,i+nb_condensates) = elements[i]->number_density;

    for (size_t j=0; j<nb_elements; ++j)
    {
      for (auto l : elements[i]->molecule_list)
         jacobian(i+nb_condensates,j+nb_condensates) += 
           molecules[l].stoichiometric_vector[elements[i]->index] 
           * molecules[l].stoichiometric_vector[elements[j]->index] 
           * molecules[l].number_density;

      for (auto l : condensates_rem)
        if (number_densities[l] > 0)
          jacobian(i+nb_condensates,j+nb_condensates) += 
            condensates[l]->stoichiometric_vector[elements[i]->index] 
            * condensates[l]->stoichiometric_vector[elements[j]->index] 
            * number_densities[l] / activity_corr[l];
    }
  }


  for (size_t i=0; i<nb_condensates; ++i)
    for (size_t j=0; j<nb_elements; ++j)
    {
      jacobian(i, j+nb_condensates) = 
        condensates[condensates_jac[i]]->stoichiometric_vector[elements[j]->index];
      
      jacobian(j+nb_condensates, i) = 
        condensates[condensates_jac[i]]->stoichiometric_vector[elements[j]->index]
        * number_densities[condensates_jac[i]];
    }


  return jacobian;
}



template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::assembleJacobianFull(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr,
  const std::vector<double_type>& number_densities,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules,
  const double_type total_element_density)
{
  const size_t nb_condensates = condensates.size();
  const size_t nb_elements = elements.size();

  Eigen::MatrixXdt<double_type> jacobian;
  jacobian.setZero(nb_elements + 2*nb_condensates, nb_elements + 2*nb_condensates);

  for (size_t i=0; i<nb_condensates; ++i)
  {
    jacobian(i, i) = activity_corr[i]*number_densities[i];
    jacobian(i, i+nb_condensates) = number_densities[i];
    jacobian(i+nb_condensates, i+nb_condensates) = 1.0;
  }


  for (size_t i=0; i<nb_elements; ++i)
  {
    jacobian(i+2*nb_condensates,i+2*nb_condensates) = elements[i]->number_density;

    for (size_t j=0; j<nb_elements; ++j)
    {
      for (auto l : elements[i]->molecule_list)
         jacobian(i+2*nb_condensates,j+2*nb_condensates) += 
           molecules[l].stoichiometric_vector[elements[i]->index] 
           * molecules[l].stoichiometric_vector[elements[j]->index] 
           * molecules[l].number_density;
    }
  }

  for (size_t i=0; i<nb_condensates; ++i)
    for (size_t j=0; j<nb_elements; ++j)
    {
      jacobian(i+nb_condensates, j+2*nb_condensates) = 
        condensates[i]->stoichiometric_vector[elements[j]->index];
      
      jacobian(j+2*nb_condensates, i) = 
        condensates[i]->stoichiometric_vector[elements[j]->index] 
        * number_densities[i];
    }

  return jacobian;
}



template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::assembleJacobianFullLog(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr,
  const std::vector<double_type>& number_densities,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<Molecule<double_type>>& molecules,
  const double_type total_element_density)
{
  const size_t nb_condensates = condensates.size();
  const size_t nb_elements = elements.size();

  Eigen::MatrixXdt<double_type> jacobian;
  jacobian.setZero(nb_elements + 2*nb_condensates, nb_elements + 2*nb_condensates);

  for (size_t i=0; i<nb_condensates; ++i)
  {
    jacobian(i, i) = 1.0;
    jacobian(i, i+nb_condensates) = 1.0;
    jacobian(i+nb_condensates, i+nb_condensates) = activity_corr[i];
  }


  for (size_t i=0; i<nb_elements; ++i)
  {
    jacobian(i+2*nb_condensates,i+2*nb_condensates) = elements[i]->number_density;

    for (size_t j=0; j<nb_elements; ++j)
    {
      for (auto l : elements[i]->molecule_list)
         jacobian(i+2*nb_condensates,j+2*nb_condensates) += 
           molecules[l].stoichiometric_vector[elements[i]->index] 
           * molecules[l].stoichiometric_vector[elements[j]->index] 
           * molecules[l].number_density;
    }
  }

  for (size_t i=0; i<nb_condensates; ++i)
    for (size_t j=0; j<nb_elements; ++j)
    {
      jacobian(i+nb_condensates, j+2*nb_condensates) = 
        condensates[i]->stoichiometric_vector[elements[j]->index];
      
      jacobian(j+2*nb_condensates, i) = 
        condensates[i]->stoichiometric_vector[elements[j]->index] 
        * number_densities[i];
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
      const double_type log_tau,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      double_type& objective_function)
{
  const size_t nb_elements = elements.size();
  const size_t nb_cond_rem = condensates_rem.size();
  const size_t nb_cond_jac = condensates_jac.size();

  Eigen::VectorXdt<double_type> rhs_vector;
  rhs_vector.setZero(nb_elements + nb_cond_jac);


  for (size_t i=0; i<nb_cond_jac; ++i)
  { 
    const int index = condensates_jac[i];

    rhs_vector(i) = - condensates[index]->log_activity - condensates[index]->tau / number_densities[index];
  }


  for (size_t i=0; i<nb_elements; ++i)
  { 
    rhs_vector(i+nb_cond_jac) = total_element_density * elements[i]->epsilon - elements[i]->number_density;
    
    for (auto j : elements[i]->molecule_list)
      rhs_vector(i+nb_cond_jac) -= molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;
    
    for (size_t j=0; j<condensates.size(); ++j)
      rhs_vector(i+nb_cond_jac) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j];

    for (auto j : condensates_rem)
        rhs_vector(i+nb_cond_jac) -= 
          condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j]
          * (condensates[j]->log_activity/activity_corr[j] 
          + condensates[j]->tau/(number_densities[j] * activity_corr[j]));
  }


  for (auto i=0; i<rhs_vector.rows(); ++i)
    rhs_vector(i) /= scaling_factors(i);

  objective_function = 0.5*rhs_vector.transpose()*rhs_vector;

  return rhs_vector;
}


template <class double_type>
Eigen::VectorXdt<double_type> CondPhaseSolver<double_type>::assembleRightHandSideLog(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_densities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const double_type total_element_density,
      const double_type log_tau,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      double_type& objective_function)
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
                     * (1.0 + condensates[index]->log_tau - std::log(number_densities[index]) 
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
        rhs_vector(i+nb_cond_jac) -= 
          condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j]
          * (condensates[j]->log_activity/activity_corr[j] 
          + condensates[j]->log_tau - std::log(number_densities[j]) -  std::log(activity_corr[j]) + 1.0);
  }


  for (auto i=0; i<rhs_vector.rows(); ++i)
    rhs_vector(i) /= scaling_factors(i);

  objective_function = 0.5*rhs_vector.transpose()*rhs_vector;

  return rhs_vector;
}


template <class double_type>
Eigen::VectorXdt<double_type> CondPhaseSolver<double_type>::assembleRightHandSideFull(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_densities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const double_type total_element_density,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      double_type& objective_function)
{
  const size_t nb_elements = elements.size();
  const size_t nb_cond = condensates.size();

  Eigen::VectorXdt<double_type> rhs_vector;
  rhs_vector.setZero(nb_elements + 2*nb_cond);


  for (size_t i=0; i<nb_cond; ++i)
  { 
    rhs_vector(i) = condensates[i]->tau - activity_corr[i]*number_densities[i]; 
    rhs_vector(i+nb_cond) = - condensates[i]->log_activity - activity_corr[i];
  }


  for (size_t i=0; i<nb_elements; ++i)
  { 
    rhs_vector(i+2*nb_cond) = total_element_density * elements[i]->epsilon - elements[i]->number_density;
    
    for (auto j : elements[i]->molecule_list)
      rhs_vector(i+2*nb_cond) -= molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;
    
    for (size_t j=0; j<condensates.size(); ++j)
      rhs_vector(i+2*nb_cond) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j];
  }


  for (auto i=0; i<rhs_vector.rows(); ++i)
    rhs_vector(i) /= scaling_factors(i);

  objective_function = 0.5*rhs_vector.transpose()*rhs_vector;

  return rhs_vector;
}


template <class double_type>
Eigen::VectorXdt<double_type> CondPhaseSolver<double_type>::assembleRightHandSideFullLog(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_densities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const double_type total_element_density,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      double_type& objective_function)
{
  const size_t nb_elements = elements.size();
  const size_t nb_cond = condensates.size();

  Eigen::VectorXdt<double_type> rhs_vector;
  rhs_vector.setZero(nb_elements + 2*nb_cond);


  for (size_t i=0; i<nb_cond; ++i)
  { 
    rhs_vector(i) = condensates[i]->log_tau - std::log(activity_corr[i]) - std::log(number_densities[i]); 
    rhs_vector(i+nb_cond) = - condensates[i]->log_activity - activity_corr[i];
  }


  for (size_t i=0; i<nb_elements; ++i)
  { 
    rhs_vector(i+2*nb_cond) = total_element_density * elements[i]->epsilon - elements[i]->number_density;
    
    for (auto j : elements[i]->molecule_list)
      rhs_vector(i+2*nb_cond) -= molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;
    
    for (size_t j=0; j<condensates.size(); ++j)
      rhs_vector(i+2*nb_cond) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j];
  }


  for (auto i=0; i<rhs_vector.rows(); ++i)
    rhs_vector(i) /= scaling_factors(i);

  objective_function = 0.5*rhs_vector.transpose()*rhs_vector;

  return rhs_vector;
}



template <class double_type>
double_type CondPhaseSolver<double_type>::objectiveFunction(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_densities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      const double_type total_element_density)
{
  double_type objective_function = 0;

  Eigen::VectorXdt<double_type> rhs = assembleRightHandSide(
    condensates,
    condensates_jac,
    condensates_rem,
    activity_corr,
    number_densities,
    elements,
    molecules,
    total_element_density,
    0,
    scaling_factors,
    objective_function);

  return objective_function;
}


template <class double_type>
double_type CondPhaseSolver<double_type>::objectiveFunctionFull(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_densities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      const double_type total_element_density)
{
  double_type objective_function = 0;

  Eigen::VectorXdt<double_type> rhs = assembleRightHandSideFullLog(
    condensates,
    activity_corr,
    number_densities,
    elements,
    molecules,
    total_element_density,
    scaling_factors,
    objective_function);

  return objective_function;
}


template <class double_type>
Eigen::MatrixXdt<double_type> CondPhaseSolver<double_type>::computePerturbedHessian(
  const Eigen::MatrixXdt<double_type>& jacobian,
  const double_type perturbation)
{
  Eigen::MatrixXdt<double_type> hessian = jacobian.transpose()*jacobian;

  const double_type norm = hessian.template lpNorm<1>();

  for (auto i=0; i<hessian.rows(); ++i)
    hessian(i,i) += perturbation * norm; 

  return hessian;
}




template <class double_type> bool CondPhaseSolver<double_type>::solveSystem(
  const Eigen::MatrixXdt<double_type>& jacobian,
  const Eigen::VectorXdt<double_type>& rhs,
  Eigen::VectorXdt<double_type>& result)
{

  if (!options.cond_use_full_pivot)
  {
    Eigen::PartialPivLU<Eigen::Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>> solver;
    solver.compute(jacobian);
    result = solver.solve(rhs);

    return true;
  }


  if (options.cond_use_full_pivot)
  {
    Eigen::FullPivLU<Eigen::Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>> solver;
    solver.compute(jacobian);
    result = solver.solve(rhs);

    if (!solver.isInvertible())
    {
      std::cout << "FastChem warning: Jacobian is (almost) singular! ";
      
      if (options.cond_use_svd)
      {
        std::cout << "Switching to Single Value Decomposition.\n";
        Eigen::BDCSVD<Eigen::Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>> solver;
        Eigen::VectorXdt<double_type> result = jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
      }
      else
      {
        std::cout << "Switching to perturbed Hessian approximation.\n";
        const double_type perturbation = std::numeric_limits<double_type>::epsilon() * 10;

        Eigen::MatrixXdt<double_type> hessian = computePerturbedHessian(
          jacobian,
          perturbation);

        solver.compute(hessian);
        Eigen::VectorXdt<double_type> rhs_pertubed = jacobian.transpose()*rhs;

        result = solver.solve(rhs_pertubed);
      }

      return false;
    }

    return true;
  }

}



template class CondPhaseSolver<double>;
template class CondPhaseSolver<long double>;
}



