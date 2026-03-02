/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2022 Daniel Kitzmann, Joachim Stock
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


bool CondPhaseSolver::newtonStep(
  const std::vector<Condensate*>& condensates,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<Element*>& elements,
  const std::vector<Molecule>& molecules,
  const double total_element_density,
  const std::vector<double>& cond_densities,
  const std::vector<double>& element_densities,
  const std::vector<double>& activity_corr,
  const std::vector<double>& log_cond_densities,
  const std::vector<double>& log_activity_corr,
  Eigen::VectorXdt& result,
  Eigen::VectorXdt& scaling_factors,
  double& objective_function)
{
  Eigen::MatrixXdt jacobian;
  scaling_factors = assembleJacobian(
    condensates,
    activity_corr,
    cond_densities,
    condensates_jac,
    condensates_rem,
    elements,
    molecules,
    total_element_density,
    jacobian);

  Eigen::VectorXdt rhs;
  objective_function = assembleRightHandSide(
    condensates,
    condensates_jac,
    condensates_rem,
    activity_corr,
    cond_densities,
    log_cond_densities,
    log_activity_corr,
    elements,
    molecules,
    total_element_density,
    scaling_factors,
    rhs);


  const bool jacobian_is_invertible = solveSystem(jacobian, rhs, result);

  return jacobian_is_invertible;
}



bool CondPhaseSolver::newtonStepFull(
  const std::vector<Condensate*>& condensates,
  const std::vector<Element*>& elements,
  const std::vector<Molecule>& molecules,
  const double total_element_density,
  const std::vector<double>& cond_densities,
  const std::vector<double>& element_densities,
  const std::vector<double>& activity_corr,
  const std::vector<double>& log_cond_densities,
  const std::vector<double>& log_activity_corr,
  Eigen::VectorXdt& result,
  Eigen::VectorXdt& scaling_factors,
  double& objective_function)
{
  Eigen::MatrixXdt jacobian;
  scaling_factors = assembleJacobianFull(
    condensates,
    activity_corr,
    cond_densities,
    elements,
    molecules,
    total_element_density,
    jacobian);

  Eigen::VectorXdt rhs;
  objective_function= assembleRightHandSideFull(
    condensates,
    activity_corr,
    cond_densities,
    log_cond_densities,
    log_activity_corr,
    elements,
    molecules,
    total_element_density,
    scaling_factors,
    rhs);


  const bool jacobian_is_invertible = solveSystem(jacobian, rhs, result);

  return jacobian_is_invertible;
}


//The reduced Jacobian, see Eq. 43, Paper III
Eigen::VectorXdt CondPhaseSolver::assembleJacobian(
  const std::vector<Condensate*>& condensates,
  const std::vector<double>& activity_corr,
  const std::vector<double>& number_densities,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<Element*>& elements,
  const std::vector<Molecule>& molecules,
  const double total_element_density,
  Eigen::MatrixXdt& jacobian)
{
  const size_t nb_condensates = condensates_jac.size();
  const size_t nb_elements = elements.size();

  jacobian.setZero(nb_elements + nb_condensates, nb_elements + nb_condensates);


  for (size_t i=0; i<nb_condensates; ++i)
  {
    jacobian(i, i) = - activity_corr[condensates_jac[i]];

    for (size_t j=0; j<nb_elements; ++j)
    {
      jacobian(i, j+nb_condensates) = 
        condensates[condensates_jac[i]]->stoichiometric_vector[elements[j]->index];
      
      jacobian(j+nb_condensates, i) = 
        condensates[condensates_jac[i]]->stoichiometric_vector[elements[j]->index]
        * number_densities[condensates_jac[i]];
    }
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


  Eigen::VectorXdt scaling_factors = jacobian.cwiseAbs().rowwise().maxCoeff();

  for (auto i = 0; i < scaling_factors.rows(); ++i)
    if (scaling_factors(i) == 0.0) scaling_factors(i) = 1.0;

  for (auto i=0; i<jacobian.rows(); ++i)
  {
    for (auto j=0; j<jacobian.rows(); ++j)
      jacobian(i,j) /= scaling_factors(i);
  }

  return scaling_factors;
}


//Full Jacobian, see Eq. 33, Paper III
Eigen::VectorXdt CondPhaseSolver::assembleJacobianFull(
  const std::vector<Condensate*>& condensates,
  const std::vector<double>& activity_corr,
  const std::vector<double>& number_densities,
  const std::vector<Element*>& elements,
  const std::vector<Molecule>& molecules,
  const double total_element_density,
  Eigen::MatrixXdt& jacobian)
{
  const size_t nb_condensates = condensates.size();
  const size_t nb_elements = elements.size();

  jacobian.setZero(nb_elements + 2*nb_condensates, nb_elements + 2*nb_condensates);

  for (size_t i=0; i<nb_condensates; ++i)
  {
    jacobian(i, i) = 1.0;
    jacobian(i, i+nb_condensates) = 1.0;
    jacobian(i+nb_condensates, i+nb_condensates) = activity_corr[i];

    for (size_t j=0; j<nb_elements; ++j)
    {
      jacobian(i+nb_condensates, j+2*nb_condensates) = 
        condensates[i]->stoichiometric_vector[elements[j]->index];
      
      jacobian(j+2*nb_condensates, i) = 
        condensates[i]->stoichiometric_vector[elements[j]->index] 
        * number_densities[i];
    }
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


  Eigen::VectorXdt scaling_factors = jacobian.cwiseAbs().rowwise().maxCoeff();

  for (auto i = 0; i < scaling_factors.rows(); ++i)
    if (scaling_factors(i) == 0.0) scaling_factors(i) = 1.0;

  for (auto i=0; i<jacobian.rows(); ++i)
  {
    for (auto j=0; j<jacobian.rows(); ++j)
      jacobian(i,j) /= scaling_factors(i);
  }

  return scaling_factors;
}


//see Eq. 43, Paper III
double CondPhaseSolver::assembleRightHandSide(
      const std::vector<Condensate*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double>& activity_corr,
      const std::vector<double>& number_densities,
      const std::vector<double>& log_number_densities,
      const std::vector<double>& log_activity_corr,
      const std::vector< Element* >& elements,
      const std::vector< Molecule >& molecules,
      const double total_element_density,
      const Eigen::VectorXdt& scaling_factors,
      Eigen::VectorXdt& rhs)
{
  const size_t nb_cond_jac = condensates_jac.size();

  rhs.setZero(elements.size() + nb_cond_jac);

  for (size_t i=0; i<nb_cond_jac; ++i)
  {
    const int index = condensates_jac[i];

    rhs(i) = - condensates[index]->log_activity - activity_corr[index]
            * (1.0 + condensates[index]->log_tau - log_number_densities[index]
            - log_activity_corr[index]);
  }


  for (size_t i=0; i<elements.size(); ++i)
  {
    rhs(i+nb_cond_jac) = total_element_density * elements[i]->epsilon - elements[i]->number_density;

    for (auto j : elements[i]->molecule_list)
      rhs(i+nb_cond_jac) -= molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;

    for (size_t j=0; j<condensates.size(); ++j)
      rhs(i+nb_cond_jac) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j];

    for (auto j : condensates_rem)
        rhs(i+nb_cond_jac) -=
          condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j]
          * (condensates[j]->log_activity/activity_corr[j]
          + condensates[j]->log_tau - log_number_densities[j] - log_activity_corr[j] + 1.0);
  }


  for (auto i=0; i<rhs.rows(); ++i)
    rhs(i) /= scaling_factors(i);


  const double objective_function = 0.5*rhs.transpose()*rhs;

  return objective_function;
}


//see Eq. 33, Paper III
double CondPhaseSolver::assembleRightHandSideFull(
  const std::vector<Condensate*>& condensates,
  const std::vector<double>& activity_corr,
  const std::vector<double>& number_densities,
  const std::vector<double>& log_number_densities,
  const std::vector<double>& log_activity_corr,
  const std::vector< Element* >& elements,
  const std::vector< Molecule >& molecules,
  const double total_element_density,
  const Eigen::VectorXdt& scaling_factors,
  Eigen::VectorXdt& rhs)
{
  const size_t nb_elements = elements.size();
  const size_t nb_cond = condensates.size();

  rhs.setZero(nb_elements + 2*nb_cond);

  for (size_t i=0; i<nb_cond; ++i)
  {
    rhs(i) = condensates[i]->log_tau - log_activity_corr[i] - log_number_densities[i];
    rhs(i+nb_cond) = - condensates[i]->log_activity - activity_corr[i];
  }


  for (size_t i=0; i<nb_elements; ++i)
  { 
    rhs(i+2*nb_cond) = total_element_density * elements[i]->epsilon - elements[i]->number_density;
    
    for (auto j : elements[i]->molecule_list)
      rhs(i+2*nb_cond) -= molecules[j].stoichiometric_vector[elements[i]->index] * molecules[j].number_density;
    
    for (size_t j=0; j<condensates.size(); ++j)
      rhs(i+2*nb_cond) -= condensates[j]->stoichiometric_vector[elements[i]->index] * number_densities[j];
  }


  for (auto i=0; i<rhs.rows(); ++i)
    rhs(i) /= scaling_factors(i);


  const double objective_function = 0.5*rhs.transpose()*rhs;

  return objective_function;
}



Eigen::MatrixXdt CondPhaseSolver::assemblePerturbedHessian(
  const Eigen::MatrixXdt& jacobian,
  const double perturbation)
{
  Eigen::MatrixXdt hessian = jacobian.transpose()*jacobian;

  const double norm = hessian.template lpNorm<1>();

  for (auto i=0; i<hessian.rows(); ++i)
    hessian(i,i) += perturbation * norm; 

  return hessian;
}




bool CondPhaseSolver::solveSystem(
  const Eigen::MatrixXdt& jacobian,
  const Eigen::VectorXdt& rhs,
  Eigen::VectorXdt& result)
{
  if (options.cond_use_lm)
  {
    bool needs_clamping = false;

    for (auto i = 0; i < jacobian.rows(); ++i)
    {
      if (std::abs(jacobian(i,i)) < lm_mu)
      {
        needs_clamping = true;
        break;
      }
    }

    if (!needs_clamping)
    {
      Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> solver;
      solver.compute(jacobian);
      result = solver.solve(rhs);
      return true;
    }

    if (options.verbose_level >= 3)
      std::cout << "FastChem: Jacobian near-singular, clamping diagonal (mu=" << lm_mu << ")\n";

    Eigen::MatrixXdt J_reg = jacobian;

    for (auto i = 0; i < J_reg.rows(); ++i)
    {
      if (std::abs(J_reg(i,i)) < lm_mu)
        J_reg(i,i) = (J_reg(i,i) >= 0) ? lm_mu : -lm_mu;
    }

    Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> lm_solver;
    lm_solver.compute(J_reg);
    result = lm_solver.solve(rhs);

    // If the clamped system still produced NaN/inf, fall back to the minimum-norm
    // SVD solution, which is always finite and bounded for any near-singular matrix.
    bool result_bad = false;
    for (int i = 0; i < result.rows(); ++i)
    {
      if (!std::isfinite(static_cast<double>(result(i))))
      {
        result_bad = true;
        break;
      }
    }

    if (result_bad)
    {
      if (options.verbose_level >= 3)
        std::cout << "FastChem: PartialPivLU on clamped Jacobian produced NaN/inf, switching to SVD\n";

      result = jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
    }

    return false;
  }

  if (!options.cond_use_full_pivot)
  {
    Eigen::PartialPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> solver;
    solver.compute(jacobian);
    result = solver.solve(rhs);

    return true;
  }

  Eigen::FullPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> solver;
  solver.compute(jacobian);
  result = solver.solve(rhs);

  if (!solver.isInvertible())
  {
    if (options.verbose_level >= 3)
      std::cout << "FastChem warning: Jacobian is (almost) singular! ";

    if (options.cond_use_svd)
    {
      if (options.verbose_level >= 3)
        std::cout << "Switching to Singular Value Decomposition.\n";

      Eigen::BDCSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> solver;
      result = jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
    }
    else
    {
      if (options.verbose_level >= 3)
        std::cout << "Switching to diagonal clamping fallback.\n";

      Eigen::MatrixXdt J_reg = jacobian;
      const double min_diag = 0.01;

      for (auto i = 0; i < J_reg.rows(); ++i)
      {
        if (std::abs(J_reg(i,i)) < min_diag)
          J_reg(i,i) = (J_reg(i,i) >= 0) ? min_diag : -min_diag;
      }

      solver.compute(J_reg);
      result = solver.solve(rhs);
    }

    return false;
  }

  return true;
}



void CondPhaseSolver::resetLM()
{
  lm_mu = 0.01;
  lm_objective_prev = 0;
  lm_has_prev_objective = false;
}


void CondPhaseSolver::adaptLM(const double objective_function)
{
  if (lm_has_prev_objective)
  {
    if (objective_function < lm_objective_prev)
      lm_mu *= 0.5;
    else
      lm_mu *= 2.0;

    if (lm_mu < 1e-4) lm_mu = 1e-4;
    if (lm_mu > 0.5) lm_mu = 0.5;
  }

  lm_objective_prev = objective_function;
  lm_has_prev_objective = true;
}



double CondPhaseSolver::backtrackStep(
  const double objective_function_0,
  const double objective_function_prev,
  const double objective_function_2prev,
  const double lambda_prev,
  const double lambda_2prev)
{
  const double object_function_deriv = -2.0*objective_function_0;

  double lambda = 0;

  if (lambda_2prev == 0)
  {
    lambda = -object_function_deriv /(2.0 * (objective_function_prev - objective_function_0 - object_function_deriv));
  }
  else
  {
    Eigen::Matrix<double, 2, 2> cubic_fit_matrix;
    Eigen::Matrix<double, 2,1> cubic_fit_vector;

    cubic_fit_matrix(0,0) = 1.0/(lambda_prev*lambda_prev);
    cubic_fit_matrix(1,0) = -lambda_2prev/(lambda_prev*lambda_prev);
    cubic_fit_matrix(0,1) = -1.0/(lambda_2prev*lambda_2prev);
    cubic_fit_matrix(1,1) = lambda_prev/(lambda_2prev*lambda_2prev);

    cubic_fit_vector(0) = objective_function_prev - objective_function_0 - object_function_deriv*lambda_prev;
    cubic_fit_vector(1) = objective_function_2prev - objective_function_0 - object_function_deriv*lambda_2prev;

    cubic_fit_vector /= lambda_prev-lambda_2prev;

    Eigen::Matrix<double, 2,1> cubic_fit = cubic_fit_matrix*cubic_fit_vector;
    const double a = cubic_fit(0);
    const double b = cubic_fit(1);

    lambda = (-b + std::sqrt(b*b - 3*a*object_function_deriv))/(3*a);
  }

  if (lambda < 0.1*lambda_prev) lambda = 0.1*lambda_prev;
  if (lambda > 0.5*lambda_prev) lambda = 0.5*lambda_prev;
  
  return lambda;
}



double CondPhaseSolver::objectiveFunction(
  const std::vector<Condensate*>& condensates,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<double>& activity_corr,
  const std::vector<double>& number_denities,
  const std::vector<double>& log_number_densities,
  const std::vector<double>& log_activity_corr,
  const std::vector< Element* >& elements,
  const std::vector< Molecule >& molecules,
  const double total_element_density,
  const Eigen::VectorXdt& scaling_factors)
{
  double objective_function = 0;

  if (options.cond_reduce_system_size == false)
  {
    Eigen::VectorXdt rhs;
    objective_function = assembleRightHandSideFull(
      condensates,
      activity_corr,
      number_denities,
      log_number_densities,
      log_activity_corr,
      elements,
      molecules,
      total_element_density,
      scaling_factors,
      rhs);
  }
  else
  {
    Eigen::VectorXdt rhs;
    objective_function = assembleRightHandSide(
      condensates,
      condensates_jac,
      condensates_rem,
      activity_corr,
      number_denities,
      log_number_densities,
      log_activity_corr,
      elements,
      molecules,
      total_element_density,
      scaling_factors,
      rhs);
  }

  return objective_function;
}


}



