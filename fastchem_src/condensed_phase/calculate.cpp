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


#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "condensed_phase.h"

#include "../fastchem_constants.h"
#include "../species_struct.h"
#include "solver.h"
#include "../../_ext/Eigen/Dense"


namespace fastchem {


//This is the main FastChem iteration for the condensed phase
template <class double_type>
bool CondensedPhase<double_type>::calculate(
  std::vector<Condensate<double_type>*>& condensates_act,
  std::vector<Element<double_type>*>& elements_cond,
  const double temperature,
  const double density,
  const double total_element_density,
  std::vector<Molecule<double_type>>& molecules,
  unsigned int& nb_iterations)
{ 
  if (condensates_act.size() == 0) return true;


  double_type tau = options.cond_tau;

  for (auto & i : condensates_act)
  {
    i->tau = tau * element_data.elements[i->reference_element].epsilon * total_element_density;
    i->log_tau = std::log(i->tau);
  }


  for (auto & i : condensates_act)
  {
    if (i->number_density == 0 && i->activity_correction == 0)
    {
       i->number_density = i->max_number_density;
       i->activity_correction = 1.0;
    }
  }

  std::vector<unsigned int> condensates_jac;
  condensates_jac.reserve(condensates_act.size());

  std::vector<unsigned int> condensates_rem;
  condensates_rem.reserve(condensates_act.size());


  std::vector<double_type> elem_densities_old(elements_cond.size(), 0.0);
  std::vector<double_type> elem_densities_new(elements_cond.size(), 0.0);

  std::vector<double_type> cond_densities_old(condensates_act.size(), 0.0);
  std::vector<double_type> cond_densities_new(condensates_act.size(), 0.0);

  std::vector<double_type> activity_corr_old(condensates_act.size(), 1.0);
  std::vector<double_type> activity_corr_new(condensates_act.size(), 0.0);

  for (size_t i=0; i<elements_cond.size(); ++i)
    elem_densities_old[i] = elements_cond[i]->number_density;


  for (size_t i=0; i<condensates_act.size(); ++i)
  { 
    if (condensates_act[i]->number_density == 0.0) condensates_act[i]->number_density = 1e-20;

    cond_densities_old[i] = condensates_act[i]->number_density;
    activity_corr_old[i] = condensates_act[i]->activity_correction;
  }

  bool cond_converged = false;

  double_type limit = options.cond_iter_change_limit;


  for (nb_iterations=0; nb_iterations<options.nb_max_cond_iter; ++nb_iterations)
  { 
    double_type objective_function_0 = 0;
    Eigen::VectorXdt<double_type> scaling_factors;
    Eigen::VectorXdt<double_type> result;


    double_type max_delta = 0;
    bool system_invertible = true;

    if (options.cond_solve_full_matrix)
    {
      system_invertible = solver.newtonStepFull(
        condensates_act,
        elements_cond,
        molecules,
        total_element_density,
        cond_densities_old,
        elem_densities_old,
        activity_corr_old,
        result,
        scaling_factors,
        objective_function_0);

      Eigen::VectorXdt<double_type> result_scaled = result;

      double_type max_value = result_scaled.cwiseAbs().maxCoeff();

      if (max_value > limit)
        result_scaled *= limit/max_value;

      max_delta = correctValuesFull(
        result_scaled,
        condensates_act,
        activity_corr_old,
        activity_corr_new,
        cond_densities_old,
        cond_densities_new,
        elements_cond,
        elem_densities_old,
        elem_densities_new);
    }
    else
    {
      selectJacobianCondensates(
        condensates_act,
        cond_densities_old,
        activity_corr_old,
        condensates_jac,
        condensates_rem);

      system_invertible = solver.newtonStep(
        condensates_act,
        condensates_jac,
        condensates_rem,
        elements_cond,
        molecules,
        total_element_density,
        cond_densities_old,
        elem_densities_old,
        activity_corr_old,
        result,
        scaling_factors,
        objective_function_0);

      Eigen::VectorXdt<double_type> result_scaled = result;

      double_type max_value = result_scaled.cwiseAbs().maxCoeff();

      if (max_value > limit)
        result_scaled *= limit/max_value;

      max_delta = correctValues(
        result_scaled,
        condensates_act,
        condensates_jac,
        condensates_rem,
        activity_corr_old,
        activity_corr_new,
        cond_densities_old,
        cond_densities_new,
        elements_cond,
        elem_densities_old,
        elem_densities_new,
        limit);
    }


    for (size_t i=0; i<elements_cond.size(); ++i)
      elements_cond[i]->number_density = elem_densities_new[i];

    for (auto & i : condensates_act)
      i->calcActivity(temperature, elements, options.cond_use_data_validity_limits);

    for (auto & i : molecules)
      i.calcNumberDensity(elements);

    // double_type max_delta1 = newtonBacktrack(
    //   objective_function_0,
    //   result,
    //   scaling_factors,
    //   condensates_act,
    //   condensates_jac,
    //   condensates_rem,
    //   activity_corr_old,
    //   activity_corr_new,
    //   cond_densities_old,
    //   cond_densities_new,
    //   elements_cond,
    //   elem_densities_old,
    //   elem_densities_new,
    //   molecules,
    //   total_element_density,
    //   temperature,
    //   limit);

    {
    //std::cout << "iter: " << nb_iterations << "  " << objective_function_0 << "\n";
    //for (size_t i=0; i<condensates_act.size(); ++i)
      //std::cout <<i << "  " << condensates_act[i]->symbol << "\t" << cond_densities_old[i] << "\t" << cond_densities_new[i] << "\t" << activity_corr_old[i] << "\t" << activity_corr_new[i] << "\t" << condensates_act[i]->log_activity << "\t" << condensates_act[i]->tau << "\n";
    }


    elem_densities_old = elem_densities_new;
    cond_densities_old = cond_densities_new;
    activity_corr_old = activity_corr_new;

    cond_converged = max_delta < options.cond_accuracy && system_invertible && objective_function_0 < 0.001;
  
    if (cond_converged) break;
  }


  for (size_t i=0; i<condensates_act.size(); ++i)
  {
    condensates_act[i]->number_density = cond_densities_new[i];
    condensates_act[i]->activity_correction = activity_corr_new[i];
  }


  double_type phi_sum = 0;

  for (auto & i : elements)
  {
    i.calcDegreeOfCondensation(condensates, total_element_density);
    phi_sum += i.phi;
  }

  for (auto & i : elements)
    i.normalisePhi(phi_sum);

  return cond_converged;
}



template <class double_type>
double_type CondensedPhase<double_type>::correctValues(
  const Eigen::VectorXdt<double_type>& result,
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<double_type>& activity_corr_old,
  std::vector<double_type>& activity_corr_new,
  const std::vector<double_type>& cond_number_dens_old,
  std::vector<double_type>& cond_number_dens_new,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<double_type>& elem_number_dens_old,
  std::vector<double_type>& elem_number_dens_new,
  const double max_change)
{
  std::vector<double_type> delta_n_cond(condensates.size(), 0);
  const size_t nb_cond_jac = condensates_jac.size();

  for (size_t i=0; i<nb_cond_jac; ++i)
    delta_n_cond[condensates_jac[i]] = result(i);

  for (size_t i=0; i<condensates_rem.size(); ++i)
  {
    const unsigned int index = condensates_rem[i];

    delta_n_cond[index] = condensates[index]->log_activity/activity_corr_old[index] 
      + condensates[index]->log_tau 
      - std::log(cond_number_dens_old[index]) 
      - std::log(activity_corr_old[index]) + 1;

    for (size_t j=0; j<elements.size(); ++j)
      delta_n_cond[index] += condensates[index]->stoichiometric_vector[elements[j]->index] 
        * result(j+nb_cond_jac) 
        / activity_corr_old[index];
  }


  double_type max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {
    if (delta_n_cond[i] > max_change) 
      delta_n_cond[i] = max_change;

    if (delta_n_cond[i] < -max_change) 
      delta_n_cond[i] = -max_change;

    cond_number_dens_new[i] = cond_number_dens_old[i] * std::exp(delta_n_cond[i]);

    if (cond_number_dens_new[i] > condensates[i]->max_number_density) 
      cond_number_dens_new[i] = condensates[i]->max_number_density;


    double_type delta_lambda = condensates[i]->log_tau 
      - std::log(activity_corr_old[i]) 
      - std::log(cond_number_dens_old[i]) 
      - delta_n_cond[i];
   
    if (delta_lambda > max_change) 
      delta_lambda = max_change;
    
    if (delta_lambda < -max_change) 
      delta_lambda = -max_change;

    activity_corr_new[i] = activity_corr_old[i] * std::exp(delta_lambda);

    double_type delta = std::fabs(cond_number_dens_new[i] - cond_number_dens_old[i])/cond_number_dens_old[i]; 
    if (delta > max_delta) max_delta = delta;

    delta = std::fabs(activity_corr_new[i] - activity_corr_old[i])/activity_corr_old[i];
    if (delta > max_delta) max_delta = delta;
  }


  for (size_t i=0; i<elements.size(); ++i)
  { 
    double_type delta_n_elem = result(i + nb_cond_jac);

    if (delta_n_elem > max_change) 
      delta_n_elem = max_change;
    
    if (delta_n_elem < -max_change) 
      delta_n_elem = -max_change;

    elem_number_dens_new[i] = elem_number_dens_old[i] * std::exp(delta_n_elem);

    const double_type delta = std::fabs(elem_number_dens_new[i] - elem_number_dens_old[i])/elem_number_dens_old[i]; 
    if (delta > max_delta) max_delta = delta;
  }

  return max_delta;
}



template <class double_type>
double_type CondensedPhase<double_type>::correctValuesFull(
  const Eigen::VectorXdt<double_type>& result,
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr_old,
  std::vector<double_type>& activity_corr_new,
  const std::vector<double_type>& cond_number_dens_old,
  std::vector<double_type>& cond_number_dens_new,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<double_type>& elem_number_dens_old,
  std::vector<double_type>& elem_number_dens_new)
{ 
  double_type max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {  
    const double_type delta_n = result(i);
    const double_type delta_lambda = result(i+condensates.size());

    cond_number_dens_new[i] = cond_number_dens_old[i] * std::exp(delta_n);
    if (cond_number_dens_new[i] > condensates[i]->max_number_density) cond_number_dens_new[i] = condensates[i]->max_number_density;

    activity_corr_new[i] = activity_corr_old[i] * std::exp(delta_lambda);

    double_type delta = std::fabs(cond_number_dens_new[i] - cond_number_dens_old[i])/cond_number_dens_old[i]; 
    if (delta > max_delta) max_delta = delta;

    delta = std::fabs(activity_corr_new[i] - activity_corr_old[i])/activity_corr_old[i];
    if (delta > max_delta) max_delta = delta;
  }


  for (size_t i=0; i<elements.size(); ++i)
  { 
    const double_type delta_n = result(i + 2*condensates.size());

    elem_number_dens_new[i] = elem_number_dens_old[i] * std::exp(delta_n);

    const double_type delta = std::fabs(elem_number_dens_new[i] - elem_number_dens_old[i])/elem_number_dens_old[i]; 
    if (delta > max_delta) max_delta = delta;
  }

  return max_delta;
}



template <class double_type>
double_type CondensedPhase<double_type>::newtonBacktrack(
  const double_type objective_function_0,
  const Eigen::VectorXdt<double_type>& result,
  const Eigen::VectorXdt<double_type>& scaling_factors,
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<double_type>& activity_corr_old,
  std::vector<double_type>& activity_corr_new,
  const std::vector<double_type>& cond_number_dens_old,
  std::vector<double_type>& cond_number_dens_new,
  const std::vector<Element<double_type>*>& elements_cond,
  const std::vector<double_type>& elem_number_dens_old,
  std::vector<double_type>& elem_number_dens_new,
  std::vector<Molecule<double_type>>& molecules,
  const double total_element_density,
  const double temperature,
  const double max_change)
{ 
  double_type objective_function_1 = solver.objectiveFunction(
    condensates,
    condensates_jac,
    condensates_rem,
    activity_corr_new,
    cond_number_dens_new,
    elements_cond,
    molecules,
    total_element_density,
    scaling_factors);

  double_type max_delta = 0;

  double_type lambda_prev = 1.0;
  double_type lambda_2prev = 0.0;

  double_type lamdba_min = 1e-5;

  if (objective_function_1 > objective_function_0)
  {
      const double_type object_function_deriv = -2.0*objective_function_0;

      double_type objective_function_prev = objective_function_1;
      double_type objective_function_2prev = 0.0;

      double_type lambda = 1;
      double_type limit1 = max_change;

      while(objective_function_prev > objective_function_0 + 1e-4*lambda*object_function_deriv && lambda_prev > lamdba_min)
      { 
        lambda = solver.backtrackStep(
          objective_function_0,
          objective_function_prev,
          objective_function_2prev,
          lambda_prev,
          lambda_2prev);

        limit1 = max_change * lambda;

        Eigen::VectorXdt<double_type> result_scaled = result;

        double_type max_value = result_scaled.cwiseAbs().maxCoeff();

        if (max_value > max_change)
          result_scaled *= limit1/max_value;

        if (options.cond_solve_full_matrix)
        { 
          max_delta = correctValuesFull(
            result_scaled,
            condensates,
            activity_corr_old,
            activity_corr_new,
            cond_number_dens_old,
            cond_number_dens_new,
            elements_cond,
            elem_number_dens_old,
            elem_number_dens_new);
        }
        else
        {
          max_delta = correctValues(
            result_scaled,
            condensates,
            condensates_jac,
            condensates_rem,
            activity_corr_old,
            activity_corr_new,
            cond_number_dens_old,
            cond_number_dens_new,
            elements_cond,
            elem_number_dens_old,
            elem_number_dens_new,
            limit1);
        }

        for (size_t i=0; i<elements_cond.size(); ++i)
          elements_cond[i]->number_density = elem_number_dens_new[i];

        for (auto & i : condensates)
          i->calcActivity(temperature, elements, options.cond_use_data_validity_limits);

        for (auto & i : molecules)
          i.calcNumberDensity(elements);

        objective_function_2prev = objective_function_prev;
        lambda_2prev = lambda_prev;
        lambda_prev = lambda;

        objective_function_prev = solver.objectiveFunction(
          condensates,
          condensates_jac,
          condensates_rem,
          activity_corr_new,
          cond_number_dens_new,
          elements_cond,
          molecules,
          total_element_density,
          scaling_factors);
      }

      if (lambda_prev <= lamdba_min)
      { 
        Eigen::VectorXdt<double_type> result_scaled = result;

        double_type max_value = result_scaled.cwiseAbs().maxCoeff();

        if (max_value > 1.0)
          result_scaled *= 1.0/max_value;
        
        if (options.cond_solve_full_matrix)
        {
          max_delta = correctValuesFull(
            result_scaled,
            condensates,
            activity_corr_old,
            activity_corr_new,
            cond_number_dens_old,
            cond_number_dens_new,
            elements_cond,
            elem_number_dens_old,
            elem_number_dens_new);
        }
        else
        {
          max_delta = correctValues(
            result_scaled,
            condensates,
            condensates_jac,
            condensates_rem,
            activity_corr_old,
            activity_corr_new,
            cond_number_dens_old,
            cond_number_dens_new,
            elements_cond,
            elem_number_dens_old,
            elem_number_dens_new,
            1.0);
        }

        for (size_t i=0; i<elements_cond.size(); ++i)
          elements_cond[i]->number_density = elem_number_dens_new[i];

        for (auto & i : condensates)
          i->calcActivity(temperature, elements, options.cond_use_data_validity_limits);

        for (auto & i : molecules)
          i.calcNumberDensity(elements);
      }

    }

  return max_delta;
}


template class CondensedPhase<double>;
template class CondensedPhase<long double>;
}