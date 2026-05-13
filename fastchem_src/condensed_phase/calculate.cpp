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
bool CondensedPhase::calculate(
  std::vector<Condensate*>& condensates_act,
  std::vector<Element*>& elements_cond,
  const double temperature,
  const double density,
  const double total_element_density,
  std::vector<Molecule>& molecules,
  unsigned int& nb_iterations)
{
  if (condensates_act.size() == 0) return true;

  double N_total = total_element_density;

  double tau = options.cond_tau;

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


  //log-space primary iteration variables
  std::vector<double> log_elem_dens_old(elements_cond.size(), static_cast<double>(LOG_DENSITY_FLOOR));
  std::vector<double> log_elem_dens_new(elements_cond.size(), static_cast<double>(LOG_DENSITY_FLOOR));

  std::vector<double> log_cond_dens_old(condensates_act.size(), static_cast<double>(LOG_DENSITY_FLOOR));
  std::vector<double> log_cond_dens_new(condensates_act.size(), static_cast<double>(LOG_DENSITY_FLOOR));

  std::vector<double> log_activity_corr_old(condensates_act.size(), 0.0);
  std::vector<double> log_activity_corr_new(condensates_act.size(), static_cast<double>(LOG_DENSITY_FLOOR));

  for (size_t i=0; i<elements_cond.size(); ++i)
    log_elem_dens_old[i] = elements_cond[i]->log_number_density;


  for (size_t i=0; i<condensates_act.size(); ++i)
  {
    if (condensates_act[i]->number_density == 0.0) condensates_act[i]->number_density = 1e-20;

    log_cond_dens_old[i] = std::log(condensates_act[i]->number_density);
    log_activity_corr_old[i] = safeLog(condensates_act[i]->activity_correction);
  }

  bool cond_converged = false;

  double limit = options.cond_iter_change_limit;

  if (options.cond_use_lm)
    solver.resetLM();
  
  for (nb_iterations=0; nb_iterations<options.nb_max_cond_iter; ++nb_iterations)
  {
    double objective_function_0 = 0;
    Eigen::VectorXdt scaling_factors;
    Eigen::VectorXdt result;


    double max_delta = 0;
    bool system_invertible = true;

    //compute temporary linear vectors from log for Jacobian and Newton step
    std::vector<double> cond_dens_linear(condensates_act.size());
    std::vector<double> activity_corr_linear(condensates_act.size());
    std::vector<double> elem_dens_linear(elements_cond.size());

    for (size_t i=0; i<condensates_act.size(); ++i)
    {
      cond_dens_linear[i] = safeExp(log_cond_dens_old[i]);
      activity_corr_linear[i] = safeExp(log_activity_corr_old[i]);
    }

    for (size_t i=0; i<elements_cond.size(); ++i)
      elem_dens_linear[i] = safeExp(log_elem_dens_old[i]);

    if (options.cond_reduce_system_size == false)
    {
      system_invertible = solver.newtonStepFull(
        condensates_act,
        elements_cond,
        molecules,
        N_total,
        cond_dens_linear,
        elem_dens_linear,
        activity_corr_linear,
        log_cond_dens_old,
        log_activity_corr_old,
        result,
        scaling_factors,
        objective_function_0);

      Eigen::VectorXdt result_scaled = result;

      double max_value = result_scaled.cwiseAbs().maxCoeff();

      if (max_value > limit)
        result_scaled *= limit/max_value;

      max_delta = correctValuesFull(
        result_scaled,
        condensates_act,
        log_activity_corr_old,
        log_activity_corr_new,
        log_cond_dens_old,
        log_cond_dens_new,
        elements_cond,
        log_elem_dens_old,
        log_elem_dens_new);
    }
    else
    {
      selectJacobianCondensates(
        condensates_act,
        cond_dens_linear,
        activity_corr_linear,
        condensates_jac,
        condensates_rem);

      system_invertible = solver.newtonStep(
        condensates_act,
        condensates_jac,
        condensates_rem,
        elements_cond,
        molecules,
        N_total,
        cond_dens_linear,
        elem_dens_linear,
        activity_corr_linear,
        log_cond_dens_old,
        log_activity_corr_old,
        result,
        scaling_factors,
        objective_function_0);

      Eigen::VectorXdt result_scaled = result;

      double max_value = result_scaled.cwiseAbs().maxCoeff();

      if (max_value > limit)
        result_scaled *= limit/max_value;

      max_delta = correctValues(
        result_scaled,
        condensates_act,
        condensates_jac,
        condensates_rem,
        log_activity_corr_old,
        log_activity_corr_new,
        log_cond_dens_old,
        log_cond_dens_new,
        elements_cond,
        log_elem_dens_old,
        log_elem_dens_new,
        limit);
    }


    for (size_t i=0; i<elements_cond.size(); ++i)
    {
      elements_cond[i]->log_number_density = log_elem_dens_new[i];
      elements_cond[i]->number_density = safeExp(log_elem_dens_new[i]);
    }

    for (auto & i : condensates_act)
      i->calcActivity(temperature, elements, options.cond_use_data_validity_limits);

    for (auto & i : molecules)
      i.calcNumberDensity(elements);

    log_elem_dens_old = log_elem_dens_new;
    log_cond_dens_old = log_cond_dens_new;
    log_activity_corr_old = log_activity_corr_new;

    if (options.cond_use_lm)
      solver.adaptLM(objective_function_0);

    cond_converged = max_delta < options.cond_accuracy
      && (system_invertible || options.cond_use_lm)
      && objective_function_0 < 0.001;
    
    if (cond_converged) break;
  }


  for (size_t i=0; i<condensates_act.size(); ++i)
  {
    condensates_act[i]->number_density = safeExp(log_cond_dens_new[i]);
    condensates_act[i]->activity_correction = safeExp(log_activity_corr_new[i]);
  }

  return cond_converged;
}



double CondensedPhase::correctValues(
  const Eigen::VectorXdt& result,
  const std::vector<Condensate*>& condensates,
  const std::vector<unsigned int>& condensates_jac,
  const std::vector<unsigned int>& condensates_rem,
  const std::vector<double>& log_activity_corr_old,
  std::vector<double>& log_activity_corr_new,
  const std::vector<double>& log_cond_dens_old,
  std::vector<double>& log_cond_dens_new,
  const std::vector<Element*>& elements,
  const std::vector<double>& log_elem_dens_old,
  std::vector<double>& log_elem_dens_new,
  const double max_change)
{
  std::vector<double> delta_n_cond(condensates.size(), 0);
  const size_t nb_cond_jac = condensates_jac.size();

  for (size_t i=0; i<nb_cond_jac; ++i)
    delta_n_cond[condensates_jac[i]] = result(i);

  for (size_t i=0; i<condensates_rem.size(); ++i)
  {
    const unsigned int index = condensates_rem[i];

    const double activity_corr_lin = safeExp(log_activity_corr_old[index]);

    if (activity_corr_lin == 0.0) continue;

    delta_n_cond[index] = condensates[index]->log_activity/activity_corr_lin
      + condensates[index]->log_tau
      - log_cond_dens_old[index]
      - log_activity_corr_old[index] + 1;

    for (size_t j=0; j<elements.size(); ++j)
      delta_n_cond[index] += condensates[index]->stoichiometric_vector[elements[j]->index]
        * result(j+nb_cond_jac)
        / activity_corr_lin;
  }


  double max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {
    if (delta_n_cond[i] > max_change)
      delta_n_cond[i] = max_change;

    if (delta_n_cond[i] < -max_change)
      delta_n_cond[i] = -max_change;

    //log-space additive update (Eq. 39)
    log_cond_dens_new[i] = log_cond_dens_old[i] + delta_n_cond[i];

    const double log_max = std::log(condensates[i]->max_number_density);
    if (log_cond_dens_new[i] > log_max)
      log_cond_dens_new[i] = log_max;


    double delta_lambda = condensates[i]->log_tau
      - log_activity_corr_old[i]
      - log_cond_dens_old[i]
      - delta_n_cond[i];

    if (delta_lambda > max_change)
      delta_lambda = max_change;

    if (delta_lambda < -max_change)
      delta_lambda = -max_change;

    log_activity_corr_new[i] = log_activity_corr_old[i] + delta_lambda;

    //convergence metric: |delta| in log-space
    double delta_conv = std::fabs(delta_n_cond[i]);
    if (delta_conv > max_delta) max_delta = delta_conv;

    delta_conv = std::fabs(delta_lambda);
    if (delta_conv > max_delta) max_delta = delta_conv;
  }


  for (size_t i=0; i<elements.size(); ++i)
  {
    double delta_n_elem = result(i + nb_cond_jac);

    if (delta_n_elem > max_change)
      delta_n_elem = max_change;

    if (delta_n_elem < -max_change)
      delta_n_elem = -max_change;

    log_elem_dens_new[i] = log_elem_dens_old[i] + delta_n_elem;

    const double delta_conv = std::fabs(delta_n_elem);
    if (delta_conv > max_delta) max_delta = delta_conv;
  }

  return max_delta;
}



double CondensedPhase::correctValuesFull(
  const Eigen::VectorXdt& result,
  const std::vector<Condensate*>& condensates,
  const std::vector<double>& log_activity_corr_old,
  std::vector<double>& log_activity_corr_new,
  const std::vector<double>& log_cond_dens_old,
  std::vector<double>& log_cond_dens_new,
  const std::vector<Element*>& elements,
  const std::vector<double>& log_elem_dens_old,
  std::vector<double>& log_elem_dens_new)
{
  double max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {
    const double delta_n = result(i);
    const double delta_lambda = result(i+condensates.size());

    //log-space additive update (Eq. 39)
    log_cond_dens_new[i] = log_cond_dens_old[i] + delta_n;

    const double log_max = std::log(condensates[i]->max_number_density);
    if (log_cond_dens_new[i] > log_max)
      log_cond_dens_new[i] = log_max;

    log_activity_corr_new[i] = log_activity_corr_old[i] + delta_lambda;

    //convergence metric: |delta| in log-space
    double delta_conv = std::fabs(delta_n);
    if (delta_conv > max_delta) max_delta = delta_conv;

    delta_conv = std::fabs(delta_lambda);
    if (delta_conv > max_delta) max_delta = delta_conv;
  }


  for (size_t i=0; i<elements.size(); ++i)
  {
    const double delta_n = result(i + 2*condensates.size());

    log_elem_dens_new[i] = log_elem_dens_old[i] + delta_n;

    const double delta_conv = std::fabs(delta_n);
    if (delta_conv > max_delta) max_delta = delta_conv;
  }

  return max_delta;
}



}
