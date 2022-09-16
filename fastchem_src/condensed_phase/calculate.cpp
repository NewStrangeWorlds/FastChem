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


//This is the main FastChem iteration for the gas phase
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

  double_type tau = 1e-25;
  double_type log_tau = std::log(tau);

  for (auto & i : condensates_act)
  {
    i->tau = tau * element_data.elements[i->reference_element].epsilon * total_element_density;
    i->log_tau = std::log(i->tau);
    //i->log_tau = std::log(tau);
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

    std::cout << condensates_act[i]->symbol << "\t" << condensates_act[i]->log_activity << "\n";
  }
  

  bool cond_converged = false;

  for (nb_iterations=0; nb_iterations<1000; ++nb_iterations)
  {
    selectJacobianCondensates(
      condensates_act,
      cond_densities_old,
      activity_corr_old,
      condensates_jac,
      condensates_rem);

    Eigen::MatrixXdt<double_type> jacobian = solver.assembleJacobian(
      condensates_act,
      activity_corr_old,
      cond_densities_old,
      condensates_jac,
      condensates_rem,
      elements_cond,
      molecules);

    Eigen::VectorXdt<double_type> rhs = solver.assembleRightHandSide(
      condensates_act,
      condensates_jac,
      condensates_rem,
      activity_corr_old,
      cond_densities_old,
      elements_cond,
      molecules,
      total_element_density,
      log_tau);

    /*Eigen::VectorXdt<double_type> rhs = solver.assembleRightHandSide(
      condensates_act,
      condensates_jac,
      condensates_rem,
      activity_corr_old,
      cond_densities_old,
      elements_cond,
      molecules,
      total_element_density);*/


    std::vector<double_type> result = solver.solveSystem(jacobian, rhs);


    double_type max_delta = correctValues(
      result,
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
      log_tau,
      1.0);

    /*double_type max_delta = correctValues(
      result,
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
      1.0);*/

    for (size_t i=0; i<elements_cond.size(); ++i)
      elements_cond[i]->number_density = elem_densities_new[i];

    for (auto & i : condensates_act)  i->calcActivity(temperature, elements);

    for (auto & i : molecules)  i.calcNumberDensity(elements);
    
    //if (nb_iterations > 8990)
    {
    std::cout << "iter: " << nb_iterations << "\n";
    for (size_t i=0; i<condensates_act.size(); ++i)
      std::cout << condensates_act[i]->symbol << "\t" << cond_densities_old[i] << "\t" << cond_densities_new[i] << "\t" << activity_corr_old[i] << "\t" << activity_corr_new[i] << "\t" << condensates_act[i]->log_activity << "\n";
    }

    elem_densities_old = elem_densities_new;
    cond_densities_old = cond_densities_new;
    activity_corr_old = activity_corr_new;
    
    std::cout << "cond delta " << max_delta << "\n";
    cond_converged = max_delta < 1e-6;

    if (cond_converged) break;
  }


  for (size_t i=0; i<condensates_act.size(); ++i)
  {
    condensates_act[i]->number_density = cond_densities_new[i];
    condensates_act[i]->activity_correction = activity_corr_new[i];

    //if (condensates_act[i]->number_density <= 1e-15) condensates_act[i]->number_density = 0.0;
    if (condensates_act[i]->log_activity < -1e-1) condensates_act[i]->number_density = 0.0;
  }

  for (size_t i=0; i<elements_cond.size(); ++i)
  {
    elements_cond[i]->number_density = elem_densities_new[i];
  }

  for (auto & i : elements)
  {
    i.calcDegreeOfCondensation(condensates, total_element_density);
  }


  double_type phi_sum = 0;

  for (auto & i : elements)
    phi_sum += i.phi;

  for (auto & i : elements)
    i.normalisePhi(phi_sum);

  return cond_converged;
}



template <class double_type>
double_type CondensedPhase<double_type>::correctValues(
  const std::vector<double_type>& result,
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
  const double_type ln_tau,
  const double max_change)
{
  std::vector<double_type> delta_n(condensates.size(), 0);

  for (size_t i=0; i<condensates_jac.size(); ++i)
    delta_n[condensates_jac[i]] = result[i];


  const size_t nb_cond_jac = condensates_jac.size();

  for (size_t i=0; i<condensates_rem.size(); ++i)
  {
    const unsigned int index = condensates_rem[i];

    for (size_t j=0; j<elements.size(); ++j)
      delta_n[index] += condensates[index]->stoichiometric_vector[elements[j]->index] * result[j+nb_cond_jac];

    delta_n[index] /= activity_corr_old[index];

    /*delta_n[index] += condensates[index]->log_activity / activity_corr_old[index] 
                    + ln_tau - std::log(activity_corr_old[index]) 
                    - std::log(cond_number_dens_old[index]) 
                    + 1;*/
    delta_n[index] += condensates[index]->log_activity / activity_corr_old[index] 
                    + condensates[index]->log_tau - std::log(activity_corr_old[index]) 
                    - std::log(cond_number_dens_old[index]) 
                    + 1;
  }


  double_type max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {
    if (delta_n[i] > max_change) delta_n[i] = max_change;
    if (delta_n[i] < -max_change) delta_n[i] = -max_change;

    if (std::fabs(delta_n[i]) > max_delta) max_delta = std::fabs(delta_n[i]);


    cond_number_dens_new[i] = cond_number_dens_old[i] * std::exp(delta_n[i]);
    if (cond_number_dens_new[i] > condensates[i]->max_number_density) cond_number_dens_new[i] = condensates[i]->max_number_density;

    //double_type delta_lambda = ln_tau - std::log(activity_corr_old[i]) - std::log(cond_number_dens_old[i]) - delta_n[i];
    double_type delta_lambda = condensates[i]->log_tau - std::log(activity_corr_old[i]) - std::log(cond_number_dens_old[i]) - delta_n[i];

    if (delta_lambda > max_change) delta_lambda = max_change;
    if (delta_lambda < -max_change) delta_lambda = -max_change;

    activity_corr_new[i] = activity_corr_old[i] * std::exp(delta_lambda);
  }


  for (size_t i=0; i<elements.size(); ++i)
  { 
    double_type delta_n = result[i + nb_cond_jac];

    if (delta_n > max_change) delta_n = max_change;
    if (delta_n < -max_change) delta_n = -max_change;

    if (std::fabs(delta_n) > max_delta) max_delta = std::fabs(delta_n);

    elem_number_dens_new[i] = elem_number_dens_old[i] * std::exp(delta_n);

    //if (elem_number_dens_new[i] < options.element_density_minlimit) elem_number_dens_new[i] = options.element_density_minlimit;
  }


  return max_delta;
}


/*template <class double_type>
double_type CondensedPhase<double_type>::correctValues(
  const std::vector<double_type>& result,
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
  const double_type tau,
  const double max_change)
{
  std::vector<double_type> delta_n(condensates.size(), 0);

  for (size_t i=0; i<condensates_jac.size(); ++i)
    delta_n[condensates_jac[i]] = result[i];


  const size_t nb_cond_jac = condensates_jac.size();

  for (size_t i=0; i<condensates_rem.size(); ++i)
  {
    const unsigned int index = condensates_rem[i];

    for (size_t j=0; j<elements.size(); ++j)
      delta_n[index] += condensates[index]->stoichiometric_vector[elements[j]->index] * result[j+nb_cond_jac];

    delta_n[index] /= activity_corr_old[index];

    delta_n[index] += condensates[index]->log_activity / activity_corr_old[index] + tau/(cond_number_dens_old[index] * activity_corr_old[index]);
  }


  double_type max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {
    if (delta_n[i] > max_change) delta_n[i] = max_change;
    if (delta_n[i] < -max_change) delta_n[i] = -max_change;

    if (std::fabs(delta_n[i]) > max_delta) max_delta = std::fabs(delta_n[i]);


    cond_number_dens_new[i] = cond_number_dens_old[i] * std::exp(delta_n[i]);
    //if (cond_number_dens_new[i] > condensates[i]->max_number_density) cond_number_dens_new[i] = condensates[i]->max_number_density;

    double_type delta_lambda = -1.0 + tau/(cond_number_dens_old[i] * activity_corr_old[i]) - delta_n[i];

    if (delta_lambda > max_change) delta_lambda = max_change;
    if (delta_lambda < -max_change) delta_lambda = -max_change;

    activity_corr_new[i] = activity_corr_old[i] * std::exp(delta_lambda);
  }


  for (size_t i=0; i<elements.size(); ++i)
  { 
    double_type delta_n = result[i + nb_cond_jac];

    if (delta_n > max_change) delta_n = max_change;
    if (delta_n < -max_change) delta_n = -max_change;

    if (std::fabs(delta_n) > max_delta) max_delta = std::fabs(delta_n);

    elem_number_dens_new[i] = elem_number_dens_old[i] * std::exp(delta_n);

    //if (elem_number_dens_new[i] < options.element_density_minlimit) elem_number_dens_new[i] = options.element_density_minlimit;
  }


  return max_delta;
}*/





template class CondensedPhase<double>;
template class CondensedPhase<long double>;
}