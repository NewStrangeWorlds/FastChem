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


  double_type tau = 1e-15;

  for (auto & i : condensates_act)
  {
    i->tau = tau * element_data.elements[i->reference_element].epsilon * total_element_density;
    i->log_tau = std::log(i->tau);
  }


  for (auto & i : condensates_act)
  {
    if (i->number_density == 0 && i->activity_correction == 0)
    {
       i->number_density = i->max_number_density; //1e-10;
       i->activity_correction = 1.0; //i->log_activity; //1; //1.0;
       //i->activity_correction = i->tau/i->number_density; //1; //1.0;

       //std::cout << i->symbol << "\t" << i->number_density << "\t" << i->activity_correction << "\n"; 
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

  double_type limit = 5.0;

  for (nb_iterations=0; nb_iterations<10000; ++nb_iterations)
  { 
    selectJacobianCondensates(
      condensates_act,
      cond_densities_old,
      activity_corr_old,
      condensates_jac,
      condensates_rem);

    double_type objective_function_0 = 0;
    Eigen::VectorXdt<double_type> scaling_factors;
    Eigen::VectorXdt<double_type> result;

    const bool system_invertible = solver.newtonStep(
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

    // std::vector<double_type> result = solver.newtonStepFull(
    //   condensates_act,
    //   elements_cond,
    //   molecules,
    //   total_element_density,
    //   cond_densities_old,
    //   elem_densities_old,
    //   activity_corr_old,
    //   scaling_factors,
    //   objective_function_0);

    Eigen::VectorXdt<double_type> result_scaled = result;

    double_type max_value = result.cwiseAbs().maxCoeff();

    if (max_value > limit)
      result *= limit/max_value;

    double_type max_delta = 0;


    max_delta = correctValues(
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
      limit);

    // max_delta = correctValuesFull(
    //   result,
    //   condensates_act,
    //   activity_corr_old,
    //   activity_corr_new,
    //   cond_densities_old,
    //   cond_densities_new,
    //   elements_cond,
    //   elem_densities_old,
    //   elem_densities_new,
    //   limit);

    for (size_t i=0; i<elements_cond.size(); ++i)
      elements_cond[i]->number_density = elem_densities_new[i];

    for (auto & i : condensates_act)  i->calcActivity(temperature, elements);

    for (auto & i : molecules)  i.calcNumberDensity(elements);

    // double_type objective_function_1 = solver.objectiveFunction(
    //   condensates_act,
    //   condensates_jac,
    //   condensates_rem,
    //   activity_corr_new,
    //   cond_densities_new,
    //   elements_cond,
    //   molecules,
    //   scaling_factors,
    //   total_element_density);

    // double_type objective_function_1 = solver.objectiveFunctionFull(
    //       condensates_act,
    //       activity_corr_new,
    //       cond_densities_new,
    //       elements_cond,
    //       molecules,
    //       scaling_factors,
    //       total_element_density);


    double_type lambda_prev = 1.0;
    double_type lambda_2prev = 0.0;

    double_type lamdba_min = 1e-5;

    //std::cout << "of_0: " <<  objective_function_0 << "\t" << max_delta << "\n";

    /*if (objective_function_1 > objective_function_0)
    {
      const double_type object_function_deriv = -2.0*objective_function_0;

      double_type objective_function_prev = objective_function_1;
      double_type objective_function_2prev = 0.0;

      double_type lambda = 1;

      double_type limit1 = limit;

      // for (unsigned int i=0; i<4; ++i)
      // {
      //   limit1 *= 0.5;

      //   max_delta = correctValues(
      //     result,
      //     condensates_act,
      //     condensates_jac,
      //     condensates_rem,
      //     activity_corr_old,
      //     activity_corr_new,
      //     cond_densities_old,
      //     cond_densities_new,
      //     elements_cond,
      //     elem_densities_old,
      //     elem_densities_new,
      //     log_tau,
      //     limit1);

      //   for (size_t i=0; i<elements_cond.size(); ++i)
      //     elements_cond[i]->number_density = elem_densities_new[i];

      //   for (auto & i : condensates_act)  i->calcActivity(temperature, elements);

      //   for (auto & i : molecules)  i.calcNumberDensity(elements);

      //   objective_function_2prev = objective_function_prev;
      //   lambda_2prev = lambda_prev;
      //   lambda_prev = lambda;

      //   objective_function_prev = solver.objectiveFunction(
      //     condensates_act,
      //     condensates_jac,
      //     condensates_rem,
      //     activity_corr_new,
      //     cond_densities_new,
      //     elements_cond,
      //     molecules,
      //     scaling_factors,
      //     total_element_density);

      //   std::cout << "l " << i << "\t" << objective_function_prev << "\t" << objective_function_0 << "\t" << "\t" << limit1 << "\n";

      //   if (objective_function_prev < objective_function_0) break;
      // }
      
      while(objective_function_prev > objective_function_0 + 1e-4*lambda*object_function_deriv && lambda_prev > lamdba_min)
      { 

        if (lambda_2prev == 0)
        {
          lambda = -object_function_deriv /(2.0 * (objective_function_prev - objective_function_0 - object_function_deriv));
        }
        else
        {
          Eigen::Matrix<double_type, 2, 2> cubic_fit_matrix;
          Eigen::Matrix<double_type, 2,1> cubic_fit_vector;

          cubic_fit_matrix(0,0) = 1.0/(lambda_prev*lambda_prev);
          cubic_fit_matrix(1,0) = -lambda_2prev/(lambda_prev*lambda_prev);
          cubic_fit_matrix(0,1) = -1.0/(lambda_2prev*lambda_2prev);
          cubic_fit_matrix(1,1) = lambda_prev/(lambda_2prev*lambda_2prev);

          cubic_fit_vector(0) = objective_function_prev - objective_function_0 - object_function_deriv*lambda_prev;
          cubic_fit_vector(1) = objective_function_2prev - objective_function_0 - object_function_deriv*lambda_2prev;

          cubic_fit_vector /= lambda_prev-lambda_2prev;

          Eigen::Matrix<double_type, 2,1> cubic_fit = cubic_fit_matrix*cubic_fit_vector;
          const double_type a = cubic_fit(0);
          const double_type b = cubic_fit(1);

          lambda = (-b + std::sqrt(b*b - 3*a*object_function_deriv))/(3*a);

          //std::cout << "cf " << a << "\t" << b << "\t" << lambda << "\n";
        }

        if (lambda < 0.1*lambda_prev) lambda = 0.1*lambda_prev;
        if (lambda > 0.5*lambda_prev) lambda = 0.5*lambda_prev;

        std::vector<double_type> result1 = result;

        //for (auto & i : result1) 
          //i *= lambda;

        limit1 = limit * lambda;

        // max_delta = correctValues(
        //   result1,
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
        //   log_tau,
        //   limit);

        max_delta = correctValuesFull(
          result1,
          condensates_act,
          activity_corr_old,
          activity_corr_new,
          cond_densities_old,
          cond_densities_new,
          elements_cond,
          elem_densities_old,
          elem_densities_new,
          limit1);

        for (size_t i=0; i<elements_cond.size(); ++i)
          elements_cond[i]->number_density = elem_densities_new[i];

        for (auto & i : condensates_act)  i->calcActivity(temperature, elements);

        for (auto & i : molecules)  i.calcNumberDensity(elements);

        objective_function_2prev = objective_function_prev;
        lambda_2prev = lambda_prev;
        lambda_prev = lambda;

        // objective_function_prev = solver.objectiveFunction(
        //   condensates_act,
        //   condensates_jac,
        //   condensates_rem,
        //   activity_corr_new,
        //   cond_densities_new,
        //   elements_cond,
        //   molecules,
        //   scaling_factors,
        //   total_element_density);

        objective_function_prev = solver.objectiveFunctionFull(
          condensates_act,
          activity_corr_new,
          cond_densities_new,
          elements_cond,
          molecules,
          scaling_factors,
          total_element_density);

        //std::cout << objective_function_prev << "\t" << objective_function_0 << "\t" << objective_function_0 + 1e-4*lambda*object_function_deriv << "\t" << lambda << "\t" << object_function_deriv << "\t" << limit1 << "\n";
      }

      if (lambda_prev <= lamdba_min)
      {
        // max_delta = correctValues(
        //   result,
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
        //   log_tau,
        //   1.0);
        //std::cout << lambda_prev << "\t" << lamdba_min << "\n";
        max_delta = correctValuesFull(
          result,
          condensates_act,
          activity_corr_old,
          activity_corr_new,
          cond_densities_old,
          cond_densities_new,
          elements_cond,
          elem_densities_old,
          elem_densities_new,
          1.0);

        for (size_t i=0; i<elements_cond.size(); ++i)
          elements_cond[i]->number_density = elem_densities_new[i];

        for (auto & i : condensates_act)  i->calcActivity(temperature, elements);

        for (auto & i : molecules)  i.calcNumberDensity(elements);
      }

    }*/
    //if (nb_iterations == 181) exit(0);
    //if (nb_iterations > 900)
    /*{
    std::cout << "iter: " << nb_iterations << "\n";
    for (size_t i=0; i<condensates_act.size(); ++i)
      std::cout <<i << "  " << condensates_act[i]->symbol << "\t" << cond_densities_old[i] << "\t" << cond_densities_new[i] << "\t" << activity_corr_old[i] << "\t" << activity_corr_new[i] << "\t" << condensates_act[i]->log_activity << "\t" << condensates_act[i]->tau << "\n";
    }*/
    //exit(0);
    elem_densities_old = elem_densities_new;
    cond_densities_old = cond_densities_new;
    activity_corr_old = activity_corr_new;
    
    //std::cout << "cond delta " << max_delta << "\t" << objective_function_0 << "\n";
    cond_converged = max_delta < 1e-9;
    //cond_converged = objective_function_0 < 1e-6;
    
    //if (nb_iterations == 2858) exit(0);
    if (cond_converged) break;
  }

  //if (nb_iterations > 9998) exit(0);

  
  for (size_t i=0; i<condensates_act.size(); ++i)
  {
    condensates_act[i]->number_density = cond_densities_new[i];
    condensates_act[i]->activity_correction = activity_corr_new[i];

    //if (condensates_act[i]->number_density <= 1e-15) condensates_act[i]->number_density = 0.0;
    //if (condensates_act[i]->log_activity < -e-1) condensates_act[i]->number_density = 0.0;
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

    if (std::fabs(delta_n_cond[i]) > max_delta) 
      max_delta = std::fabs(delta_n_cond[i]);

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
    // if (activity_corr_old[i] + delta_lambda > 0)
    //   activity_corr_new[i] = activity_corr_old[i] + delta_lambda;
    // else
    //   activity_corr_new[i] = activity_corr_old[i] * (1.0 - 0.9999);
    //activity_corr_new[i] = condensates[i]->tau / cond_number_dens_new[i];

    //if (cond_number_dens_new[i] < 1e-5 && activity_corr_new[i] < 1e-5)
      //activity_corr_new[i] = std::fabs(condensates[i]->log_activity);

    //if (activity_corr_new[i] < condensates[i]->tau / cond_number_dens_new[i]) activity_corr_new[i] = condensates[i]->tau / cond_number_dens_new[i];
  }


  for (size_t i=0; i<elements.size(); ++i)
  { 
    double_type delta_n_elem = result(i + nb_cond_jac);

    if (delta_n_elem > max_change) 
      delta_n_elem = max_change;
    
    if (delta_n_elem < -max_change) 
      delta_n_elem = -max_change;

    if (std::fabs(delta_n_elem) > max_delta) 
      max_delta = std::fabs(delta_n_elem);

    elem_number_dens_new[i] = elem_number_dens_old[i] * std::exp(delta_n_elem);
  }


  return max_delta;
}



template <class double_type>
double_type CondensedPhase<double_type>::correctValuesFull(
  const std::vector<double_type>& result2,
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& activity_corr_old,
  std::vector<double_type>& activity_corr_new,
  const std::vector<double_type>& cond_number_dens_old,
  std::vector<double_type>& cond_number_dens_new,
  const std::vector<Element<double_type>*>& elements,
  const std::vector<double_type>& elem_number_dens_old,
  std::vector<double_type>& elem_number_dens_new,
  const double max_change)
{ 
  std::vector<double_type> result = result2;

  double_type max_value = 0;
  
  for (auto & i : result)
    if (std::abs(i) > max_value) max_value = std::abs(i);

  //std::cout << "max value " <<  max_value << "\n";
  
  if (max_value > max_change)
    for (auto & i : result)
      i = i / max_value * max_change;

  std::vector<double_type> delta_n(condensates.size(), 0);
  std::vector<double_type> delta_lambda(condensates.size(), 0);

  for (size_t i=0; i<condensates.size(); ++i)
  {
    delta_n[i] = result[i];
    delta_lambda[i] = result[i+condensates.size()];
  }


  double_type max_delta = 0;

  for (size_t i=0; i<condensates.size(); ++i)
  {
    if (delta_n[i] > max_change) delta_n[i] = max_change;
    if (delta_n[i] < -max_change) delta_n[i] = -max_change;

    if (std::fabs(delta_n[i]) > max_delta) max_delta = std::fabs(delta_n[i]);


    cond_number_dens_new[i] = cond_number_dens_old[i] * std::exp(delta_n[i]);
    if (cond_number_dens_new[i] > condensates[i]->max_number_density) cond_number_dens_new[i] = condensates[i]->max_number_density;

    if (delta_lambda[i] > max_change) delta_lambda[i] = max_change;
    if (delta_lambda[i] < -max_change) delta_lambda[i] = -max_change;

    activity_corr_new[i] = activity_corr_old[i] * std::exp(delta_lambda[i]);
    //activity_corr_new[i] = activity_corr_old[i] + delta_lambda[i];

    // if ( (activity_corr_old[i] + delta_lambda[i]) > 0)
    //   activity_corr_new[i] = activity_corr_old[i] + delta_lambda[i];
    // else
    //   activity_corr_new[i] = activity_corr_old[i] * (1.0 - 0.9999);

    //std::cout << condensates[i]->symbol << "\t" << delta_lambda << "\t" << delta_n[i] << "\t" << condensates[i]->tau << "\n";
  }


  for (size_t i=0; i<elements.size(); ++i)
  { 
    double_type delta_n = result[i + 2*condensates.size()];

    if (delta_n > max_change) delta_n = max_change;
    if (delta_n < -max_change) delta_n = -max_change;

    if (std::fabs(delta_n) > max_delta) max_delta = std::fabs(delta_n);

    elem_number_dens_new[i] = elem_number_dens_old[i] * std::exp(delta_n);

    //std::cout << i << "\t" << i + 2*condensates.size() << "\t" << elements[i]->symbol << "\n";

    //if (elem_number_dens_new[i] < options.element_density_minlimit) elem_number_dens_new[i] = options.element_density_minlimit;
  }


  return max_delta;
}



template class CondensedPhase<double>;
template class CondensedPhase<long double>;
}