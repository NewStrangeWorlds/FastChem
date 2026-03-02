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


#include "solver.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//Solver for an element that is not part of other species
//See Paper I, Eq. (2.32)
template <class double_type>
void GasPhaseSolver<double_type>::intertSol(
  Element<double_type>& species,
  std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules,
  const double_type gas_density)
{
  double_type n = species.phi * gas_density - species.number_density_min - species.number_density_maj;

  if (n > 0)
  {
    species.log_number_density = std::log(n);
    species.number_density = n;
  }
  else
  {
    species.log_number_density = static_cast<double_type>(LOG_DENSITY_FLOOR);
    species.number_density = 0.0;
  }
}


//Analytic solution for linear equation in log-space
//See Paper I, Sect. 2.4.2 and Eq. (2.32)
//Computes y = ln(R) - ln(A1) where R = phi*n_gas - n_min - n_maj
//and ln(A1) is computed via logSumExp for numerical stability
template <class double_type>
void GasPhaseSolver<double_type>::linSol(
  Element<double_type>& species,
  std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules,
  const double_type gas_density)
{
  double_type R = species.phi * gas_density - species.number_density_min - species.number_density_maj;

  if (R <= 0)
  {
    species.log_number_density = static_cast<double_type>(LOG_DENSITY_FLOOR);
    species.number_density = 0.0;
    return;
  }

  //Build ln(A1) = logSumExp of molecule terms + free-atom term
  //A1 = 1 + Sigma_i kappa_i * exp(mac_i + Sigma_{l!=j} nu_il * y_l)
  std::vector<double_type> log_terms;
  std::vector<double_type> coeffs;

  //Free atom contribution: exp(0) * 1 = 1
  log_terms.push_back(0.0);
  coeffs.push_back(1.0);

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[species.index] != 1) continue;
    if (molecules[i].abundance != species.abundance) continue;

    double_type log_term = molecules[i].mass_action_constant;

    for (auto & l : molecules[i].element_indices)
    {
      if (l != species.index && molecules[i].stoichiometric_vector[l] != 0)
        log_term += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    const double_type kappa = 1.0 + species.phi * molecules[i].sigma;

    log_terms.push_back(log_term);
    coeffs.push_back(kappa);
  }

  double_type ln_A1 = logSumExp(log_terms, coeffs);

  species.log_number_density = std::log(R) - ln_A1;
  species.number_density = safeExp(species.log_number_density);
}


//Analytic solution for quadratic equation in log-space
//See Paper I, Sect. 2.4.2 and Eq. (2.32)
//Computes A0 + A1*n + A2*n^2 = 0 entirely in log-space
//Falls back to linSol if A2 underflows
template <class double_type>
void GasPhaseSolver<double_type>::quadSol(
  Element<double_type>& species,
  std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules,
  const double_type gas_density)
{
  double_type R = species.phi * gas_density - species.number_density_min - species.number_density_maj;

  if (R <= 0)
  {
    species.log_number_density = static_cast<double_type>(LOG_DENSITY_FLOOR);
    species.number_density = 0.0;
    return;
  }

  //Compute ln(A1) and ln(A2) via logSumExp
  std::vector<double_type> log_terms_A1, coeffs_A1;
  std::vector<double_type> log_terms_A2, coeffs_A2;

  //Free atom contribution to A1
  log_terms_A1.push_back(0.0);
  coeffs_A1.push_back(1.0);

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].abundance != species.abundance) continue;

    const int nu_j = molecules[i].stoichiometric_vector[species.index];

    if (nu_j != 1 && nu_j != 2) continue;

    double_type log_term = molecules[i].mass_action_constant;

    for (auto & l : molecules[i].element_indices)
    {
      if (l != species.index && molecules[i].stoichiometric_vector[l] != 0)
        log_term += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    const double_type kappa = nu_j + species.phi * molecules[i].sigma;

    if (nu_j == 1)
    {
      log_terms_A1.push_back(log_term);
      coeffs_A1.push_back(kappa);
    }
    else
    {
      log_terms_A2.push_back(log_term);
      coeffs_A2.push_back(kappa);
    }
  }

  double_type ln_A1 = logSumExp(log_terms_A1, coeffs_A1);
  double_type ln_A2 = logSumExp(log_terms_A2, coeffs_A2);

  //If A2 underflows, fall back to linSol
  if (ln_A2 <= static_cast<double_type>(LOG_DENSITY_FLOOR))
  {
    if (options.verbose_level >= 3)
      std::cout << "FastChem: Underflow in QuadSol for species "
                << species.symbol << " : switching to LinSol.\n";

    linSol(species, elements, molecules, gas_density);
    return;
  }

  double_type ln_R = std::log(R);

  //discriminant = A1^2 + 4*A2*R (always positive)
  std::vector<double_type> disc_log = {2.0 * ln_A1, std::log(4.0) + ln_A2 + ln_R};
  std::vector<double_type> disc_coeff = {1.0, 1.0};
  double_type ln_disc = logSumExp(disc_log, disc_coeff);

  //Q = -0.5 * (A1 + sqrt(disc))
  //ln|Q| = ln(0.5) + logSumExp([ln_A1, 0.5*ln_disc], [1, 1])
  std::vector<double_type> q_log = {ln_A1, 0.5 * ln_disc};
  std::vector<double_type> q_coeff = {1.0, 1.0};
  double_type ln_Q = std::log(0.5) + logSumExp(q_log, q_coeff);

  //n = R / |Q|, so y = ln(R) - ln|Q|
  species.log_number_density = ln_R - ln_Q;
  species.number_density = safeExp(species.log_number_density);
}


template class GasPhaseSolver<double>;
template class GasPhaseSolver<long double>;
}


