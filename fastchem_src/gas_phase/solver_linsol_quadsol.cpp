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
void GasPhaseSolver::inertSol(
  Element& species,
  std::vector< Element >& elements,
  const std::vector< Molecule >& molecules,
  const double gas_density)
{
  double n = species.phi * gas_density - species.number_density_min - species.number_density_maj;

  if (n > 0)
  {
    species.log_number_density = std::log(n);
    species.number_density = n;
  }
  else
  {
    species.log_number_density = static_cast<double>(LOG_DENSITY_FLOOR);
    species.number_density = 0.0;
  }
}


//Analytic solution for linear equation in log-space
//See Paper I, Sect. 2.4.2 and Eq. (2.32)
//Computes y = ln(R) - ln(A1) where R = phi*n_gas - n_min - n_maj
//and ln(A1) is computed via logSumExp for numerical stability
void GasPhaseSolver::linSol(
  Element& species,
  std::vector< Element >& elements,
  const std::vector< Molecule >& molecules,
  const double gas_density)
{
  double R = species.phi * gas_density - species.number_density_min - species.number_density_maj;

  if (R <= 0)
  {
    species.log_number_density = static_cast<double>(LOG_DENSITY_FLOOR);
    species.number_density = 0.0;
    return;
  }

  //Build ln(A1) = logSumExp of molecule terms + free-atom term
  //A1 = 1 + Sigma_i kappa_i * exp(mac_i + Sigma_{l!=j} nu_il * y_l)
  scratch_log_terms_.clear();
  scratch_coeffs_P_.clear();

  //Free atom contribution: exp(0) * 1 = 1
  scratch_log_terms_.push_back(0.0);
  scratch_coeffs_P_.push_back(1.0);

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[species.index] != 1) continue;
    if (molecules[i].abundance != species.abundance) continue;

    double log_term = molecules[i].mass_action_constant;

    for (auto & l : molecules[i].element_indices)
    {
      if (l != species.index && molecules[i].stoichiometric_vector[l] != 0)
        log_term += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    const double kappa = 1.0 + species.phi * molecules[i].sigma;

    scratch_log_terms_.push_back(log_term);
    scratch_coeffs_P_.push_back(kappa);
  }

  double ln_A1 = logSumExp(scratch_log_terms_, scratch_coeffs_P_);

  species.log_number_density = std::log(R) - ln_A1;
  species.number_density = safeExp(species.log_number_density);
}


//Analytic solution for quadratic equation in log-space
//See Paper I, Sect. 2.4.2 and Eq. (2.32)
//Computes A0 + A1*n + A2*n^2 = 0 entirely in log-space
//Falls back to linSol if A2 underflows
void GasPhaseSolver::quadSol(
  Element& species,
  std::vector< Element >& elements,
  const std::vector< Molecule >& molecules,
  const double gas_density)
{
  double R = species.phi * gas_density - species.number_density_min - species.number_density_maj;

  if (R <= 0)
  {
    species.log_number_density = static_cast<double>(LOG_DENSITY_FLOOR);
    species.number_density = 0.0;
    return;
  }

  //Compute ln(A1) and ln(A2) via logSumExp
  scratch_log_terms_.clear();
  scratch_coeffs_P_.clear();
  scratch_log_terms_2_.clear();
  scratch_coeffs_dP_.clear();

  //Free atom contribution to A1
  scratch_log_terms_.push_back(0.0);
  scratch_coeffs_P_.push_back(1.0);

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].abundance != species.abundance) continue;

    const int nu_j = molecules[i].stoichiometric_vector[species.index];

    if (nu_j != 1 && nu_j != 2) continue;

    double log_term = molecules[i].mass_action_constant;

    for (auto & l : molecules[i].element_indices)
    {
      if (l != species.index && molecules[i].stoichiometric_vector[l] != 0)
        log_term += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    const double kappa = nu_j + species.phi * molecules[i].sigma;

    if (nu_j == 1)
    {
      scratch_log_terms_.push_back(log_term);
      scratch_coeffs_P_.push_back(kappa);
    }
    else
    {
      scratch_log_terms_2_.push_back(log_term);
      scratch_coeffs_dP_.push_back(kappa);
    }
  }

  double ln_A1 = logSumExp(scratch_log_terms_, scratch_coeffs_P_);
  double ln_A2 = logSumExp(scratch_log_terms_2_, scratch_coeffs_dP_);

  //If A2 underflows, fall back to linSol
  if (ln_A2 <= static_cast<double>(LOG_DENSITY_FLOOR))
  {
    if (options.verbose_level >= 3)
      std::cout << "FastChem: Underflow in QuadSol for species "
                << species.symbol << " : switching to LinSol.\n";

    linSol(species, elements, molecules, gas_density);
    return;
  }

  double ln_R = std::log(R);

  //discriminant = A1^2 + 4*A2*R (always positive)
  std::vector<double> disc_log = {2.0 * ln_A1, std::log(4.0) + ln_A2 + ln_R};
  std::vector<double> disc_coeff = {1.0, 1.0};
  double ln_disc = logSumExp(disc_log, disc_coeff);

  //Q = -0.5 * (A1 + sqrt(disc))
  //ln|Q| = ln(0.5) + logSumExp([ln_A1, 0.5*ln_disc], [1, 1])
  std::vector<double> q_log = {ln_A1, 0.5 * ln_disc};
  std::vector<double> q_coeff = {1.0, 1.0};
  double ln_Q = std::log(0.5) + logSumExp(q_log, q_coeff);

  //n = R / |Q|, so y = ln(R) - ln|Q|
  species.log_number_density = ln_R - ln_Q;
  species.number_density = safeExp(species.log_number_density);
}


}


