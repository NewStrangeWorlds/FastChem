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
#include "../species_struct.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//Log-space residual for the element conservation equation.
//Returns ln_P = ln(P(y_j)), ln_dP = ln(dP/dy_j), and R,
//where P is the positive-definite sum of free-atom and molecule contributions
//and R is the required density.
//The caller computes the Newton step from these without overflow:
//  delta = -(P - R) / dP = -exp(ln_P - ln_dP) + R * exp(-ln_dP)
//When use_all_molecules is false (standard mode), only major molecules (same abundance)
//are included, and R = phi*n_gas - n_min - n_maj.
//When use_all_molecules is true (alternative/backup mode), all molecules in the element's
//molecule_list are used, and R = phi*n_gas - n_exc.
void GasPhaseSolver::logSpaceResidual(
  const Element& species,
  const std::vector<Element>& elements,
  const std::vector<Molecule>& molecules,
  const double gas_density,
  const double y_j,
  double& ln_P,
  double& ln_dP,
  double& R,
  const bool use_all_molecules)
{
  if (use_all_molecules)
  {
    //Alternative mode: R = phi*n_gas - n_exc
    //n_exc = phi * Sigma_{i: nu_ij=0} sigma_i * n_i (molecules not containing element j)
    double n_exc = 0.0;

    for (size_t i = 0; i < molecules.size(); ++i)
      if (molecules[i].stoichiometric_vector[species.index] == 0)
        n_exc += molecules[i].sigma * molecules[i].number_density;

    n_exc *= species.phi;
    R = species.phi * gas_density - n_exc;
  }
  else
  {
    //Standard mode: R = phi*n_gas - n_min - n_maj
    R = species.phi * gas_density - species.number_density_min - species.number_density_maj;
  }

  //Build the positive sum P(y_j) and its derivative dP/dy_j
  //P(y_j) = exp(y_j) + Sigma_i kappa_ij * exp(mac_i + Sigma_l nu_il * y_l)
  //dP/dy_j = exp(y_j) + Sigma_i nu_ij * kappa_ij * exp(mac_i + Sigma_l nu_il * y_l)
  std::vector<double> log_terms;
  std::vector<double> coeffs_P;
  std::vector<double> coeffs_dP;

  //Free atom term
  log_terms.push_back(y_j);
  coeffs_P.push_back(1.0);
  coeffs_dP.push_back(1.0);

  for (auto & i : species.molecule_list)
  {
    const int nu_j = molecules[i].stoichiometric_vector[species.index];

    if (nu_j < 1) continue;

    if (!use_all_molecules && molecules[i].abundance != species.abundance)
      continue;

    //log(n_i) = mac_i + Sigma_l nu_il * y_l
    double log_n = molecules[i].mass_action_constant;

    for (auto & l : molecules[i].element_indices)
    {
      if (l == species.index)
        log_n += molecules[i].stoichiometric_vector[l] * y_j;
      else if (molecules[i].stoichiometric_vector[l] != 0)
        log_n += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    const double kappa = nu_j + species.phi * molecules[i].sigma;

    log_terms.push_back(log_n);
    coeffs_P.push_back(kappa);
    coeffs_dP.push_back(static_cast<double>(nu_j) * kappa);
  }

  //Compute ln(P) and ln(dP/dy) using logSumExp for numerical stability
  ln_P = logSumExp(log_terms, coeffs_P);
  ln_dP = logSumExp(log_terms, coeffs_dP);
}



//Am coefficients for electron density computation
//Uses log_number_density to avoid log(0) for underflowed densities
double GasPhaseSolver::AmCoeffElectron(
  const Element& electron,
  const std::vector< Element >& elements,
  const std::vector< Molecule >& molecules,
  const int order)
{
  double Am = 0.0;

  for (auto & i : electron.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[electron.index] == order)
    {
      double sum = 0;

      for (auto & j : molecules[i].element_indices)
      {
        if (j != electron.index && molecules[i].stoichiometric_vector[j] != 0)
          sum += molecules[i].stoichiometric_vector[j] * elements[j].log_number_density;
      }

      Am += safeExp(molecules[i].mass_action_constant + sum) * order;
    }
  }


  return Am;
}


}


