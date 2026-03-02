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

#include "gas_phase.h"

#include "solver.h"


namespace fastchem {



//Calculation of the electron density
template <class double_type>
void GasPhase<double_type>::calculateElectronDensities(
  Element<double_type>& electron,
  const double_type& old_number_density,
  const double_type gas_density)
{
  //Am I the electron? 
  if (electron.symbol != "e-") return;

  //Electron log density floor: use log(element_density_minlimit) instead of LOG_DENSITY_FLOOR
  //because cations have negative electron stoichiometry (nu_e < 0), and -nu_e * LOG_DENSITY_FLOOR
  //would create astronomically large ion contributions in the element conservation equation.
  const double_type electron_log_floor = std::log(options.element_density_minlimit);

  //no ions present
  if (electron.molecule_list.size() == 0)
  {
    electron.number_density = 0.0;
    electron.log_number_density = electron_log_floor;
    return;
  } 


  //if we have't determined the maximum order of cations and anions, we do so now
  if (solver.order_anion == -999 && solver.order_cation == -999)
  {
    solver.order_cation = 0;

    for (auto & i : electron.molecule_list)
      if (molecules[i].stoichiometric_vector[electron.index] < solver.order_cation )
        solver.order_cation = molecules[i].stoichiometric_vector[electron.index];

    solver.order_cation = std::abs(solver.order_cation);


    solver.order_anion = 0;

    for (auto & i : electron.molecule_list)
      if (molecules[i].stoichiometric_vector[electron.index] > solver.order_anion )
        solver.order_anion = molecules[i].stoichiometric_vector[electron.index];

    solver.order_anion = std::abs(solver.order_anion);
  }


  //for singly-ionised species we use the analytic solution
  if (electron.solver_order == 1)
    calculateSinglyIonElectrons(electron, old_number_density);
  else
    calculateMultIonElectrons(electron, old_number_density, gas_density);
}


//Calculation of the electron density for at most singly-ionised species
//Uses the analytical solution from Paper 1, Appendix B
template <class double_type>
void GasPhase<double_type>::calculateSinglyIonElectrons(
  Element<double_type>& electron, const double_type& old_number_density)
{
  const unsigned int index = electron.index;

  //Accumulate alpha and (1+beta) in log space to avoid underflow in sqrt(alpha/(1+beta))
  //log(0) = -inf for alpha; log(1) = 0 for (1+beta)
  const double_type neg_inf = -std::numeric_limits<double_type>::infinity();
  double_type log_alpha = neg_inf;
  double_type log_one_plus_beta = 0.0;

  for (auto & i : electron.molecule_list)
  {
    double_type sum = 0;

    for (auto & j : molecules[i].element_indices)
      if (j != electron.index && molecules[i].stoichiometric_vector[j] != 0)
        sum += molecules[i].stoichiometric_vector[j] * elements[j].log_number_density;

    double_type log_term = molecules[i].mass_action_constant + sum;

    //the anions, Eq. (B3) in Paper 1
    if (molecules[i].stoichiometric_vector[index] == 1)
      log_one_plus_beta = logAddExp(log_one_plus_beta, log_term);
    else if (molecules[i].stoichiometric_vector[index] == -1)  //the cations, Eq. (B4) in Paper 1
      log_alpha = logAddExp(log_alpha, log_term);
  }

  //Eq. (B2) in Paper 1: n_e = sqrt(alpha / (1 + beta))
  //In log space: log(n_e) = 0.5 * (log(alpha) - log(1 + beta))
  const double_type electron_log_floor = std::log(options.element_density_minlimit);

  double_type electron_density = 0.0;
  double_type log_electron_density = electron_log_floor;

  if (log_alpha > neg_inf)
  {
    log_electron_density = 0.5 * (log_alpha - log_one_plus_beta);
    electron_density = std::exp(log_electron_density);
  }
  
  elements[e_].number_density = electron_density;
  elements[e_].log_number_density = (electron_density > 0) ? log_electron_density : electron_log_floor;
}




//Calculation of the electron density, based on charge conservation
//This approach is used for multi-ionised species
//First tries to estimate the electron density via Paper 1, Eq. (2.35).
//In case that fails (electron density not sufficiently high enough), it switches to a 1D Newton's method.
//See Sect. 2.4.3 for details.
template <class double_type>
void GasPhase<double_type>::calculateMultIonElectrons(
  Element<double_type>& electron, const double_type& old_number_density, const double_type& gas_density)
{
  electron.number_density = 0.0;


  double_type positive_ion_density = 0;
  double_type negative_ion_density = 0;

  for (auto & i : electron.molecule_list)
    if (molecules[i].stoichiometric_vector[electron.index] > 0)
      negative_ion_density += molecules[i].stoichiometric_vector[e_] * molecules[i].number_density;
    else
      positive_ion_density -= molecules[i].stoichiometric_vector[e_] * molecules[i].number_density;


  double_type electron_density = positive_ion_density - negative_ion_density;


  const double_type electron_log_floor = std::log(options.element_density_minlimit);

  double_type delta = 0.9;

  if (electron_density > delta*positive_ion_density)
  {
    electron.number_density = std::sqrt(electron_density * old_number_density);
    electron.log_number_density = (electron.number_density > 0)
      ? std::log(electron.number_density) : electron_log_floor;
  }
  else
  {
    //switching to Newton's method
    solver.newtonSolElectron(electron, elements, molecules, gas_density);
  }
}



template class GasPhase<double>;
template class GasPhase<long double>;
}



