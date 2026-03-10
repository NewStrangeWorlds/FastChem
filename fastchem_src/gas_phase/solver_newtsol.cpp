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


//Newton's method in log-space to solve for element densities
//Iterates on y_j = ln(n_j) using the log-space residual.
//logSpaceResidual returns ln_P and ln_dP (log of P and dP/dy),
//and R (the required density). The Newton delta is computed as:
//  delta = -(P - R) / dP = -exp(ln_P - ln_dP) + R * exp(-ln_dP)
//This avoids overflow when mass action constants are huge (low T).
void GasPhaseSolver::newtonSol(
  Element& species,
  std::vector<Element>& elements,
  const std::vector<Molecule>& molecules,
  const double gas_density,
  const bool use_alternative)
{ 
  const double max_step = 10.0;

  //Compute R: the required density from the element conservation equation
  double R;

  if (use_alternative)
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

  double y;

  //Use the current log_number_density as initial guess if it's reasonable.
  //This avoids recomputing an initial guess that may be dominated by molecules
  //with extreme mass action constants at low temperatures.
  const bool have_previous = (species.log_number_density > static_cast<double>(LOG_DENSITY_FLOOR) + 1.0
                              && species.log_number_density < std::log(gas_density) + 1.0);

  if (R > 0)
  {
    if (have_previous)
    {
      y = species.log_number_density;
    }
    else
    {
      //Initial guess: find the tightest upper bound on y across all molecules.
      //From each molecule i with stoichiometric coefficient nu_j:
      //  kappa_i * exp(mac_i + nu_j * y) <= R
      //  y <= (ln(R) - ln(kappa_i) - mac_i_other) / nu_j
      //The free atom term gives y <= ln(R).
      //Taking the minimum gives a starting point close to the root.
      double ln_R = std::log(R);
      y = ln_R;

      for (auto & i : species.molecule_list)
      {
        const int nu_j = molecules[i].stoichiometric_vector[species.index];
        if (nu_j < 1) continue;
        if (!use_alternative && molecules[i].abundance != species.abundance) continue;

        double mac_other = molecules[i].mass_action_constant;

        for (auto & l : molecules[i].element_indices)
        {
          if (l != species.index && molecules[i].stoichiometric_vector[l] != 0)
            mac_other += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
        }

        const double kappa = static_cast<double>(nu_j) + species.phi * molecules[i].sigma;
        double y_upper = (ln_R - std::log(kappa) - mac_other) / static_cast<double>(nu_j);

        if (y_upper < y) y = y_upper;
      }

      //Clamp to a reasonable range
      if (y < static_cast<double>(LOG_DENSITY_FLOOR))
        y = static_cast<double>(LOG_DENSITY_FLOOR);
    }
  }
  else
  {
    //R <= 0: start from log(gas_density) and converge from above
    y = have_previous ? species.log_number_density : std::log(gas_density);
  }


  //Newton iteration in log-space using logSpaceResidual
  bool converged = false;

  for (unsigned int mu = 0; mu < options.nb_max_newton_iter; ++mu)
  {
    double ln_P, ln_dP, R_res;
    logSpaceResidual(species, elements, molecules, gas_density, y, ln_P, ln_dP, R_res, use_alternative);

    //Compute delta = -(P - R) / dP_dy in a numerically stable way
    //delta = -exp(ln_P - ln_dP) + R * exp(-ln_dP)
    //The first term (P/dP) is always moderate (~1/nu for dominant species).
    //The second term (R/dP) vanishes when far from the root (P >> R).
    double delta;

    if (ln_dP < static_cast<double>(LOG_DENSITY_FLOOR) + 100)
    {
      //dP/dy is essentially zero — can't compute a meaningful step.
      //Use sign of (ln_P - ln(R)) to determine direction.
      if (R_res > 0 && ln_P < std::log(R_res))
        delta = max_step;   //P < R, need to increase y
      else
        delta = -max_step;  //P > R, need to decrease y
    }
    else
    {
      double P_over_dP = safeExp(ln_P - ln_dP);
      double R_over_dP = R_res * safeExp(-ln_dP);
      delta = -(P_over_dP - R_over_dP);
    }

    //Damp large steps
    if (delta > max_step) delta = max_step;
    if (delta < -max_step) delta = -max_step;

    double y_new = y + delta;

    if (std::fabs(delta) < options.newton_err)
    {
      y = y_new;
      converged = true;
      break;
    }

    y = y_new;

    if (std::isnan(y)) break;
  }

  if (converged)
  {
    species.log_number_density = y;
    species.number_density = safeExp(y);
  }


  //In case the normal Newton solver does not converge, switch to other solvers
  if (!converged)
  {
    if (!use_alternative)
    {
      newtonSol(species, elements, molecules, gas_density, true);

      if (options.verbose_level >= 3)
        std::cout << "FastChem: WARNING: NewtSol failed for species "
          << species.symbol << " switched to Backup "
          << y << "\t"
          << species.log_number_density << "\n";
    }
    else
    {
      bisection(species, elements, molecules, gas_density, true);

      if (options.verbose_level >= 3)
        std::cout << "FastChem: WARNING: NewtSol Alt failed for species "
          << species.symbol << " switched to Bisection as backup "
          << y << "\t"
          << species.log_number_density << "\n";
    }
  }
}



//Newton's method for the electrons
//Instead of element conservation, solves for charge balance
void GasPhaseSolver::newtonSolElectron(
  Element& species,
  std::vector<Element>& elements,
  const std::vector<Molecule>& molecules,
  const double gas_density)
{
  //Calculation of the polynomial coefficients
  std::vector<double> Aj_cation(order_cation+1, 0.0);
  std::vector<double> Aj_anion(order_anion+1, 0.0);

  for (int k=1; k<order_cation+1; ++k)
    Aj_cation[k] = AmCoeffElectron(species, elements, molecules, -k);

  for (int k=1; k<order_anion+1; ++k)
    Aj_anion[k] = AmCoeffElectron(species, elements, molecules, k);


  //Newton's method
  bool converged = false;
  double x = order_cation/(1.0 + order_cation) * gas_density; //Initial guess ensures monotonous convergence.


  //one Newton step as lambda function
  auto newton_step = [&] (const double &x)
    {
      //Horner's method for the anions
      double P_anion = Aj_anion[order_anion];
      double P_prime_anion = order_anion*Aj_anion[order_anion];

      for (int k = order_anion-1; k >= 1; --k)
      {
        P_anion = Aj_anion[k] + x * P_anion;
        P_prime_anion = k * Aj_anion[k] + x * P_prime_anion;
      }

      //The cations
      double P_cation = 0.0;
      double P_prime_cation = 0.0;

      for (int k=1; k<order_cation+1; k++)
      {
        P_cation += Aj_cation[k] * std::pow(x, -k);
        P_prime_cation += -k * Aj_cation[k] * std::pow(x, -k-1);
      }

      const double P_j = x  + x * P_anion + P_cation;  //this is the charge balance
      const double P_j_prime = 1.0 + P_prime_cation + P_prime_anion; //derivative

      return x - P_j/P_j_prime; //Newton step
    };


  //Newton iteration
  for (unsigned int mu=0; mu<options.nb_max_newton_iter; ++mu)
  {
    double x_new = newton_step(x);

    if (std::fabs(x_new - x) <= options.newton_err * std::fabs(x_new))  //root found?
    {
      x = x_new;
      converged = true;

      break;
    }

    //prevent x to become negative due to numerical underflow
    if (x_new < 1.e-8*x) x_new = 1.e-8*x;

    x = x_new;


    if (std::isnan(x)) break;
  }


  // Test if root is in (max(0,x*(1-newton_err)),x*(1+newton_err))
  const double x_lower = std::fmax(0., x * (1. - options.newton_err));
  const double x_upper = x * (1. + options.newton_err);


  double P_anion_lower = Aj_anion[order_anion];
  double P_anion_upper = Aj_anion[order_anion];

  for (int k = order_anion-1; k >= 1; --k)
  {
    P_anion_lower = Aj_anion[k] + x_lower * P_anion_lower;
    P_anion_upper = Aj_anion[k] + x_upper * P_anion_upper;
  }

  double P_cation_lower = 0.0;
  double P_cation_upper = 0.0;

  for (int k=1; k<order_cation+1; k++)
  {
    P_cation_lower += Aj_cation[k] * std::pow(x_lower, -k);
    P_cation_upper += Aj_cation[k] * std::pow(x_upper, -k);
  }

  const double P_j_lower = x_lower  + x_lower * P_anion_lower + P_cation_lower;
  const double P_j_upper = x_upper  + x_upper * P_anion_upper + P_cation_upper;


  if (converged && x > 0 && P_j_lower*P_j_upper <= 0.)
  {
    species.number_density = x;
    species.log_number_density = std::log(x);
  }
  else
  {
    species.number_density = x;
    //Use log(element_density_minlimit) as floor for electrons instead of LOG_DENSITY_FLOOR
    //to prevent cation contributions from blowing up via -nu_e * LOG_DENSITY_FLOOR
    const double electron_log_floor = std::log(options.element_density_minlimit);
    species.log_number_density = (x > 0) ? std::log(x) : electron_log_floor;
  }


  //in case something went wrong again, we try to use another backup
  if (x < 0 || !converged || P_j_lower*P_j_upper > 0.)
  {
    const double init = std::log(std::fabs(x));
    nelderMeadElectron(species, elements, molecules, init, 0.0);

    if (options.verbose_level >= 3)
      std::cout << "FastChem: WARNING: NewtSol failed for electrons, switching to Nelder-Mead Backup "
        << x << "\t"
        << species.number_density << "\n";
  }
}


}
