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

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>


namespace fastchem {


//Bisection method in log-space
//Bisects on y in [LOG_DENSITY_FLOOR, ln(gas_density)] using the log-space residual.
//The sign of F = P - R is determined by comparing ln_P with ln(R),
//avoiding any conversion to linear space that could overflow.
bool GasPhaseSolver::bisection(
  Element& species,
  std::vector<Element>& elements,
  const std::vector<Molecule>& molecules,
  const double gas_density,
  const bool use_all_molecules)
{
  double y_lo = static_cast<double>(LOG_DENSITY_FLOOR);
  double y_hi = std::log(gas_density);

  bool converged = false;

  //Helper: compute the sign of F = P - R entirely in log-space
  //Returns +1 if P > R (y too large), -1 if P < R (y too small), 0 if equal
  auto signF = [&](double y) -> int {
    double ln_P, ln_dP, R;
    logSpaceResidual(species, elements, molecules, gas_density, y, ln_P, ln_dP, R, use_all_molecules);

    if (R <= 0) return 1;  //P > 0 >= R, so F = P - R > 0

    double ln_R = std::log(R);
    if (ln_P > ln_R + 1e-12) return 1;
    if (ln_P < ln_R - 1e-12) return -1;
    return 0;
  };

  int sign_lo = signF(y_lo);
  int sign_hi = signF(y_hi);

  //Check that the root is bracketed
  if (sign_lo * sign_hi > 0)
  {
    //Root not bracketed, use midpoint as best guess
    species.log_number_density = 0.5 * (y_lo + y_hi);
    species.number_density = safeExp(species.log_number_density);

    if (options.verbose_level >= 3)
      std::cout << "FastChem: WARNING: Bisection root not bracketed for "
                << species.symbol << "\n";

    return false;
  }


  for (unsigned int iter_step = 0; iter_step < options.nb_max_bisection_iter; ++iter_step)
  {
    const double y_mid = 0.5 * (y_lo + y_hi);

    int sign_mid = signF(y_mid);

    if (sign_mid == 0)
    {
      y_lo = y_mid;
      y_hi = y_mid;
      converged = true;
      break;
    }

    if (sign_mid * sign_lo < 0)
      y_hi = y_mid;
    else
      y_lo = y_mid;

    //Convergence test
    if (std::fabs(y_hi - y_lo) < options.chem_accuracy * 1e-3)
    {
      converged = true;
      break;
    }
  }


  species.log_number_density = 0.5 * (y_lo + y_hi);
  species.number_density = safeExp(species.log_number_density);


  if (!converged && options.verbose_level >= 3)
    std::cout << "Bisection iteration limit reached, result may not be optimal."
              << "\t" << y_lo << "\t" << y_hi
              << "\t" << std::fabs(y_hi - y_lo) << "\t" << options.chem_accuracy * 1e-3  << "\n";


  return converged;
}


}
