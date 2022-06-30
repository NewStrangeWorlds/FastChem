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
  const double temperature,
  const double density,
  std::vector<Molecule<double_type>>& molecules,
  unsigned int& nb_iterations)
{
  std::vector<Condensate<double_type>*> condensates_act;
  std::vector<Element<double_type>*> elements_cond;

  selectActiveCondensates(condensates_act, elements_cond);

  Eigen::MatrixXdt<double_type> jacobian = solver.assembleJacobian(
    condensates_act, 0, elements_cond, molecules);
  
  nb_iterations = 1;

  return true;
} 



template class CondensedPhase<double>;
template class CondensedPhase<long double>;
}