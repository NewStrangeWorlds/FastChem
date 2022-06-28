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


//This is the backup solver for the standard analytic FastChem solvers
//Currently configured to use a 1D Newton solver, that employs Eq. (2.34)
template <class double_type>
void GasPhaseSolver<double_type>::backupSol(
  Element<double_type>& species,
  std::vector<Element<double_type>>& elements,
  const std::vector<Molecule<double_type>>& molecules, 
  const double_type gas_density)
{

  newtonSol(species, elements, molecules, gas_density, true);
  
}


template class GasPhaseSolver<double>;
template class GasPhaseSolver<long double>;
}



