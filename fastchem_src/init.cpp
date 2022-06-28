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


#include <cmath>
#include <limits>

#include "fastchem.h"


namespace fastchem {


template <class double_type>
void FastChem<double_type>::init()
{
  //init accuracy limit based on the numerical precision
  if (std::numeric_limits<double_type>::max_exponent10 > 1000)
  {
    options.element_density_minlimit = 1e-512L;
    options.molecule_density_minlimit = 1e-512L;
  }
  else
  {
    options.element_density_minlimit = 1e-155;
    options.molecule_density_minlimit = 1e-155;
  }

}



template class FastChem<double>;
template class FastChem<long double>;
}
