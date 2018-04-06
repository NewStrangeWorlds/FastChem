/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
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


#include "fastchem.h"


#include <algorithm>
#include <vector>
#include <cmath>


namespace fastchem {


//Initial guess for n_<H>, see Sect. 2.3
template <class double_type>
double_type FastChem<double_type>::setInitialHDensity(const double_type total_density, const unsigned int grid_index)
{
  unsigned int H2_ = getMoleculeIndex("H2");
  unsigned int He_ = getElementIndex("He");

  //set initial total total H density
  //we have to treat different cases since H2 or He could be absent
  double_type h_density = 0.0;

  //general case
  if (H2_ != FASTCHEM_UNKNOWN_SPECIES && He_ != FASTCHEM_UNKNOWN_SPECIES)
    h_density = 1./(2. * std::pow(1 + 2.*elements[He_].abundance,2) * std::exp(molecules[H2_].mass_action_constant[grid_index]))
                    * ( (1. + elements[He_].abundance)
                         + 4. * (1 + 2*elements[He_].abundance) * std::exp(molecules[H2_].mass_action_constant[grid_index]) * total_density
                         - std::sqrt(std::pow(1. + elements[He_].abundance, 2)
                         + 4.*(1. + 2.*elements[He_].abundance)*std::exp(molecules[H2_].mass_action_constant[grid_index])*total_density));

  //H2 is present but He is not
  if (H2_ != FASTCHEM_UNKNOWN_SPECIES && He_ == FASTCHEM_UNKNOWN_SPECIES)
    h_density = 1. + 4. * std::exp(molecules[H2_].mass_action_constant[grid_index]) * total_density
                   + std::sqrt( 1. + 4. * std::exp(molecules[H2_].mass_action_constant[grid_index]) * total_density )
                     / (2. * std::exp(molecules[H2_].mass_action_constant[grid_index]));

  //He is present but H2 is not
  if (H2_ == FASTCHEM_UNKNOWN_SPECIES && He_ != FASTCHEM_UNKNOWN_SPECIES)
    h_density = total_density / (1. + elements[He_].abundance);


  //He and H2 are both not present
  if (H2_ == FASTCHEM_UNKNOWN_SPECIES && He_ == FASTCHEM_UNKNOWN_SPECIES)
    h_density = total_density;


  return h_density;
}


template class FastChem<double>;
template class FastChem<long double>;

}


