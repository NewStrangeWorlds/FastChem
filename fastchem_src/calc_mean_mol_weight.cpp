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

#include <vector>
#include <cmath>



namespace fastchem {


//Calculate the mean molecular weight of the converged network
//See Eq. (2.10)
template <class double_type>
double FastChem<double_type>::calcMeanMolecularWeight(const double total_density, const unsigned int grid_index)
{
   double mean_molecular_weight = 0.0;

   for(size_t i=0; i<nb_species; ++i)
     mean_molecular_weight += species[i]->molecular_weight * species[i]->number_density[grid_index];

   mean_molecular_weight /= total_density;


   return mean_molecular_weight;
}



template class FastChem<double>;
template class FastChem<long double>;


}
