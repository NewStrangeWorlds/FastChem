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
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>


namespace fastchem {


//Mass action law for a single element to be used in the bisection method
template <class double_type>
double_type FastChem<double_type>::bisectionFunction(Element<double_type>& species, const double_type x, const double_type h_density, const unsigned int grid_index)
{
  double_type f_i = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    unsigned int molecule_index = species.molecule_list[j];


    double_type f_j = 0.0;


    for (size_t l=0; l<molecules[molecule_index].element_indices.size(); ++l)
    {
      unsigned int element_index = molecules[molecule_index].element_indices[l];

      double_type species_density = std::log(elements[element_index].number_density[grid_index]);

      //just in case...
      if (std::isnan(species_density) || std::isinf(species_density)) species_density = 0.0;

      if (element_index != species.index)
        f_j += molecules[molecule_index].stoichometric_vector[element_index] * species_density;
      else
        f_j += molecules[molecule_index].stoichometric_vector[species.index] * x;
    }


    f_j += molecules[molecule_index].mass_action_constant[grid_index];

    f_j = std::exp(f_j) * molecules[molecule_index].stoichometric_vector[species.index];

    f_i += f_j;
  }



  if (species.symbol == "e-")
    f_i = (- f_i - std::exp(x));
  else
    f_i = (species.abundance * h_density - f_i - std::exp(x));

  return f_i;
}



//Bisection method in one dimension
template <class double_type>
bool FastChem<double_type>::bisectionSolve(Element<double_type>& species, const double h_density, const unsigned int grid_index)
{
  //initial density interval
  std::vector<double_type> x(2, 0.0);

  x[1] = std::log(species.abundance * h_density);
  x[0] = std::log(element_density_minlimit);


  unsigned int nb_iterations = nb_max_bisection_iter;
  bool converged = false;


  for (unsigned int iter_step = 0; iter_step < nb_iterations; ++iter_step)
  {
    double_type x_n = std::log((std::exp(x[1]) - std::exp(x[0])) * 0.5 + std::exp(x[0]));

    double_type f_n = bisectionFunction(species, x_n, h_density, grid_index);

    if (f_n < 0)
      x[1] = x_n;
    else
      x[0] = x_n;


    //Convergence test. We need to be a little more accurate than the required accuracy.
    //Otherwise FastChem doesn't converge to the desired accuracy.
    if ( std::fabs(std::exp(x[0]) - std::exp(x[1]))/std::exp(x[1]) < accuracy * 1e-3  )
    {
      converged = true;
      break;
    }

  }


  species.number_density[grid_index] = std::exp(x[0]);


  if (!converged && verbose_level >= 3)
    std::cout << "Bisection iteration limit reached, result may not be optimal." << "\t" << x[0] << "\t" << x[1]
              << "\t" << std::fabs(std::exp(x[0]) - std::exp(x[1]))/std::exp(x[1]) << "\t" << accuracy * 1e-3  << "\n";


  return converged;
}


template class FastChem<double>;
template class FastChem<long double>;

}
