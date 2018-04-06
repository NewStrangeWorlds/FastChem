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


//Calls FastChem with p_<H> instead of p_gas, i.e. no pressure iterations are done; with diagnostic output
//Input: vector of temperatures, vector of p_<H>
//Output: number densities
//        n_<H>, mean molecular weights
//        element conservations
//        FastChem output flags (for each T,p individually)
//        number of pressure iterations & number of chemistry iterations (for each T,p individually)
//Function return: highest (i.e. most problematic) state of any of the FastChem calculations (max of output flag values)
template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const std::vector<double>& temperature, const std::vector<double>& hydrogen_pressure,
                                                  std::vector < std::vector<double> >& density_out,
                                                  std::vector<double>& mean_molecular_weight_out,
                                                  std::vector< std::vector<unsigned int> >& element_conserved_out,
                                                  std::vector<unsigned int>& fastchem_flags,
                                                  std::vector<unsigned int>& nb_chemistry_iterations_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  unsigned nb_grid_points = temperature.size();


  for (auto & i : species) i->number_density.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.sum.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.mass_action_constant.assign(nb_grid_points, 0);
  for (auto & i : elements) i.element_conserved.assign(nb_grid_points, false);


  element_conserved_out.resize(nb_grid_points);
  nb_chemistry_iterations_out.assign(nb_grid_points, 0);


  mean_molecular_weight_out.assign(nb_grid_points, 0);
  density_out.resize(nb_grid_points);



  std::vector<unsigned int> state(nb_grid_points, 0);


  for (unsigned int i=0; i<nb_grid_points; i++)
    state[i] = calcDensity(temperature[i], hydrogen_pressure[i],
                           i,
                           density_out[i], mean_molecular_weight_out[i],
                           element_conserved_out[i],
                           nb_chemistry_iterations_out[i]);


  fastchem_flags = state;


  return *std::max_element(state.begin(),state.end());
}



template class FastChem<double>;
template class FastChem<long double>;


}



