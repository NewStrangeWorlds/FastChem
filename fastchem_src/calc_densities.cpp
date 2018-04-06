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


//Calls FastChem with a vector of temperatures and pressures
//Returns no detailed diagnostic output
//Input: vector of temperatures, vector of pressures
//Output: number densities
//        n_<H>, mean molecular weights
//Function return: highest (i.e. most problematic) state of any of the FastChem calculations
template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const std::vector<double>& temperature, const std::vector<double>& pressure,
                                                  std::vector < std::vector<double> >& density_out,
                                                  std::vector<double>& h_density_out, std::vector<double>& mean_molecular_weight_out)
{
  std::vector< std::vector<unsigned int> > element_conservation_out;
  std::vector<unsigned int> pressure_iteration_steps_out;
  std::vector<unsigned int> chemistry_iteration_steps_out;
  std::vector<unsigned int> fastchem_flags;


  return calcDensities(temperature, pressure,
                       density_out,
                       h_density_out, mean_molecular_weight_out,
                       element_conservation_out,
                       fastchem_flags,
                       pressure_iteration_steps_out, chemistry_iteration_steps_out);
}




//Calls FastChem with an array of temperatures and pressures, returns also diagnostic output
//Input: vector of temperatures, vector of pressures
//Output: number densities
//        n_<H>, mean molecular weights
//        element conservations
//        FastChem output flags (for each T,p individually)
//        number of pressure iterations & number of chemistry iterations (for each T,p individually)
//Function return: highest (i.e. most problematic) state of any of the FastChem calculations (max of output flag values)
template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const std::vector<double>& temperature, const std::vector<double>& pressure,
                                                  std::vector < std::vector<double> >& density_out,
                                                  std::vector<double>& h_density_out, std::vector<double>& mean_molecular_weight_out,
                                                  std::vector< std::vector<unsigned int> >& element_conserved_out,
                                                  std::vector<unsigned int>& fastchem_flags,
                                                  std::vector<unsigned int>& nb_pressure_iterations_out, std::vector<unsigned int>& nb_chemistry_iterations_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  unsigned nb_grid_points = temperature.size();


  for (auto & i : species) i->number_density.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.sum.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.mass_action_constant.assign(nb_grid_points, 0);
  for (auto & i : elements) i.element_conserved.assign(nb_grid_points, false);


  element_conserved_out.resize(nb_grid_points);
  nb_pressure_iterations_out.assign(nb_grid_points, 0);
  nb_chemistry_iterations_out.assign(nb_grid_points, 0);


  h_density_out.assign(nb_grid_points, 0);
  mean_molecular_weight_out.assign(nb_grid_points, 0);
  density_out.resize(nb_grid_points);



  std::vector<unsigned int> state(nb_grid_points, 0);


  for (unsigned int i=0; i<nb_grid_points; i++)
    state[i] = calcDensity(temperature[i], pressure[i],
                           i,
                           density_out[i], h_density_out[i], mean_molecular_weight_out[i],
                           element_conserved_out[i],
                           nb_pressure_iterations_out[i], nb_chemistry_iterations_out[i]);


  fastchem_flags = state;


  return *std::max_element(state.begin(),state.end());
}




//Calls FastChem for a single temperature and pressure
//Input: temperature, pressure
//Output: number densities, n_<H>, mean molecular weight
//        element conservation diagnostic
//        number of pressure iterations, number of chemistry iterations
//Function return: final state of FastChem calculation
template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const double temperature, const double pressure,
                                                  std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out,
                                                  std::vector<unsigned int>& element_conserved_out,
                                                  unsigned int& nb_pressure_iterations_out, unsigned int& nb_chemistry_iterations_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  std::vector<double> temperature_vector(1, temperature);
  std::vector<double> pressure_vector(1, pressure);

  std::vector<double> h_density_vector, mean_molecular_weight_vector;
  std::vector< std::vector<double> > density_tensor;

  std::vector< std::vector<unsigned int> > element_conserved_tensor;
  std::vector<unsigned int> nb_pressure_iterations_vector;
  std::vector<unsigned int> nb_chemistry_iterations_vector;
  std::vector<unsigned int> fastchem_flag_vector;



  unsigned int state = calcDensities(temperature_vector, pressure_vector,
                                     density_tensor, h_density_vector, mean_molecular_weight_vector,
                                     element_conserved_tensor,
                                     fastchem_flag_vector, nb_pressure_iterations_vector, nb_chemistry_iterations_vector);

  density_n_out = density_tensor[0];
  h_density_out = h_density_vector[0];
  mean_molecular_weight_out = mean_molecular_weight_vector[0];

  element_conserved_out = element_conserved_tensor[0];
  nb_pressure_iterations_out = nb_pressure_iterations_vector[0];
  nb_chemistry_iterations_out = nb_chemistry_iterations_vector[0];


  return state;
}




//Calls FastChem with a single temperature and pressure, does not return diagnostic output, except for the final state
//Input: temperature, pressure
//Output: number densities, n_<H>, mean molecular weight
//Function return: final state of FastChem calculation
template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const double temperature, const double pressure,
                                                  std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  std::vector<unsigned int> element_conserved;
  unsigned int nb_pressure_iterations, nb_chemistry_iterations;


  unsigned int state = calcDensities(temperature, pressure,
                                     density_n_out, h_density_out, mean_molecular_weight_out,
                                     element_conserved, nb_pressure_iterations, nb_chemistry_iterations);



  return state;
}



template class FastChem<double>;
template class FastChem<long double>;


}


