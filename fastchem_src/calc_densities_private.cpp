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


//Solve the chemistry for a single temperature and a single pressure
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
template <class double_type>
unsigned int FastChem<double_type>::calcDensity(const double temperature, const double pressure, const unsigned int grid_index,
                                                std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out,
                                                std::vector<unsigned int>& element_conserved_out,
                                                unsigned int& nb_pressure_iterations_out, unsigned int& nb_chemistry_iterations_out)
{
  for (auto & i : molecules)  i.calcMassActionConstant(temperature, grid_index);

  //this value will be fixed.
  double_type total_density = pressure/(CONST_K * temperature);

  //initial electron density
  unsigned int e_ = getElementIndex("e-");

  if (e_ != FASTCHEM_UNKNOWN_SPECIES)
    elements[e_].number_density[grid_index] = element_density_minlimit;


  double_type h_density = setInitialHDensity(total_density, grid_index);


  double_type density_iteration_lambda = 0.99; //initial value for lambda
  double_type density_iteration_error = 1.0;   //initial error


  bool fastchem_converged = false;
  bool pressure_converged = false;


  double_type muH = 0;
  double_type amu=1.66055e-24;
  for (size_t i=0; i<nb_elements; ++i)
    muH += elements[i].molecular_weight * chemical_elements[elements[i].element_index].abundance * amu;


  unsigned int nb_iterations = 0;
  unsigned int nb_fastchem_iterations = 0;


  for (nb_iterations=0; nb_iterations<nb_max_pressure_iter; ++nb_iterations)
  {

    fastchem_converged = solveFastchem(temperature, h_density, grid_index, nb_fastchem_iterations);
    pressure_converged = calcTotalHydrogenDensityAlt(temperature, pressure, grid_index,
                                                     h_density, muH, density_iteration_error);
    //pressure_converged = calcTotalHydrogenDensity(temperature, pressure, grid_index,
    //                                              h_density, density_iteration_lambda, density_iteration_error);

    if (std::isnan(h_density)) break;

    if (pressure_converged) break;
  }


  if (!pressure_converged && verbose_level >= 1) std::cout << "Pressure convergence problem in FastChem. :(\n";
  if (!fastchem_converged && verbose_level >= 1) std::cout << "FastChem convergence problem in FastChem. :(\n";


  //return output
  h_density_out = h_density;
  density_n_out.assign(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    density_n_out[i] = species[i]->number_density[grid_index];

  mean_molecular_weight_out = calcMeanMolecularWeight(total_density, grid_index);


  for (size_t i=0; i<nb_elements; i++)
    checkElementConservation(elements[i], h_density, grid_index);


  checkChargeConservation(grid_index);



  element_conserved_out.assign(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    element_conserved_out[i] = elements[i].element_conserved[grid_index];



  nb_pressure_iterations_out = nb_iterations;
  nb_chemistry_iterations_out = nb_fastchem_iterations;



  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!pressure_converged) return_state = FASTCHEM_NO_PRESSURE_CONVERGENCE;
  if (!fastchem_converged) return_state = FASTCHEM_NO_FASTCHEM_CONVERGENCE;

  if (!fastchem_converged && !pressure_converged) return_state = FASTCHEM_NO_CONVERGENCE;

  return return_state;
}




//Calls FastChem with p_<H> instead of p_gas, i.e. no pressure iterations are done
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by the public calcDensity function for p_<H>
template <class double_type>
unsigned int FastChem<double_type>::calcDensity(const double temperature, const double hydrogen_pressure, const unsigned int grid_index,
                                                std::vector<double>& density_n_out, double& mean_molecular_weight_out,
                                                std::vector<unsigned int>& element_conserved_out,
                                                unsigned int& nb_chemistry_iterations_out)
{
  for (auto & i : molecules)  i.calcMassActionConstant(temperature, grid_index);


  //initial electron density
  unsigned int e_ = getElementIndex("e-");

  if (e_ != FASTCHEM_UNKNOWN_SPECIES)
    elements[e_].number_density[grid_index] = element_density_minlimit;


  double_type h_density = hydrogen_pressure / (CONST_K * temperature);


  bool fastchem_converged = false;


  unsigned int nb_fastchem_iterations = 0;


  fastchem_converged = solveFastchem(temperature, h_density, grid_index, nb_fastchem_iterations);

  if (!fastchem_converged && verbose_level >= 1) std::cout << "FastChem convergence problem in FastChem. :(\n";


  //return output
  density_n_out.assign(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    density_n_out[i] = species[i]->number_density[grid_index];


  double_type total_density = 0;

  for (size_t i=0; i<nb_species; ++i)
    total_density += species[i]->number_density[grid_index];

  mean_molecular_weight_out = calcMeanMolecularWeight(total_density, grid_index);


  for (size_t i=0; i<nb_elements; i++)
    checkElementConservation(elements[i], h_density, grid_index);


  checkChargeConservation(grid_index);



  element_conserved_out.assign(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    element_conserved_out[i] = elements[i].element_conserved[grid_index];

  nb_chemistry_iterations_out = nb_fastchem_iterations;



  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!fastchem_converged) return_state = FASTCHEM_NO_FASTCHEM_CONVERGENCE;


  return return_state;
}



template class FastChem<double>;
template class FastChem<long double>;

}


