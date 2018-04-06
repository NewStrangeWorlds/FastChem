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

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>


namespace fastchem {


//This is the main FastChem iteration
template <class double_type>
bool FastChem<double_type>::solveFastchem(const double temperature_gas, const double_type h_density, const unsigned int grid_index, unsigned int& nb_iterations)
{
  std::vector<double_type> number_density_old(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    number_density_old[i] = species[i]->number_density[grid_index];


  //starting values for contribution of minor species
  std::vector<double_type> number_density_min(nb_elements, 0.0);


  unsigned int O_ = getSpeciesIndex("O");
  unsigned int C_ = getSpeciesIndex("C");
  unsigned int e_ = getSpeciesIndex("e-");


  if (elements[O_].abundance > elements[C_].abundance)
    number_density_min[O_] = elements[C_].abundance * h_density;
  else
    number_density_min[C_] = elements[O_].abundance * h_density;


  bool converged = false;
  unsigned int iter_step = 0;
  bool use_backup_solver = false;


  unsigned int max_iter = nb_max_fastchem_iter;

  for (iter_step=0; iter_step<max_iter; ++iter_step)
  {
    for (std::vector<unsigned int>::iterator it = element_calculation_order.begin(); it<element_calculation_order.end(); it++)
      calculateElementDensities(elements[*it], h_density, number_density_min[*it], grid_index, use_backup_solver);

    calculateMoleculeDensities(h_density, grid_index);

    calculateMinorSpeciesDensities(number_density_min, grid_index);

    if (e_ != FASTCHEM_UNKNOWN_SPECIES)  //only calculate electrons if they are present in the element list
      calculateElectronDensities(number_density_old[e_], h_density, grid_index);



    //check if n_j_min are small enough, if not use backup solver
    for (unsigned int i=0; i<nb_elements; i++)
      if (number_density_min[i] > elements[i].abundance * h_density)
      {
        use_backup_solver = true;

        if (verbose_level >= 4)
        std::cout << "Too large n_j_min. Switching to backup. Grid index: " << grid_index << "\t Iteration step: " << iter_step << "\n";

        break;
      }


    //convergence check
    if (iter_step > 0)
    {
      converged = true;

      for (size_t i=0; i<nb_species; ++i)
        if (std::fabs((species[i]->number_density[grid_index] - number_density_old[i])) > accuracy*number_density_old[i]
             && species[i]->number_density[grid_index]/h_density > 1.e-155)
        {
          converged = false;
          break;
        }
    }


    if (converged)
      break;


    //in case the standard FastChem iteration doesn't converge, switch to the backup solver
    if (iter_step == max_iter-1 && !converged && use_backup_solver == false)
    {
      if (verbose_level >= 4)
        std::cout << "Standard FastChem iteration failed. Switching to backup. Grid index " << grid_index << "\n";

      use_backup_solver = true;
      max_iter += nb_max_fastchem_iter;
    }


    for (size_t i=0; i<nb_species; ++i)
      number_density_old[i] = species[i]->number_density[grid_index];
  }


  if (!converged && verbose_level >= 3) std::cout << "Intermediate convergence problem in FastChem. :(\n";


  nb_iterations = iter_step;

  return converged;
}



template class FastChem<double>;
template class FastChem<long double>;

}
