/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2024 Daniel Kitzmann, Joachim Stock
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

#include "gas_phase.h"

#include "../fastchem_constants.h"
#include "../species_struct.h"


namespace fastchem {


//This is the main FastChem iteration for the gas phase
template <class double_type>
bool GasPhase<double_type>::calculate(
  const double temperature_gas, const double gas_density, unsigned int& nb_iterations)
{
  for (auto & i : elements) i.number_density_maj = 0.0;

  
  //starting values for contribution of minor species
  for (auto & i : elements) i.calcMinorSpeciesDensities(molecules);


  std::vector<double_type> number_density_old(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    number_density_old[i] = species[i]->number_density;


  bool converged = false;
  unsigned int iter_step = 0;
  bool use_backup_solver = options.chem_use_backup_solver;

  unsigned int max_iter = options.nb_max_fastchem_iter;

  for (iter_step=0; iter_step<max_iter; ++iter_step)
  {
    double_type n_maj = 0.0;

    //check if n_j_min are small enough, if not use backup solver
    for (auto & i : elements)
      if ( (i.number_density_min + i.number_density_maj > i.phi * gas_density) && use_backup_solver == false)
      {
        use_backup_solver = true;
       
        if (options.verbose_level >= 4)
          std::cout << "Too large n_min and n_maj for species " 
            << i.symbol << ". Switching to backup.  Iteration step: " 
            << iter_step << "\n";

        break;
      }

    //calculate the element densities in their respective order
    for (auto it = element_calculation_order.begin(); it<element_calculation_order.end(); it++)
      calculateElementDensities(elements[*it], gas_density, use_backup_solver, n_maj);


    //calculate the molecule densities
    for (auto & i : elements) i.calcMinorSpeciesDensities(molecules);


    //only calculate electrons if they are present
    if (e_ != FASTCHEM_UNKNOWN_SPECIES) 
      calculateElectronDensities(elements[e_], number_density_old[e_], gas_density);
 
    
    //check if n_j_min are small enough, if not use backup solver
    for (auto & i : elements)
      if ( (i.number_density_min + i.number_density_maj > i.phi * gas_density) && use_backup_solver == false)
      {
        use_backup_solver = true;
       
        if (options.verbose_level >= 4)
          std::cout << "Too large n_min and n_maj for species " 
            << i.symbol << ". Switching to backup.  Iteration step: " 
            << iter_step << "\n";

        break;
      }


    //convergence check
    if (iter_step > 0)
    {
      converged = true;
      
      for (size_t i=0; i<nb_species; ++i)
        if (std::fabs((species[i]->number_density - number_density_old[i])) > options.chem_accuracy*number_density_old[i]
             && species[i]->number_density/gas_density > 1.e-155)
        { 
          // std::cout << iter_step << "\t" << species[i]->symbol << "\t" 
          //           << std::fabs((species[i]->number_density - number_density_old[i]))/number_density_old[i] << "\t" 
          //           << options.chem_accuracy*number_density_old[i] << "\t" 
          //           << species[i]->number_density << "\t" 
          //           <<  number_density_old[i] << "\t" 
          //           << use_backup_solver << "\n";
          converged = false;
          break;
        }
    }


    //sanity check
    for (auto & e : elements)
    {
      if (std::isnan(e.number_density) || std::isinf(e.number_density))
      {
        if (options.verbose_level >= 4)
          std::cout << "Encountered NaN or Inf number density for element " 
                    << e.symbol 
                    << ". Stopping calculation.\n";

        nb_iterations = iter_step;

        return false;
      }
    }


    if (converged)
      break;


    //switch to multi-dimensional Newton solver if the system
    //hasn't converged yet
    if (iter_step > options.nb_switch_to_newton)
    {
      if (options.verbose_level >= 4)
        std::cout << "Standard FastChem iteration failed. Switching to multi-dimensional Newton. " << "\n";

      std::vector<Element<double_type>*> newton_elements;

      solver.selectNewtonElements(
        elements,
        molecules,
        number_density_old,
        gas_density,
        newton_elements);

      if (newton_elements.size() > 0)
      {
        double total_element_density = totalElementDensity();

        solver.newtonSolMult(
          newton_elements,
          elements,
          molecules,
          total_element_density);

        //update densities
        for (auto & e : elements) 
          calculateMoleculeDensities(e, gas_density);

        for (auto & e : elements) 
          e.calcMinorSpeciesDensities(molecules);
      }

    }


    //in case the standard FastChem iteration doesn't converge, switch to the backup solver
    if (iter_step == 390 && !converged && use_backup_solver == false)
    {
      if (options.verbose_level >= 4)
        std::cout << "Standard FastChem iteration failed. Switching to backup. " << "\n";

      use_backup_solver = true;
    }


    for (size_t i=0; i<nb_species; ++i)
      number_density_old[i] = species[i]->number_density;
  }
  

  nb_iterations = iter_step;
  
  return converged;
} 



template class GasPhase<double>;
template class GasPhase<long double>;

} 