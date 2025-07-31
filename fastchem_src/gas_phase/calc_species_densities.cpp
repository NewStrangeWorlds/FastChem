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


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

#include "gas_phase.h"

#include "solver.h"


namespace fastchem {


//Selects the appropriate solver for each element
//See Sect. 2.4.2 and Eq. (2.32)
template <class double_type>
void GasPhase<double_type>::calculateElementDensities(
  Element<double_type>& species,
  const double_type gas_density,
  bool use_backup_solver,
  double_type& n_major)
{
  if (species.symbol == "e-") return; //electrons have their own, special solver

  species.number_density_maj = n_major * species.phi;

  if (species.fixed_by_condensation == false && species.epsilon > 0)
  {
    //calculate the scaling factor if required
    if (options.use_scaling_factor)
      species.calcSolverScalingFactor(elements, molecules, options.additional_scaling_factor);
    else
      species.solver_scaling_factor = 0.0;


    //in case the usual FastChem iterations failed to converge, we switch to a backup
    if (use_backup_solver)
    {
      if (species.solver_order == 0 && (species.minor_molecules.size() == 0 || species.molecule_list.size() == 0))
        solver.intertSol(species, elements, molecules, gas_density);
      else
        solver.backupSol(species, elements, molecules, gas_density);
    }
    else
    {
      //selection of the solver for each element, see Eq. (2.32)
      switch (species.solver_order)
      {
        case 0 : solver.intertSol(species, elements, molecules, gas_density); break;
        case 1 : solver.linSol(species, elements, molecules, gas_density); break;
        case 2 : solver.quadSol(species, elements, molecules, gas_density); break;
        default : solver.newtonSol(species, elements, molecules, gas_density, false);
      }
    }
  }

  if (species.epsilon == 0)
    species.number_density = 0.0;

  species.checkN(options.element_density_minlimit, gas_density);

  n_major += calculateMoleculeDensities(species, gas_density);;
}




//Calculates the number density of species, based on previously computed element densities
template <class double_type>
double_type GasPhase<double_type>::calculateMoleculeDensities(
  Element<double_type>& species, const double_type gas_density)
{
  double_type n_major = 0.0;
  
  for (size_t ii=0; ii<species.major_molecules_inc.size(); ++ii)
  {
    const unsigned int i = species.major_molecules_inc[ii];
    double_type sum = 0.0;


    for (size_t ll=0; ll<molecules[i].element_indices.size(); ++ll)
    {
      const unsigned int l = molecules[i].element_indices[ll];

      sum += molecules[i].stoichiometric_vector[l] * std::log(elements[l].number_density);
    }
    
    molecules[i].number_density = std::exp(sum + molecules[i].mass_action_constant);
    molecules[i].checkN(options.molecule_density_minlimit, gas_density);
    
    n_major += molecules[i].number_density * molecules[i].sigma;
  }

  return n_major;
}



template <class double_type>
double GasPhase<double_type>::totalElementDensity()
{
  double n_tot = 0.0;

  //first we count the elements locked in molecules and ions
  for (size_t i=0; i<nb_molecules; ++i)
  {
    for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
    {
      const unsigned int element_index = molecules[i].element_indices[j];

      n_tot += molecules[i].number_density * molecules[i].stoichiometric_vector[element_index];
    }
  }


  //then we add the free atoms
  for (auto & i : elements) n_tot += i.number_density;

  return n_tot;
}



template class GasPhase<double>;
template class GasPhase<long double>;
}



