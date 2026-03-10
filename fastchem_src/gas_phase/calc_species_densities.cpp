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
void GasPhase::calculateElementDensities(
  Element& species,
  const double gas_density,
  bool use_backup_solver,
  double& n_major,
  const double log_gas_density)
{
  if (species.symbol == "e-") return; //electrons have their own, special solver

  species.number_density_maj = n_major * species.phi;
  
  if (species.fixed_by_condensation == false && species.epsilon > 0)
  {
    //in case the usual FastChem iterations failed to converge, we switch to a backup
    if (use_backup_solver)
    {
      if (species.solver_order == 0 && (species.minor_molecules.size() == 0 || species.molecule_list.size() == 0))
        solver.inertSol(species, elements, molecules, gas_density);
      else
        solver.backupSol(species, elements, molecules, gas_density);
    }
    else
    {
      //selection of the solver for each element, see Eq. (2.32)
      switch (species.solver_order)
      {
        case 0 : solver.inertSol(species, elements, molecules, gas_density); break;
        case 1 : solver.linSol(species, elements, molecules, gas_density); break;
        case 2 : solver.quadSol(species, elements, molecules, gas_density); break;
        default : solver.newtonSol(species, elements, molecules, gas_density, false);
      }
    }
  }

  if (species.epsilon == 0)
  {
    species.number_density = 0.0;
    species.log_number_density = static_cast<double>(LOG_DENSITY_FLOOR);
  }

  species.checkN(options.element_density_minlimit, log_gas_density);

  n_major += calculateMoleculeDensities(species, log_gas_density);
}




//Calculates the number density of species, based on previously computed element densities
double GasPhase::calculateMoleculeDensities(
  Element& species, const double log_gas_density)
{
  double n_major = 0.0;

  for (size_t ii=0; ii<species.major_molecules_inc.size(); ++ii)
  {
    const unsigned int i = species.major_molecules_inc[ii];
    double sum = 0.0;


    for (size_t ll=0; ll<molecules[i].element_indices.size(); ++ll)
    {
      const unsigned int l = molecules[i].element_indices[ll];

      sum += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    molecules[i].log_number_density = sum + molecules[i].mass_action_constant;
    molecules[i].number_density = safeExp(molecules[i].log_number_density);
    molecules[i].checkN(options.molecule_density_minlimit, log_gas_density);

    n_major += molecules[i].number_density * molecules[i].sigma;
  }

  return n_major;
}


//Calculates the number density of species, based on previously computed element densities
void GasPhase::updateMoleculeDensities()
{
  for (size_t i=0; i<nb_molecules; ++i)
  {
    double sum = 0.0;

    for (size_t ll=0; ll<molecules[i].element_indices.size(); ++ll)
    {
      const unsigned int l = molecules[i].element_indices[ll];

      sum += molecules[i].stoichiometric_vector[l] * elements[l].log_number_density;
    }

    molecules[i].log_number_density = sum + molecules[i].mass_action_constant;
    molecules[i].number_density = safeExp(molecules[i].log_number_density);
  }
}



double GasPhase::totalElementDensity()
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



}



