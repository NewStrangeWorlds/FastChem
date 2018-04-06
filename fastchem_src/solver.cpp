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
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//Selects the appropriate solver for each element
//See Sect. 2.4.2 and Eq. (2.32)
template <class double_type>
void FastChem<double_type>::calculateElementDensities(Element<double_type>& species, const double_type h_density, const double_type number_density_min,
                                                      const unsigned int grid_index, bool use_backup_solver)
{

  if (species.symbol == "e-") return; //electrons have their own, special solver

  //in case the usual FastChem iterations failed to converge, we switch to a backup
  if (use_backup_solver)
  {
    if (species.solver_order == 0)
      intertSol(species, h_density, number_density_min, grid_index);
    else
      backupSol(species, h_density, number_density_min, grid_index);

    return;
  }

  //selection of the solver for each element, see Eq. (2.32)
  switch (species.solver_order)
  {
    case 0 : intertSol(species, h_density, number_density_min, grid_index); break;
    case 1 : linSol(species, h_density, number_density_min, grid_index); break;
    case 2 : quadSol(species, h_density, number_density_min, grid_index); break;
    default : newtSol(species, h_density, number_density_min, grid_index);
  }

}



//This is the backup solver for the standard analytic FastChem solvers
//Currently configured to use a 1D Newton solver, that employs Eq. (2.34)
template <class double_type>
void FastChem<double_type>::backupSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min,
                                      const unsigned int grid_index)
{

  newtonSolveAlt(species, h_density, grid_index);

  checkN(species, h_density, grid_index);

}



//Solver for an element that is not part of other species
//See Eq. (2.32)
template <class double_type>
void FastChem<double_type>::intertSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index)
{

  species.number_density[grid_index] = species.abundance * h_density - number_density_min;

  checkN(species, h_density, grid_index);

}




//Calculates the number density of species, based on previously computed element densities
template <class double_type>
void FastChem<double_type>::calculateMoleculeDensities(const double_type& h_density, const unsigned int grid_index)
{

  //for performance reasons, we first calculate the logarithms of the element densities once here
  //and use those values later in the calculation of the molecule densities
  std::vector<double_type> log_element_densities(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    log_element_densities[i] = std::log(elements[i].number_density[grid_index]);



  for (size_t i=0; i<nb_molecules; ++i)
  {
    molecules[i].sum[grid_index] = 0.0;


    for (size_t l=0; l<molecules[i].element_indices.size(); ++l)
    {
      unsigned int index = molecules[i].element_indices[l];

      molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[index] * log_element_densities[index];
    }


    molecules[i].number_density[grid_index] = std::exp(molecules[i].sum[grid_index] + molecules[i].mass_action_constant[grid_index]);

    checkN(molecules[i], h_density, grid_index);
  }

}




//Determination of the correction factor n_j_min, that contains the number density of an element j contained in molecules
//with elements that are less abundant than j
//See Eq. (2.23)
template <class double_type>
void FastChem<double_type>::calculateMinorSpeciesDensities(std::vector<double_type>& number_density_min, const unsigned int grid_index)
{

  for (size_t j=0; j<nb_elements; ++j)
  {
    number_density_min[j] = 0.0;

    for (size_t i=0; i<elements[j].molecule_list.size(); ++i)
    {
      unsigned int index = elements[j].molecule_list[i];

      if (molecules[index].abundance < elements[j].abundance)
        number_density_min[j] += molecules[index].stoichometric_vector[j] * molecules[index].number_density[grid_index];
    }

  }


}




//Calculation of the electron density, based on charge conservation
//First tries to estimate the electron density via Eq. (2.35).
//In case that fails (electron density not sufficiently high enough), it switches to a 1D Nelder-Mead method.
//See Sect. 2.4.3 for details.
template <class double_type>
void FastChem<double_type>::calculateElectronDensities(const double_type& old_number_density, const double_type& h_density, const unsigned int grid_index)
{

  if (e_ == FASTCHEM_UNKNOWN_SPECIES) return;  //no electrons in the system


  elements[e_].number_density[grid_index] = 0.0;


  if (elements[e_].molecule_list.size() == 0) return; //no ions present


  double_type positive_ion_density = 0;
  double_type negative_ion_density = 0;


  for (size_t i=0; i<elements[e_].molecule_list.size(); ++i)
    if (molecules[elements[e_].molecule_list[i]].stoichometric_vector[e_] > 0)
      negative_ion_density += molecules[elements[e_].molecule_list[i]].stoichometric_vector[e_] * molecules[elements[e_].molecule_list[i]].number_density[grid_index];
    else
      positive_ion_density -= molecules[elements[e_].molecule_list[i]].stoichometric_vector[e_] * molecules[elements[e_].molecule_list[i]].number_density[grid_index];


  double_type electron_density = positive_ion_density - negative_ion_density;


  double_type delta = 0.9;


  if (electron_density > delta*positive_ion_density)
  {

    elements[e_].number_density[grid_index] = std::sqrt(electron_density * old_number_density);

  }
  else
  {
    //switching to Nelder-Mead method
    double_type init = std::log(std::fabs(old_number_density));

    nelderMeadSimplexSolve(elements[e_], init, h_density, grid_index);

  }


}



template class FastChem<double>;
template class FastChem<long double>;

}



