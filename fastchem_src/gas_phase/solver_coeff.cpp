/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2021 Daniel Kitzmann, Joachim Stock
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


#include "solver.h"
#include "../species_struct.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


template <class double_type>
double_type GasPhaseSolver<double_type>::A0Coeff(
  const Element<double_type>& species, const double_type gas_density)
{
  double_type A0 = 0.0;

  if (options.use_scaling_factor)
    A0 = std::exp(-species.solver_scaling_factor) * (species.number_density_maj + species.number_density_min - gas_density * species.phi);
  else
    //calculation of coefficient A_j0, see Eq. (2.27)
    A0 = species.number_density_maj + species.number_density_min - gas_density * species.phi;

  return A0;
}


template <class double_type>
double_type GasPhaseSolver<double_type>::A1Coeff(
  const Element<double_type>& species,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules)
{
  //calculation of coefficient A_1, see Eq. (2.28)
  double_type A1 = 0.0;

  for(auto & i : species.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[species.index] == 1 && molecules[i].abundance == species.abundance)
    {
      double_type sum = 0;
      
      for (auto & j : molecules[i].element_indices)
      {
        if (j != species.index && molecules[i].stoichiometric_vector[j] != 0)
          sum += molecules[i].stoichiometric_vector[j] * std::log(elements[j].number_density);
      }

      const double_type kappa = 1.0 + species.phi * molecules[i].sigma;
      
      A1 += std::exp(molecules[i].mass_action_constant + sum - species.solver_scaling_factor) * kappa;
    }
  }


  if (options.use_scaling_factor)
    A1 += std::exp(-species.solver_scaling_factor);
  else
    A1 += 1.0;


  return A1;
}



template <class double_type>
double_type GasPhaseSolver<double_type>::A2Coeff(
  const Element<double_type>& species,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules)
{
  double_type A2 = 0.0;

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[species.index] == 2 && molecules[i].abundance == species.abundance)
    {
      double_type sum = 0;
      
      for (auto & j : molecules[i].element_indices)
      {
        if (j != species.index && molecules[i].stoichiometric_vector[j] != 0)
          sum += molecules[i].stoichiometric_vector[j] * std::log(elements[j].number_density);
      }

      const double_type kappa = 2.0 + species.phi * molecules[i].sigma; 
      
      A2 += std::exp(molecules[i].mass_action_constant + sum - species.solver_scaling_factor) * kappa;
    }
  }
  
  
  return A2;
}



template <class double_type>
double_type GasPhaseSolver<double_type>::AmCoeff(
  const Element<double_type>& species,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules, 
  const unsigned int order)
{
  double_type Am = 0.0;

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[species.index] == int(order) && molecules[i].abundance == species.abundance)
    {
      double_type sum = 0;

      for (auto & j : molecules[i].element_indices)
      {
        if (j != species.index && molecules[i].stoichiometric_vector[j] != 0)
          sum += molecules[i].stoichiometric_vector[j] * std::log(elements[j].number_density);
      }

      const double_type kappa = order + species.phi * molecules[i].sigma; 
      
      Am += std::exp(molecules[i].mass_action_constant + sum - species.solver_scaling_factor) * kappa;
    }
  }
  

  if (order == 1)
  {
    if (options.use_scaling_factor)
      Am += std::exp(-species.solver_scaling_factor);
    else
      Am += 1.0;
  }

  return Am;
}



//Alternative description of the Am coefficients
//Takes the full law of mass action into account, i.e. doesn't stop at minor species as the regular calculation
template <class double_type>
double_type GasPhaseSolver<double_type>::AmCoeffAlt(
  const Element<double_type>& species,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules, 
  const unsigned int order)
{
  double_type Am = 0.0;

  for (auto & i : species.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[species.index] == int(order))
    {
      double_type sum = 0;
      
      for (auto & j : molecules[i].element_indices)
      {
        if (j != species.index && molecules[i].stoichiometric_vector[j] != 0)
          sum += molecules[i].stoichiometric_vector[j] * std::log(elements[j].number_density);
      }

      const double_type kappa = order + species.phi * molecules[i].sigma; 
      
      Am += std::exp(molecules[i].mass_action_constant + sum - species.solver_scaling_factor) * kappa;
    }
  }


  if (order == 1) 
  {
    if (options.use_scaling_factor)
      Am += std::exp(-species.solver_scaling_factor);
    else
      Am += 1.0;
  }

  return Am;
}



//Alternative description of the Am coefficients
//Takes the full law of mass action into account, i.e. doesn't stop at minor species as the regular calculation
template <class double_type>
double_type GasPhaseSolver<double_type>::AmCoeffElectron(
  const Element<double_type>& electron,
  const std::vector< Element<double_type> >& elements,
  const std::vector< Molecule<double_type> >& molecules, 
  const int order)
{
  double_type Am = 0.0;

  for (auto & i : electron.molecule_list)
  {
    if (molecules[i].stoichiometric_vector[electron.index] == order)
    {
      double_type sum = 0;
      
      for (auto & j : molecules[i].element_indices)
      {
        if (j != electron.index && molecules[i].stoichiometric_vector[j] != 0)
          sum += molecules[i].stoichiometric_vector[j] * std::log(elements[j].number_density);
      }
      
      Am += std::exp(molecules[i].mass_action_constant + sum) * order;
    }
  }


  return Am;
}


template class GasPhaseSolver<double>;
template class GasPhaseSolver<long double>;
}



