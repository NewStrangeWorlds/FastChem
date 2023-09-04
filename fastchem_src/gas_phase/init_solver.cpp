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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <limits> 
#include <cmath>

#include "gas_phase.h"


namespace fastchem {



template <class double_type>
void GasPhase<double_type>::determineSolverOrder()
{
  for (auto & i : elements)
    i.solver_order = determineSolverOrder(i);
}



//Determine the solver order for an element
//For the electrons this is the highest stage of ionisation
template <class double_type>
unsigned int GasPhase<double_type>::determineSolverOrder(const Element<double_type>& species)
{
  unsigned int solver_order = 0;

  if (species.symbol != "e-") //first the normal elements
  {
    for (auto & i : species.molecule_list)
      if (static_cast<unsigned int>(molecules[i].stoichiometric_vector[species.index]) > solver_order && molecules[i].abundance == species.abundance)
        solver_order = molecules[i].stoichiometric_vector[species.index];
  }
  else //then the electrons
  {
    for (auto & i : species.molecule_list)
      if (static_cast<unsigned int>(std::abs(molecules[i].stoichiometric_vector[species.index])) > solver_order)
        solver_order = static_cast<unsigned int>(std::abs(molecules[i].stoichiometric_vector[species.index])); 
  }

  return solver_order;
}



//Sort the elements according to their abundances
//elements will later be calculated in descending order
template <class double_type>
void GasPhase<double_type>::determineElementCalculationOrder()
{
  //make sure that there is a unique abundance for each element in case they are initially identical
  for (auto & i : element_data.elements_wo_e)
    for (auto & j : element_data.elements_wo_e)
    {
      if (i == j) continue;
  
      if (i->abundance == j->abundance)
        j->abundance += std::numeric_limits<double_type>::epsilon()*j->abundance;
    }

  //sort the element list according to their abundance
  std::sort(
    std::begin(element_data.elements_wo_e),
    std::end(element_data.elements_wo_e),
    [&] (Element<double_type>* a, Element<double_type>* b) {
      return a->abundance > b->abundance;});
  
  element_calculation_order.assign(element_data.elements_wo_e.size(), 0);
  //for (auto & i : element_data.elements_wo_e) std::cout << i->symbol << "\t" << std::setprecision(20) << i->abundance << "\n"; 
  for (size_t i=0; i<element_calculation_order.size(); ++i)
    element_calculation_order[i] = element_data.elements_wo_e[i]->index;
}



template class GasPhase<double>;
template class GasPhase<long double>;
}
