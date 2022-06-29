/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2022 Daniel Kitzmann, Joachim Stock
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

#include "condensed_phase.h"

#include "../fastchem_constants.h"
#include "../species_struct.h"


namespace fastchem {

template <class double_type>
CondensedPhase<double_type>::CondensedPhase(
  FastChemOptions<double_type>& options_,
  ElementData<double_type>& element_data_)
    : options(options_)
    , element_data(element_data_)
    , elements(element_data.elements)
    , solver(options)
{
  nb_elements = elements.size();

  std::string file_path = "input/condensate_test_small.dat";

  //is_initialised = readCondensateData(options.condensates_data_file);
  is_initialised = readCondensateData(file_path);

  if (is_initialised) init();
}



template <class double_type>
CondensedPhase<double_type>::CondensedPhase(
  const CondensedPhase &obj,
  FastChemOptions<double_type>& options_,
  ElementData<double_type>& element_data_)
    : condensates(obj.condensates)
    , nb_condensates(obj.nb_condensates)
    , nb_elements(obj.nb_elements)
    , is_initialised(obj.is_initialised)
    , options(options_)
    , element_data(element_data_)
    , elements(element_data.elements)
    , solver(options)
{
   
}



template <class double_type>
void CondensedPhase<double_type>::init()
{
  for (auto & i : condensates)
    std::cout << i.symbol << "\n";

  for (auto & i : condensates) 
    i.findReferenceElement(elements);
}



template <class double_type>
void CondensedPhase<double_type>::selectActiveCondensates(
  std::vector< Condensate<double_type>* >& condensates_act,
  std::vector< Element<double_type>* >& elements_cond)
{
  condensates_act.resize(0);
  condensates_act.reserve(nb_condensates);

  for (auto & i : condensates)
    if (i.log_activity != -10 && i.log_activity > 0) 
      condensates_act.push_back(&i);
    //if (i.ln_activity != -10) active_condensates.push_back(&i);

  //for (auto & i : condensates)
    //active_condensates.push_back(&i);
  
  condensates_act.shrink_to_fit();


  elements_cond.resize(0);
  elements_cond.reserve(nb_elements);

  for (auto & i : element_data.elements)
    for (auto & j : condensates_act)
      if (j->stoichiometric_vector[i.index] != 0)
      {
        elements_cond.push_back(&i);
        break;
      }

  elements_cond.shrink_to_fit();
}



template class CondensedPhase<double>;
template class CondensedPhase<long double>;

}
