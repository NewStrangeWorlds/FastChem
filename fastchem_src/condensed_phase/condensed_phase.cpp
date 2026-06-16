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
#include <cmath>
#include <algorithm>

#include "condensed_phase.h"

#include "../fastchem_constants.h"
#include "../species_struct.h"


namespace fastchem {

CondensedPhase::CondensedPhase(
  FastChemOptions& options_,
  ElementData& element_data_)
    : options(options_)
    , element_data(element_data_)
    , elements(element_data.elements)
    , solver(options)
{
  nb_elements = elements.size();

  if (options.condensates_data_file != "none")
    is_initialised = readCondensateData(options.condensates_data_file);
  else
  {
    is_initialised = false;
    nb_condensates = 0;

    return;
  }


  if (options.verbose_level >= 4)
  {
    std::cout << "\nCondensate species list: \n";
    for (size_t i=0; i<condensates.size(); ++i)
    {
      std::cout << "  " << condensates[i].name << "\t" << condensates[i].symbol << "\n";
      
      std::cout << "    lnK coeff: ";
      for (size_t j=0; j<condensates[i].fit_coeff.size(); ++j)
      {
        for (auto& c : condensates[i].fit_coeff[j])
          std::cout << c << "\t";

       std::cout << "\n";
      }

      std::cout << "    stoichiometry: ";
      for (size_t j=0; j<condensates[i].stoichiometric_vector.size(); ++j)
        std::cout << condensates[i].stoichiometric_vector[j] << " ";
      std::cout << "\n";

      std::cout << "    elements: ";
      for (size_t j=0; j<condensates[i].element_indices.size(); ++j)
        std::cout << elements[condensates[i].element_indices[j]].symbol << ", index: " << condensates[i].element_indices[j] << "; ";
      std::cout << "\n";

      std::cout << "    phase: " << phase_state_string[condensates[i].phase] << "\n";
    }
  }

  if (is_initialised) init();
}



CondensedPhase::CondensedPhase(
  const CondensedPhase &obj,
  FastChemOptions& options_,
  ElementData& element_data_)
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



void CondensedPhase::init()
{
  for (auto & i : condensates) 
    i.findReferenceElement(elements);
}



//A condensate is a "ghost" if its maximum possible density (set by its most
//depleted constituent element) is negligible compared to the total element
//density. This happens in rainout calculations once an element has essentially
//fully condensed out at higher temperatures: its remaining abundance is so small
//that the condensate cannot meaningfully form, yet its activity can still be
//strongly supersaturated. Feeding such condensates into the Newton solver creates
//a degenerate, ill-conditioned system that produces NaNs.
bool CondensedPhase::isGhostCondensate(
  const Condensate& condensate,
  const double total_element_density) const
{
  if (!std::isfinite(total_element_density) || total_element_density <= 0)
    return false;

  // for (auto & e : condensate.element_indices)
  //   std::cout << condensate.symbol << "  " << elements[e].symbol << "  " << elements[e].phi << "\t" << elements[e].degree_of_condensation << "\n";

  // for (auto & e : condensate.element_indices)
  //   if (elements[e].phi < 1e-10) return true;

  return condensate.max_number_density
    < options.cond_limiting_density_ratio * total_element_density;
}



void CondensedPhase::selectActiveCondensates(
  std::vector< Condensate* >& condensates_act,
  std::vector< Element* >& elements_cond,
  const double total_element_density)
{
  if (condensates_act.capacity() == 0)
    condensates_act.reserve(nb_condensates);

  for (auto & i : condensates)
  {
    if (i.log_activity >= 0 && i.is_calculated == false
        && !isGhostCondensate(i, total_element_density))
    {
      condensates_act.push_back(&i);
      i.is_calculated = true;

      i.number_density = 0.0;
      i.activity_correction = 0.0;
    }
  }


  if (elements_cond.capacity() == 0)
    elements_cond.reserve(nb_elements);

  for (auto & i : element_data.elements)
  {
    if (std::find(elements_cond.begin(), elements_cond.end(), &i) != elements_cond.end())
      continue;

    for (auto & j : condensates_act)
      if (j->stoichiometric_vector[i.index] != 0)
      {
        elements_cond.push_back(&i);
        elements_cond.back()->fixed_by_condensation = true;
        break;
      }
  }

  for (auto & e : elements_cond)
    e->fixed_by_condensation = true;
}



void CondensedPhase::selectJacobianCondensates(
  const std::vector<Condensate*>& condensates_act,
  const std::vector<double>& number_density_cond,
  const std::vector<double>& activity_corr,
  std::vector<unsigned int>& condensates_jac,
  std::vector<unsigned int>& condensates_rem)
{
  condensates_jac.resize(0);
  condensates_rem.resize(0);

  for (size_t i=0; i<condensates_act.size(); ++i)
  {
    if (condensates_act[i]->log_activity > -0.1 || options.cond_reduce_system_size == false)
      condensates_jac.push_back(i);
    else
      condensates_rem.push_back(i);
  }
}



double CondensedPhase::totalElementDensity()
{
  double n_tot = 0.0;

  for (size_t i=0; i<nb_condensates; ++i)
  {
    for (size_t j=0; j<condensates[i].element_indices.size(); ++j)
    {
      const unsigned int element_index = condensates[i].element_indices[j];

      n_tot += condensates[i].number_density * condensates[i].stoichiometric_vector[element_index];
    }
  }

  return n_tot;
}




}
