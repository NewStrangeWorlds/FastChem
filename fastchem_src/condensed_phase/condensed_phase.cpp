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
  options.condensates_data_file = "input/condensate_all.dat";
  is_initialised = readCondensateData(options.condensates_data_file);

  if (options.verbose_level >= 4)
  {
    std::cout << "\nCondensate species list: \n";
    for (size_t i=0; i<condensates.size(); ++i)
    {
      std::cout << "  " << condensates[i].name << "\t" << condensates[i].symbol << "\n";
      
      std::cout << "    lnK coeff: ";
      for (size_t j=0; j<condensates[i].fit_coeff.size(); ++j)
        std::cout << condensates[i].fit_coeff[j] << "\t";
      std::cout << "\n";
      
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
  if (condensates_act.capacity() == 0)
    condensates_act.reserve(nb_condensates);

  for (auto & i : condensates)
  {
    if (i.log_activity >= 0 && i.is_calculated == false) 
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
        break;
      }
  }

}



template <class double_type>
void CondensedPhase<double_type>::selectJacobianCondensates(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& number_density_cond,
  const std::vector<double_type>& activity_corr,
  std::vector<unsigned int>& condensates_jac,
  std::vector<unsigned int>& condensates_rem)
{
  condensates_jac.resize(0);
  condensates_rem.resize(0);

  for (size_t i=0; i<condensates.size(); ++i)
  {
    if (condensates[i]->log_activity > -0.1 || options.cond_reduce_system_size == false)
      condensates_jac.push_back(i);
    else
      condensates_rem.push_back(i);
  }
}


template <class double_type>
void CondensedPhase<double_type>::selectJacobianCondensates2(
  const std::vector<Condensate<double_type>*>& condensates,
  const std::vector<double_type>& number_density_cond,
  const std::vector<double_type>& activity_corr,
  std::vector<unsigned int>& condensates_jac,
  std::vector<unsigned int>& condensates_rem,
  Eigen::MatrixXdt<double_type>& jacobian)
{
  condensates_jac.resize(0);
  condensates_rem.resize(0);

  const size_t nb_condensates = condensates.size();
  const size_t nb_elements = jacobian.rows() - nb_condensates;

  for (size_t i=0; i<condensates.size(); ++i)
  { 
    double element_min = 0;
    double element_max = 0;

    for (size_t j=0; j<nb_elements; ++j)
    {
      if (jacobian(j+nb_condensates, i) == 0) continue;

      for (size_t n=0; n<nb_elements; ++n)
      {
        if (jacobian(j+nb_condensates, n+nb_condensates) == 0) continue;

        if (jacobian(j+nb_condensates, n+nb_condensates) > element_max)
          element_max = jacobian(j+nb_condensates, n+nb_condensates);
        
        if (element_min == 0) element_min = element_max;
        else
          if (jacobian(j+nb_condensates, n+nb_condensates) < element_min)
            element_min = jacobian(j+nb_condensates, n+nb_condensates);
      }
    }

    if (activity_corr[i] < 1 && element_min < number_density_cond[i]*100)
      condensates_jac.push_back(i);
    else
      condensates_rem.push_back(i);
  }
}


template <class double_type>
double CondensedPhase<double_type>::totalElementDensity()
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



template class CondensedPhase<double>;
template class CondensedPhase<long double>;

}
