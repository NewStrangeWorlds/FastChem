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


#include "fastchem.h"
#include "options.h"

#include <string>
#include <vector>
#include <algorithm>



namespace fastchem {


//Query for a species index (both, elements and molecules) with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getGasSpeciesIndex(const std::string symbol)
{
  auto it = std::find_if(
    gas_phase.species.begin(),
    gas_phase.species.end(),
    [&] (const ChemicalSpecies<double_type>* a) {
      return a->symbol == symbol;});

  if (it == gas_phase.species.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - gas_phase.species.begin();
}


//Query for a condensate species index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getElementIndex(const std::string symbol)
{
  auto it = std::find_if(
    element_data.elements.begin(),
    element_data.elements.end(),
    [&] (const Element<double_type> a) {
      return a.symbol == symbol;});

  if (it == element_data.elements.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - element_data.elements.begin();
}


//Query for a condensate species index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getCondSpeciesIndex(const std::string symbol)
{
  auto it = std::find_if(
    condensed_phase.condensates.begin(),
    condensed_phase.condensates.end(),
    [&] (const Condensate<double_type> a) {
      return a.symbol == symbol;});

  if (it == condensed_phase.condensates.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - condensed_phase.condensates.begin();
}



//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getGasSpeciesName(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->name;
  else
    return "";
}


//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getElementName(const unsigned int species_index)
{
  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].name;
  else
    return "";
}


//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getCondSpeciesName(const unsigned int species_index)
{
  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].name;
  else
    return "";
}




//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getGasSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->symbol;
  else
    return "";
}



//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getElementSymbol(const unsigned int species_index)
{

  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].symbol;
  else
    return "";
}



//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getCondSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].symbol;
  else
    return "";
}



//Query for a the molecular weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getGasSpeciesWeight(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->weight;
  else
    return 0.;
}


//Query for a the weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getElementWeight(const unsigned int species_index)
{
  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].weight;
  else
    return 0.;
}


//Query for a the weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getCondSpeciesWeight(const unsigned int species_index)
{
  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].weight;
  else
    return 0.;
}


//Get the element abundance for a specific element
template <class double_type>
double FastChem<double_type>::getElementAbundance(const unsigned int species_index)
{
  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].abundance;
  else
    return 0.;
}



//Get the element abundandes for all elements
template <class double_type>
std::vector<double> FastChem<double_type>::getElementAbundances()
{
  std::vector<double> abundances(element_data.nb_elements, 0.0);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    abundances[i] = element_data.elements[i].abundance;

  return abundances;
}


//Set the element abundances for all elements
template <class double_type>
void FastChem<double_type>::setElementAbundances(std::vector<double> abundances)
{
  element_data.setAbundances(abundances);

  gas_phase.reInitialise();
}


//Set an internal FastChem parameter (for double values)
template <class double_type>
void FastChem<double_type>::setParameter(const std::string& parameter, const double_type value)
{
  auto param = options.resolveParameter(parameter);

  switch (param)
  {
    case ParameterFloat::cond_tau:
      options.cond_tau = value;
      break;
    
    case ParameterFloat::cond_limit_change:
      options.cond_iter_change_limit = value;
      break;

    case ParameterFloat::chem_accuracy:
      options.chem_accuracy = value;
      break;

    case ParameterFloat::element_conserve_accuracy:
      options.element_conserve_accuracy = value;
      break;

    case ParameterFloat::cond_accuracy:
      options.cond_accuracy = value;
      break;

    case ParameterFloat::newton_err:
      options.newton_err = value;
      break;

    case ParameterFloat::element_minlimit:
      options.element_density_minlimit = std::pow(10.0, value);
      break;

    case ParameterFloat::molecule_minlimit:
      options.molecule_density_minlimit = std::pow(10.0, value);
      break;

    case ParameterFloat::logK_limit:
      options.logK_limit = value;
      break;

    case ParameterFloat::additional_scaling_factor:
      options.additional_scaling_factor = value;
      break;
  
    default:
      std::cout << "Unknown parameter \"" << parameter << "\"  with a floatint-point value!\n";
      break;
  }

}


//Set an internal FastChem parameter (for boolean values)
template <class double_type>
void FastChem<double_type>::setParameter(const std::string& parameter, const bool value)
{
  auto param = options.resolveParameterBool(parameter);

  switch (param)
  {
    case ParameterBool::cond_solve_full_system:
      options.cond_solve_full_matrix = value;
      break;

    case ParameterBool::cond_reduce_system_size:
      options.cond_reduce_system_size = value;
      break;

    case ParameterBool::cond_use_full_pivot:
      options.cond_use_full_pivot = value;
      break;

    case ParameterBool::cond_use_svd:
      options.cond_use_svd = value;
      break;

    case ParameterBool::use_scaling_factor:
      options.use_scaling_factor = value;
      break;

     case ParameterBool::cond_use_data_validity_limits:
      options.cond_use_data_validity_limits = value;
      break;
  
    default:
      std::cout << "Unknown parameter \"" << parameter << "\"  with a boolean value!\n";
      break;
  }
}



//Set an internal FastChem parameter (for integer values)
template <class double_type>
void FastChem<double_type>::setParameter(const std::string& parameter, const unsigned int value)
{
  auto param = options.resolveParameterInt(parameter);

  switch (param)
  {
    case ParameterInt::nb_max_bisection_iter:
      options.nb_max_bisection_iter = value;
      break;

    case ParameterInt::nb_max_comb_iter:
      options.nb_chem_cond_iter = value;
      break;

    case ParameterInt::nb_max_cond_iter:
      options.nb_max_cond_iter = value;
      break;

    case ParameterInt::nb_max_fastchem_iter:
      options.nb_max_fastchem_iter = value;
      break;

    case ParameterInt::nb_max_neldermead_iter:
      options.nb_max_neldermead_iter = value;
      break;

    case ParameterInt::nb_max_newton_iter:
      options.nb_max_newton_iter = value;
      break;
  
    default:
      std::cout << "Unknown parameter \"" << parameter << "\"  with an integer value!\n";
      break;
  }
}



template class FastChem<double>;
template class FastChem<long double>;
}
