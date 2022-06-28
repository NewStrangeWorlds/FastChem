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


#include "fastchem.h"

#include <string>
#include <vector>
#include <algorithm>



namespace fastchem {


//Query for a species index (both, elements and molecules) with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getSpeciesIndex(const std::string symbol)
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



//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getSpeciesName(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->name;
  else
    return "";
}



//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->symbol;
  else
    return "";

}



//Query for a the molecular weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getSpeciesMolecularWeight(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->weight;
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



template class FastChem<double>;
template class FastChem<long double>;
}
