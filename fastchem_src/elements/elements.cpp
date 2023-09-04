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
#include <algorithm>

#include "elements.h"

#include "../fastchem_constants.h"
#include "chemical_element_data.h"
#include "../species_struct.h"


namespace fastchem {


template <class double_type>
ElementData<double_type>::ElementData(const std::string& abundance_file, const std::string& chemical_element_file)
{
  if (readElementList(chemical_element_file) == false || readElementAbundances(abundance_file) == false)
    is_initialised = false;
  else
    is_initialised = true;

  //setup the element list without the electron
  elements_wo_e.resize(0);
  elements_wo_e.reserve(nb_elements);

  for (auto & i : elements)
    if (i.symbol != "e-") elements_wo_e.push_back(&i);

  elements_wo_e.shrink_to_fit();

  e_ = elementIndex("e-");
}



template <class double_type>
ElementData<double_type>::ElementData(const ElementData &obj)
  : elements(obj.elements)
  , nb_elements(obj.nb_elements)
  , e_(obj.e_)
  , is_initialised(obj.is_initialised)
  , chemical_element_data(obj.chemical_element_data)
  , nb_chemical_element_data(obj.nb_chemical_element_data)
{
  elements_wo_e.reserve(nb_elements);

  for (auto & i : elements)
    if (i.symbol != "e-") elements_wo_e.push_back(&i);
}



template <class double_type>
void ElementData<double_type>::init(double_type initial_density)
{
  for (auto & i : elements)
  {
    i.number_density = initial_density;
    i.degree_of_condensation = 0;
    i.phi = i.epsilon;
    i.fixed_by_condensation = false;
  }

}



//Add a new atom to the system
template <class double_type>
void ElementData<double_type>::add(const std::string& symbol)
{
  Element<double_type> element;

  element.symbol = symbol;
  element.element_data_index = chemicalElementIndex(symbol);


  if (element.element_data_index == FASTCHEM_UNKNOWN_SPECIES)
    std::cout << "Element " << symbol << " from element abundance file not found in element data file. Neglected!\n";
  else
  {
    element.name = chemical_element_data[element.element_data_index].name;
    element.weight = chemical_element_data[element.element_data_index].atomic_weight;
    element.abundance = chemical_element_data[element.element_data_index].abundance;

    element.phase = PhaseState::gas;

    elements.push_back(element);
    elements.back().index = elements.size()-1;
  }
}



//Set the element abundances for all elements
template <class double_type>
void ElementData<double_type>::setAbundances(const std::vector<double>& abundances)
{
  if (abundances.size() != nb_elements)
  {
    std::cout << "Setting element abundances with an incorrect vector size\n";

    return;
  }

  for (size_t i=0; i<nb_elements; ++i)
  { 
    if (i == e_) continue; //the abundance of the electron remains at zero

    chemical_element_data[elements[i].element_data_index].abundance = abundances[i];
    elements[i].abundance = abundances[i];
  }
}



template <class double_type>
void ElementData<double_type>::setRelativeAbundances()
{
  double_type phi_tot = 0;

  for (auto & i : elements)
  {
    i.phi = i.epsilon * (1.0 - i.degree_of_condensation);
    phi_tot += i.phi;
  }

  for (auto & i : elements)
    i.phi /= phi_tot;
}



//Set an element abundance that was found in the input file
template <class double_type>
void ElementData<double_type>::setAbundance(
  const std::string& symbol, const double abundance)
{
  unsigned int index = chemicalElementIndex(symbol);

  if (index == FASTCHEM_UNKNOWN_SPECIES)
    std::cout << "Element " << symbol << " for setting abundances not found. Neglected!\n";
  else
    chemical_element_data[index].abundance = abundance;

  if (symbol == "e-") chemical_element_data[index].abundance = 0.0;
}



//Query for a element's index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the element does not exist
template <class double_type>
unsigned int ElementData<double_type>::elementIndex(const std::string& symbol)
{
  auto it = std::find_if(
    elements.begin(),
    elements.end(),
    [&](const Element<double_type>& a) { return a.symbol == symbol;});

  if (it == elements.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - elements.begin();
}



//Query for a basic chemical element index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the element does not exist
template <class double_type>
unsigned int ElementData<double_type>::chemicalElementIndex(const std::string& symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_chemical_element_data; ++i)
    if (symbol == chemical_element_data[i].symbol)
    {
      index = i;
      break;
    }

  return index;
}



template class ElementData<double>;
template class ElementData<long double>;
} 
