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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <cmath>



namespace fastchem {


//Add a new atom to the system
template <class double_type>
void FastChem<double_type>::addAtom(std::string symbol)
{
  Element<double_type> species;

  species.symbol = symbol;
  species.element_data_index = getChemicalElementIndex(symbol);


  if (species.element_data_index == FASTCHEM_UNKNOWN_SPECIES)
    std::cout << "Element " << symbol << " from element abundance file not found in element data file. Neglected!\n";
  else
  {
    species.name = chemical_element_data[species.element_data_index].name;
    species.weight = chemical_element_data[species.element_data_index].atomic_weight;
    species.abundance = chemical_element_data[species.element_data_index].abundance;
    elements.push_back(species);

    elements.back().index = elements.size()-1;
  }

}



//Set an element abundance that was found in the input file
template <class double_type>
void FastChem<double_type>::setElementAbundance(const std::string symbol, const double abundance)
{
  unsigned int index = getChemicalElementIndex(symbol);

  if (index == FASTCHEM_UNKNOWN_SPECIES)
    std::cout << "Element " << symbol << " for setting abundances not found. Neglected!\n";
  else
    chemical_element_data[index].abundance = abundance;


  if (symbol == "e-") chemical_element_data[index].abundance = 0.0;
}



//Add a molecule to the system and update all of its elements
template <class double_type>
void FastChem<double_type>::addMolecule(const std::string name, const std::string symbol,
                                        const std::vector<std::string> species_elements, const std::vector<int> stoichiometric_coeff,
                                        const std::vector<double_type> mass_action_coeff, const std::string mass_action_coeff_file,
                                        const int charge)
{
  Molecule<double_type> species;

  species.name = name;
  species.symbol = symbol;

  species.mass_action_coeff = mass_action_coeff;

  bool read_mass_action_file = true;

  if (mass_action_coeff_file.empty() == false)
    read_mass_action_file = species.readMassActionConstants(mass_action_coeff_file);


  species.stoichiometric_vector.assign(nb_elements, 0);


  bool is_stoichiometry_complete = true;
  unsigned int nb_species_elements = 0;

  for (size_t i=0; i<species_elements.size(); ++i)
  {
    unsigned int index = getElementIndex(species_elements[i]);

    if (index == FASTCHEM_UNKNOWN_SPECIES)
      is_stoichiometry_complete = false;
    else
    {
      species.stoichiometric_vector[index] = stoichiometric_coeff[i];
      species.element_indices.push_back(index);
    }

    nb_species_elements += stoichiometric_coeff[i];
  }


  if (is_stoichiometry_complete && read_mass_action_file)
  {
    for (size_t j=0; j<nb_elements; ++j)
     species.sigma += species.stoichiometric_vector[j];

    species.sigma = 1 - species.sigma;
    species.charge = charge;
    
    for (auto & j : species.element_indices)
      species.weight += elements[j].weight * std::fabs(species.stoichiometric_vector[j]);

    molecules.push_back(species);

    //add the current molecule index to their respective elements
    for (auto & j : species.element_indices)
      elements[j].molecule_list.push_back(molecules.size()-1);
  }
  else 
    std::cout << "Stoichiometry or mass action data of species " << symbol << " incomplete. Neglected!\n";

}



template class FastChem<double>;
template class FastChem<long double>;
}
