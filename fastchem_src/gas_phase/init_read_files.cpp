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
#include <fstream>
#include <sstream>
#include <cmath>

#include "gas_phase.h"

#include "../fastchem_constants.h"
#include "../species_struct.h"


namespace fastchem {


//Read the chemical elements file
template <class double_type>
bool GasPhase<double_type>::readSpeciesData(const std::string& file_path)
{ 
  std::fstream file(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open species data file " << file_path << "\n";
    return false;
  }


  molecules.reserve(10000);

  std::string line;

  //header
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string symbol, name, element_string, stoichiometric_coeff_string;

    input >> symbol; // >> name;

    std::string name_part;

    while (input >> name_part)
    {
      if (name_part == ":") break;

      if (name == "") 
        name = name_part;
      else
        name = name + " " + name_part; 
    }


    std::vector<std::string> species_elements;
    std::vector<int> stoichiometric_coeff;
    int charge = 0;


    while ( (input >> element_string >> stoichiometric_coeff_string) && (element_string != "#") )
    {
      species_elements.push_back(element_string);
      stoichiometric_coeff.push_back(std::stoi( stoichiometric_coeff_string ));

      if (species_elements.back() == "e-") charge = -stoichiometric_coeff.back();
    }


    std::string line;
    std::getline(file, line);

    std::istringstream ma_input(line);

    std::vector<double_type> mass_action_coeff;


    double ma_coefficient;

    while (ma_input >> ma_coefficient)
      mass_action_coeff.push_back(ma_coefficient);

    addMolecule(name, symbol, species_elements, stoichiometric_coeff, mass_action_coeff, charge);

    //blank separation line
    std::getline(file, line);
  }

  file.close();


  molecules.shrink_to_fit();

  nb_molecules = molecules.size();

  return true;
}



//Add a molecule to the system and update all of its elements
template <class double_type>
void GasPhase<double_type>::addMolecule(
  const std::string& name,
  const std::string& symbol,
  const std::vector<std::string>& species_elements,
  const std::vector<int>& stoichiometric_coeff,
  const std::vector<double_type>& mass_action_coeff,
  const int charge)
{
  Molecule<double_type> species;

  species.name = name;
  species.symbol = symbol;

  species.mass_action_coeff = mass_action_coeff;


  species.stoichiometric_vector.assign(nb_elements, 0);


  bool is_stoichiometry_complete = true;
  unsigned int nb_species_elements = 0;

  for (size_t i=0; i<species_elements.size(); ++i)
  {
    unsigned int index = element_data.elementIndex(species_elements[i]);

    if (index == FASTCHEM_UNKNOWN_SPECIES)
      is_stoichiometry_complete = false;
    else
    {
      species.stoichiometric_vector[index] = stoichiometric_coeff[i];
      species.element_indices.push_back(index);
    }

    nb_species_elements += stoichiometric_coeff[i];
  }


  if (is_stoichiometry_complete)
  {
    for (size_t j=0; j<nb_elements; ++j)
     species.sigma += species.stoichiometric_vector[j];

    species.sigma = 1 - species.sigma;
    species.charge = charge;
    species.phase = PhaseState::gas;

    for (auto & j : species.element_indices)
      species.weight += elements[j].weight * std::fabs(species.stoichiometric_vector[j]);

    molecules.push_back(species);

    //add the current molecule index to their respective elements
    for (auto & j : species.element_indices)
      elements[j].molecule_list.push_back(molecules.size()-1);
  }
  else 
    std::cout << "Stoichiometry of species " << symbol << " incomplete. Neglected!\n";
}



template class GasPhase<double>;
template class GasPhase<long double>;

} 
