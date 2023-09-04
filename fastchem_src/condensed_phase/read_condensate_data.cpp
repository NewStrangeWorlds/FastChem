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


#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "condensed_phase.h"

#include "../species_struct.h"


namespace fastchem {


template <class double_type>
bool CondensedPhase<double_type>::readCondensateData(const std::string& species_data_file)
{
  std::fstream file(species_data_file.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open condensate data file " << species_data_file << "\n";
    return false;
  }

  std::string line;

  //header
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string symbol, name, element_string, stoichometric_coeff_string;

    input >> symbol;

    std::string name_part;

    while (input >> name_part)
    {
      if (name_part == ":") break;

      if (name == "") 
        name = name_part;
      else
        name = name + " " + name_part; 
    }


    std::vector<std::string> elements;
    std::vector<int> stoichometric_coeff;

    while ( (input >> element_string >> stoichometric_coeff_string) && (element_string != "#") )
    {
      elements.push_back(element_string);
      stoichometric_coeff.push_back(std::stoi( stoichometric_coeff_string ));
    }

    //read in the phase and density
    std::string line;
    std::getline(file, line);
    std::istringstream basic_data_input(line);

    std::string phase = "";
    basic_data_input >> phase;


    //read in the phase temperature limits
    std::getline(file, line);
    std::istringstream fit_limits_input(line);
    std::vector<double> fit_coeff_limits;

    double fit_coeff_limit;

    while (fit_limits_input >> fit_coeff_limit)
      fit_coeff_limits.push_back(fit_coeff_limit);


    std::vector<std::vector<double_type>> fit_coeff(fit_coeff_limits.size(), std::vector<double_type>{});

    for (size_t i=0; i<fit_coeff_limits.size(); ++i)
    {
      std::getline(file, line);
      std::istringstream fit_input(line);

      double fit_coefficient;
      
      while (fit_input >> fit_coefficient)
        fit_coeff[i].push_back(fit_coefficient);
    }

    addCondensate(name, symbol, elements, stoichometric_coeff, phase, fit_coeff_limits, fit_coeff);

    //blank separation line
    std::getline(file, line);
  }

  file.close();

  nb_condensates = condensates.size();

  return true;
}



//Add a condensate to the system and update all of its elements
template <class double_type>
void CondensedPhase<double_type>::addCondensate(
  const std::string name,
  const std::string symbol,
  const std::vector<std::string> species_elements,
  const std::vector<int> stoichiometric_coeff,
  const std::string phase, 
  const std::vector<double>& fit_coeff_limits,
  const std::vector<std::vector<double_type>>& fit_coeff)
{
  Condensate<double_type> species;

  species.name = name;
  species.symbol = symbol;

  species.fit_coeff_limits = fit_coeff_limits;
  species.fit_coeff = fit_coeff;

  species.stoichiometric_vector.assign(nb_elements, 0);


  bool is_stoichiometry_complete = true;

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

  }


  if (is_stoichiometry_complete)
  {
    PhaseState phase_state;

    if (phase == "s" || phase == "solid")
      phase_state = PhaseState::solid;
    else if (phase == "l" || phase == "liquid")
      phase_state = PhaseState::liquid;
    else if (phase == "sl" || phase == "solid_liquid")
      phase_state = PhaseState::solid_liquid;
    else
    {
      std::cout << "Phase state " << phase << " of species " << symbol << " not recognised! Setting to solid.\n";
      phase_state = PhaseState::solid;
    }

    species.phase = phase_state;

    condensates.push_back(species);

    //add the current molecule index to their respective elements
    for (size_t i=0; i<species.element_indices.size(); ++i)
      elements[species.element_indices[i]].condensate_list.push_back(condensates.size()-1);
  }
  else 
    std::cout << "Stoichometry of species " << symbol << " incomplete. Neglected!\n";
}



template class CondensedPhase<double>;
template class CondensedPhase<long double>;
}
