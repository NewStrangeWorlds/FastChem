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

#include "elements.h"

#include "../fastchem_constants.h"
#include "chemical_element_data.h"
#include "../species_struct.h"


namespace fastchem {


//Read the chemical elements file
template <class double_type>
bool ElementData<double_type>::readElementList(const std::string& file_path)
{ 
  //if no file for the element data has been found in the options file
  //use the standard set
  if (file_path == "")
  { 
    chemical_element_data = standard_chemical_element_data;
    nb_chemical_element_data = chemical_element_data.size();
    return true;
  }


  std::fstream file(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open chemical element file " << file_path << "\n";

    return false;
  }


  chemical_element_data.reserve(200);

  std::string line;
  
  //ignore header
  std::getline(file, line);

  while (std::getline(file, line))
  {
    ChemicalElementData element;
    std::istringstream input(line);

    input >> element.symbol >> element.name >> element.atomic_weight;

    chemical_element_data.push_back(element);

    chemical_element_data.back().abundance = 0;
  }

  file.close();
  

  chemical_element_data.shrink_to_fit();

  nb_chemical_element_data = chemical_element_data.size();

  return true;
}


//Read the elemental abundances file
template <class double_type>
bool ElementData<double_type>::readElementAbundances(const std::string& file_path)
{
  std::fstream file(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open element abundances file " << file_path << "\n";

    return false;
  }


  elements.reserve(200);

  std::string line;

  //header
  std::getline(file, line);

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string symbol;
    double abundance;

    input >> symbol >> abundance;

    abundance = std::pow(10., abundance - 12.);

    setAbundance(symbol, abundance);

    add(symbol);
  }

  file.close();


  elements.shrink_to_fit();
  nb_elements = elements.size();

  return true;
}



template class ElementData<double>;
template class ElementData<long double>;

} 
