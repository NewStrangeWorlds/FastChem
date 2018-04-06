/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
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


//Read the FastChem parameter file
template <class double_type>
bool FastChem<double_type>::readParameterFile(const std::string& model_parameter_file)
{
  bool initialization_status = false;

  std::fstream file;

  file.open(model_parameter_file.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open FastChem parameter file " << model_parameter_file << "\n";

    return initialization_status;
  }


  std::string file_name, line;
  double parameter_value;


  std::getline(file, line);
  file >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    element_abundances_file = file_name;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    chemical_element_file = file_name;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    species_data_file = file_name;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    accuracy = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    accuracy_delta = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    newton_err = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    nb_max_fastchem_iter = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    nb_max_pressure_iter = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    nb_max_newton_iter = parameter_value;


  initialization_status = true;


  if (verbose_level >= 3)
  {
    std::cout << "Parameter values read: \n";
    std::cout << "element_abundances " << element_abundances_file << "\n";
    std::cout << "element_file " << chemical_element_file << "\n";
    std::cout << "species_data " << species_data_file << "\n";
    std::cout << "accuracy " << accuracy << "\n";
    std::cout << "accuracy_delta " << accuracy_delta << "\n";
    std::cout << "newton_err " << newton_err << "\n";
    std::cout << "iter_max " << nb_max_fastchem_iter << "\n";
    std::cout << "k_max " << nb_max_pressure_iter << "\n";
    std::cout << "mu_max " << nb_max_newton_iter << "\n";
  }


  return initialization_status;
}



//Read the chemical elements file
template <class double_type>
bool FastChem<double_type>::readElementList()
{
  std::fstream file(chemical_element_file.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open chemical element file " << chemical_element_file << "\n";

    return false;
  }


  std::string line;

  while (std::getline(file, line))
  {
    ChemicalElement<double_type> element;
    std::istringstream input(line);

    input >> element.symbol >> element.name >> element.atomic_weight;

    chemical_elements.push_back(element);

    chemical_elements.back().abundance = 0;
  }


  file.close();


  nb_chemical_elements = chemical_elements.size();


  return true;
}



template <class double_type>
bool FastChem<double_type>::readSpeciesData()
{
  std::fstream file(species_data_file.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "Unable to open species data file " << species_data_file << "\n";

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

    input >> symbol >> name;


    if (name == ":")
      name = "";
    else
      input >> element_string;


    std::vector<std::string> elements;
    std::vector<int> stoichometric_coeff;
    int charge = 0;


    while ( (input >> element_string >> stoichometric_coeff_string) && (element_string != "#") )
    {
      elements.push_back(element_string);
      stoichometric_coeff.push_back(std::stoi( stoichometric_coeff_string ));

      if (elements.back() == "e-") charge = -stoichometric_coeff.back();
    }


    std::string line;
    std::getline(file, line);

    std::istringstream ma_input(line);

    std::vector<double_type> mass_action_coeff;


    double ma_coefficient;

    while (ma_input >> ma_coefficient)
      mass_action_coeff.push_back(ma_coefficient);

    addMolecule(name, symbol, elements, stoichometric_coeff, mass_action_coeff, charge);

    //blank separation line
    std::getline(file, line);
  }




  file.close();



  nb_molecules = molecules.size();

  nb_species = nb_molecules + nb_elements;


  return true;
}



//Read the elemental abundances file
template <class double_type>
bool FastChem<double_type>::readElementAbundances()
{

  std::fstream file(element_abundances_file.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "Unable to open element abundances file " << element_abundances_file << "\n";

    return false;
  }



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


    setElementAbundance(symbol, abundance);

    addAtom(symbol);
  }


  file.close();


  nb_elements = elements.size();


  return true;
}


template class FastChem<double>;
template class FastChem<long double>;


}
