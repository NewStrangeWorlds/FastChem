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


#include "options.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>



namespace fastchem {


//Read the FastChem parameter file
template <class double_type>
bool FastChemOptions<double_type>::readParameterFile(const std::string& model_parameter_file)
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
  int nb_value;

  //element abundance file
  std::getline(file, line);
  std::getline(file, line);
  std::istringstream input(line);

  input >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    element_abundances_file = file_name;

  std::getline(file, line);

  //gas phase species data file
  std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    species_data_file = file_name;
  
  file_name = "";
  
   //condensate species data file
  input >> file_name;

  if (file_name == "")
    file_name = "none";

  condensates_data_file = file_name;

  std::getline(file, line);

  //chemistry accuracy
  std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    chem_accuracy = parameter_value;

  std::getline(file, line);
  
  newton_err = chem_accuracy;
  cond_accuracy = chem_accuracy;

  //element conservation accuracy
  std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    element_conserve_accuracy = parameter_value;

  std::getline(file, line);

  //chemistry iteration number
  std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> nb_value;

  if (!(nb_value > 0))
    return initialization_status;
  else
    nb_max_fastchem_iter = nb_value;

  nb_max_cond_iter = nb_value;
  nb_chem_cond_iter = nb_value;

  std::getline(file, line);


  //solver iteration number
  std::getline(file, line); std::getline(file, line);
  
  input.str(line); input.clear();
  input >> nb_value;

  if (!(nb_value > 0))
    return initialization_status;
  else
    nb_max_newton_iter = nb_value;

  nb_max_neldermead_iter = nb_value;
  nb_max_bisection_iter = nb_value;

  std::getline(file, line);

  //optional element data file
  std::getline(file, line); std::getline(file, line);
  
  file_name = "";

  input.str(line); input.clear();
  input >> file_name;

  if (file_name == "")
    chemical_element_file = "";
  else
    chemical_element_file = file_name;


  initialization_status = true;


  if (verbose_level >= 3)
  {
    std::cout << "\nParameter values read from file " << model_parameter_file <<  "\n";
    std::cout << "element abundances file: " << element_abundances_file << "\n";
    std::cout << "elements data file (optional): " << chemical_element_file << "\n";
    std::cout << "species data file: " << species_data_file << "\n";
    std::cout << "chemistry accuracy: " << chem_accuracy << "\n";
    std::cout << "Newton's method error: " << newton_err << "\n";
    std::cout << "max number of chemistry iterations: " << nb_max_fastchem_iter << "\n";
    std::cout << "max number of condensation iterations: " << nb_max_cond_iter << "\n";
    std::cout << "max number of coupled gas phase-condensation iterations: " << nb_chem_cond_iter << "\n";
    std::cout << "max number of Newton iterations: " << nb_max_newton_iter << "\n";
    std::cout << "max number of Nelder-Mead iterations: " << nb_max_neldermead_iter << "\n";
    std::cout << "max number of bisection iterations: " << nb_max_bisection_iter << "\n";
    std::cout << "\n";
  }

  return initialization_status;
}


template struct FastChemOptions<double>;
template struct FastChemOptions<long double>;
}
