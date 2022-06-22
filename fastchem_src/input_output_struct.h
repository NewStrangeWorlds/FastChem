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


#ifndef _input_output_struct_h
#define _input_output_struct_h


#include <vector>


namespace fastchem {


struct FastChemInput
{
  std::vector<double> temperature; 
  std::vector<double> pressure;
};



struct FastChemOutput
{
  std::vector<std::vector<double>> number_densities;
  std::vector<double> total_element_density;
  std::vector<double> mean_molecular_weight;

  //diagnostic output
  std::vector<std::vector<unsigned int>> element_conserved;
  std::vector<unsigned int> nb_chemistry_iterations;
  std::vector<unsigned int> fastchem_flag;
};


}


#endif
