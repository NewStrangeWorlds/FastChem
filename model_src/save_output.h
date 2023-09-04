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


#ifndef _save_output_h
#define _save_output_h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include "../fastchem_src/input_output_struct.h"
#include "../fastchem_src/fastchem.h"



bool saveChemistryOutput(
  std::string &file_path,
  bool output_mixing_ratios,
  fastchem::FastChemInput &input,
  fastchem::FastChemOutput &output,
  std::vector<std::string> species_symbols,
  std::vector<double> gas_number_density)
{
  size_t nb_species = species_symbols.size();
  

  std::fstream file(file_path.c_str(), std::ios::out);
  

  if (!file.fail())
  {
    file << std::setw(16) << std::left << "#p (bar)" << "\t"
         << std::setw(16) << std::left << "T (K)" << "\t"
         << std::setw(16) << std::left << "n_<tot> (cm-3)" << "\t"
         << std::setw(16) << std::left << "n_g (cm-3)" << "\t"
         << std::setw(16) << std::left << "m (u)";

    for (size_t i=0; i<nb_species; i++)
      file << "\t" << std::setw(16) << std::left << species_symbols[i];

    file << "\n";

    for (size_t i=0; i<input.pressure.size(); i++)
    {
      file << std::setprecision(10) << std::scientific
           << input.pressure[i] << "\t"
           << input.temperature[i] << "\t"
           << output.total_element_density[i] << "\t"
           << gas_number_density[i] << "\t"
           << output.mean_molecular_weight[i];

      for (size_t j=0; j<nb_species; j++)
        if (!output_mixing_ratios) 
          file << "\t" << output.number_densities[i][j];
        else 
          file << "\t" << output.number_densities[i][j] /gas_number_density[i];

      file << "\n";
    }

    file.close();
  }
  else
  {
    std::cout << "Unable to open chemistry outout file " << file_path << "\n";
    return false;
  }
    
  return true;
}



bool saveCondOutput(
  const std::string &file_path,
  const fastchem::FastChemInput &input,
  const fastchem::FastChemOutput &output,
  const std::vector<std::string> element_symbols,
  const std::vector<std::string> species_symbols)
{
  size_t nb_species = species_symbols.size();
  size_t nb_elements = element_symbols.size();
  

  std::fstream file(file_path.c_str(), std::ios::out);
  

  if (!file.fail())
  {
    file << std::setw(16) << std::left << "#p (bar)" << "\t"
         << std::setw(16) << std::left << "T (K)";

    for (size_t i=0; i<nb_elements; i++)
      file << "\t" << std::setw(16) << std::left << element_symbols[i];
    
    for (size_t i=0; i<nb_species; i++)
      file << "\t" << std::setw(16) << std::left << species_symbols[i];

    file << "\n";

    for (size_t i=0; i<input.pressure.size(); i++)
    {
      file << std::setprecision(10) << std::scientific
           << input.pressure[i] << "\t"
           << input.temperature[i];

      for (size_t j=0; j<nb_elements; ++j)
        file << "\t" << output.element_cond_degree[i][j];

      for (size_t j=0; j<nb_species; j++)
        file << "\t" << output.number_densities_cond[i][j];

      file << "\n";
    }

    file.close();
  }
  else
  {
    std::cout << "Unable to open condensation outout file " << file_path << "\n";
    return false;
  }
    
  return true;
}




bool saveMonitorOutput(
  std::string &file_path,
  fastchem::FastChemInput &input,
  fastchem::FastChemOutput &output,
  std::vector<std::string> element_symbols,
  std::vector<double> gas_number_density)
{
  size_t nb_elements = element_symbols.size();


  std::fstream file(file_path.c_str(), std::ios::out);


  if (!file.fail())
  {
    file << std::setw(16) << std::left << "#grid point" << "\t"
         << std::setw(16) << std::left << "iterations" << "\t"
         << std::setw(16) << std::left << "chem_iter" << "\t"
         << std::setw(16) << std::left << "cond_iter" << "\t"
         << std::setw(16) << std::left << "converged" << "\t"
         << std::setw(16) << std::left << "elem_conserved" << "\t"
         << std::setw(16) << std::left << "p (bar)" << "\t"
         << std::setw(16) << std::left << "T (K)" << "\t"
         << std::setw(16) << std::left << "n_<tot> (cm-3)" << "\t"
         << std::setw(16) << std::left << "n_g (cm-3)" << "\t"
         << std::setw(16) << std::left << "m (u)";
    
    for (unsigned int i=0; i<nb_elements; i++)
      file << "\t" << std::setw(5) << std::left << element_symbols[i];

    file << "\n";

    std::vector<std::string> output_flags {"fail", "ok"};
    std::vector<std::string> convergence_flags {"yes", "no conv", "init fail"};


    std::vector<unsigned int> all_elements_conserved(input.pressure.size(), 1);

    for (size_t i=0; i<input.pressure.size(); ++i)
      if (std::any_of(output.element_conserved[i].begin(), output.element_conserved[i].end(), [](unsigned int i){return i==0;}))
        all_elements_conserved[i] = 0;

    for (unsigned int i=0; i<input.pressure.size(); i++)
    {
      std::string c_conv;

      if (output.fastchem_flag[i] == fastchem::FASTCHEM_SUCCESS)
        c_conv = output_flags[1];
      else
        c_conv = output_flags[0];


      file << std::setw(16) << std::left << i << "\t"
           << std::setw(16) << std::left << output.nb_iterations[i] << "\t"
           << std::setw(16) << std::left << output.nb_chemistry_iterations[i] << "\t"
           << std::setw(16) << std::left << output.nb_cond_iterations[i] << "\t"
           << std::setw(16) << std::left << c_conv << "\t"
           << std::setw(16) << std::left << output_flags[all_elements_conserved[i]] << "\t";

      file << std::setprecision(10) << std::scientific
           << input.pressure[i] << "\t" << input.temperature[i] << "\t"
                                << output.total_element_density[i] << "\t"
                                << gas_number_density[i] << "\t"
                                << output.mean_molecular_weight[i];

      for (unsigned int j=0; j<nb_elements; j++)
        file << "\t" << std::setw(5) << output_flags[output.element_conserved[i][j]];

      file << "\n";
    }

    file.close();
  }
  else
  {
    std::cout << "Unable to open monitor outout file " << file_path << "\n";
    return false;
  }


  return true;
}


#endif
