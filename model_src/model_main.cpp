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

#include "../fastchem_src/fastchem.h"
#include "../fastchem_src/input_output_struct.h"
#include "../fastchem_src/fastchem_constants.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <functional>

#include <cmath>
#include "read_config.h"
#include "save_output.h"



int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "config file command line parameter missing!\n";

    return 1;
  }

  //first we read in the config file
  std::string config_file_name = argv[1];

  Config config;

  if (!readConfigFile(config_file_name, config))
    return 1;



  //create a FastChem object and set the read-in values from the config file
  fastchem::FastChem<long double> fastchem(
    config.element_abundance_file, 
    config.species_data_file, 
    config.cond_species_data_file,
    config.verbose_level);

  //we could also create a FastChem object using a separate parameter file
  //fastchem::FastChem<long double> fastchem("input/parameters.dat", config.verbose_level); 


  //set the FastChem parameters from the config file
  fastchem.setParameter(std::string("accuracyChem"), (long double) config.chemistry_accuracy);
  fastchem.setParameter(std::string("accuracyElementConservation"), (long double) config.element_conservation_accuracy);
  fastchem.setParameter(std::string("accuracyNewton"), (long double) config.newton_error);
  fastchem.setParameter(std::string("nbIterationsChem"), config.nb_chemistry_iterations);
  fastchem.setParameter(std::string("nbIterationsNewton"), config.nb_newton_iterations);
  fastchem.setParameter(std::string("nbIterationsNelderMead"), config.nb_nelder_mead_iterations);
  fastchem.setParameter(std::string("nbIterationsBisection"), config.nb_bisection_iterations);


  //read in the temperature-pressure structure
  std::vector<double> temperature;
  std::vector<double> pressure;

  std::string line;
  std::fstream file(config.atmosphere_file.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Couldn't open input file: " << config.atmosphere_file << "\n";
    return 0;
  }

  while (std::getline(file, line))
  {
    std::stringstream  line_stream(line);

    double temperature_in;
    double pressure_in;

    if (!(line_stream >> pressure_in >> temperature_in)) continue;

    pressure.push_back(pressure_in);  
    temperature.push_back(temperature_in);
  }

  file.close();

  std::cout << "temperature-pressure profile read in from file:\n";
  for (size_t i=0; i<pressure.size(); i++)
    std::cout << "  " << i << "\t" << pressure[i] << "\t" << temperature[i] << "\n";
  std::cout << "\n";

  //set up the input & output structures and run the chemistry
  fastchem::FastChemInput input;
  fastchem::FastChemOutput output;

  input.temperature = temperature;
  input.pressure = pressure;
  input.equilibrium_condensation = config.calc_condensation;
  input.rainout_condensation = config.rainout_condensation;

  unsigned int fastchem_flag = fastchem.calcDensities(input, output);

  
  if (fastchem_flag == fastchem::FASTCHEM_INITIALIZATION_FAILED)
  {
    std::cout << "FastChem initialisation failed!\n"; return 0;
  }


  //check if all elements have been conserved
  bool elements_conserved = true;
  
  for (auto & i : output.element_conserved)
  {
    if (std::any_of(i.begin(), i.end(), [](unsigned int v){return v==0;}))
      elements_conserved = false;
    
    if (elements_conserved == false) break;
  }

  std::string conserved = "yes";
  if (elements_conserved == false) conserved = "no";


  //print out the status message for the returned FastChem flag
  std::cout << "\nFastChem reports:\n";
  std::cout << "  - " << fastchem::FASTCHEM_MSG[fastchem_flag] << "\n";
  std::cout << "  - elements conserved: " << conserved << "\n\n";


  unsigned int nb_grid_points = pressure.size();
  unsigned int nb_species = fastchem.getGasSpeciesNumber();
  unsigned int nb_elements = fastchem.getElementNumber();
  unsigned int nb_condensates = fastchem.getCondSpeciesNumber();

  //calculate the gas number density via the ideal gas law for the output
  std::vector<double> gas_number_density(pressure);

  for (unsigned int i=0; i<nb_grid_points; i++)
    gas_number_density[i] /= fastchem::CONST_K * temperature[i] * 1e-6;


  std::vector<std::string> species_symbols(nb_species);

  for (size_t i=0; i<nb_species; ++i)
    species_symbols[i] = fastchem.getGasSpeciesSymbol(i);

  std::vector<std::string> element_symbols(species_symbols.begin(), species_symbols.begin()+nb_elements);


  saveMonitorOutput(config.monitor_output_file, input, output, element_symbols, gas_number_density);
  saveChemistryOutput(config.chem_output_file, config.output_mixing_ratios, input, output, species_symbols, gas_number_density);

  if (config.calc_condensation)
  {
    std::vector<std::string> cond_species_symbols(nb_condensates);

    for (size_t i=0; i<nb_condensates; ++i)
      cond_species_symbols[i] = fastchem.getCondSpeciesSymbol(i);

    saveCondOutput(config.cond_output_file, input, output, element_symbols, cond_species_symbols);
  }


  std::cout << "FastChem finished!" << std::endl;

  return 0;
}
