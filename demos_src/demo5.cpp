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


#include "../fastchem_src/fastchem.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <functional>

#include <cmath>


//forward declaration
bool read_config_file(std::string file_path, std::string& fastchem_options_file, std::string& atmosphere_file,
                      std::string& chem_output_file, std::string& monitor_output_file,
                      unsigned int& verbose_level, bool& output_mixing_ratios);




int main(int argc, char *argv[])
{

  //check if command line parameter is present
  if (argc == 1)
  {
    std::cout << "config file command line parameter missing!\n";

    return 1;
  }

  
  //set up the config parameters
  std::string config_file_name = argv[1];

  std::string fastchem_options_file;
  std::string atmosphere_file;
  std::string chem_output_file;
  std::string monitor_output_file;
  unsigned int verbose_level = 4;
  bool output_mixing_ratios = false;


  //read the config file
  if (!read_config_file(config_file_name, fastchem_options_file,atmosphere_file,chem_output_file, monitor_output_file, verbose_level, output_mixing_ratios))
    return 1;


  //read the temperature&pressure input structure
  std::vector<double> temperature;
  std::vector<double> pressure;


  std::fstream file(atmosphere_file.c_str(), std::ios::in);

  //check if file is present
  if (file.fail())
  {
    std::cout << "Unable to temperature-pressure file: " << atmosphere_file << "\n";
    return false;
  }


  //read the structure
  //this assumes that the columns in the file are: temperature pressure
  //temperature in K, pressure in bar
  std::string line;
  
  //read in the comment header
  std::getline(file, line);

  while(std::getline(file, line))
  {
    std::stringstream  line_stream(line);

    double temperature_in;
    double pressure_in;


    //if you want to read a file with differently ordered columns, you have to adapt the following lines
    if (!(line_stream >> temperature_in >> pressure_in)) break;
    //if (!(line_stream >> pressure_in >> temperature_in)) break;


    pressure.push_back(pressure_in*1.e+6);  //FastChem uses pressure in dyn cm-2
    temperature.push_back(temperature_in);
  }

  file.close();



  //set up the data structures for the output
  unsigned int nb_grid_points = pressure.size();

  std::vector<double> total_density(pressure);

  for (unsigned int i=0; i<nb_grid_points; i++)
    total_density[i] /= fastchem::CONST_K * temperature[i];


  for (unsigned int i=0; i<nb_grid_points; i++)
    std::cout << i << "\t" << pressure[i] << "\t" << temperature[i] << std::endl;


  std::vector< std::vector<double> > densities;
  densities.resize(nb_grid_points);

  std::vector<double> h_densities(nb_grid_points, 0.0);
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);


  std::vector< std::vector<unsigned int> > element_conserved;
  element_conserved.resize(nb_grid_points);

  std::vector<unsigned int> nb_pressure_iterations (nb_grid_points, 0);
  std::vector<unsigned int> nb_chemistry_iterations (nb_grid_points, 0);
  std::vector<unsigned int> fastchem_flags (nb_grid_points, 0);




  //now, we create a FastChem instance - in this case the double version
  //for the constructor, we need to supply the location of the FastChem parameter file and a verbose level
  fastchem::FastChem<double> fastchem(fastchem_options_file, verbose_level);

  //set terminal output of FastChem
  fastchem.setVerboseLevel(verbose_level);

  
  //here, we call FastChem point-wise, looping through the T-p structure manually
  //we use the version without the full diagnostic output
  for (size_t i=0; i<nb_grid_points; ++i)
    fastchem_flags[i] = fastchem.calcDensities(temperature[i],pressure[i],densities[i],h_densities[i], mean_molecular_weights[i]);
  

  //now we write the output
  //here, we only have the chemistry output, since no diagnostics were reported back by FastChem
  file.open(chem_output_file.c_str(), std::ios::out);

  unsigned int nb_species = fastchem.getSpeciesNumber();  //ask FastChem how many species it contains


  file << std::setw(16) << std::left << "#P (bar)" << "\t"
       << std::setw(16) << std::left << "T(k)" << "\t"
       << std::setw(16) << std::left << "n_<H> (cm-3)" << "\t"
       << std::setw(16) << std::left << "n_g (cm-3)" << "\t"
       << std::setw(16) << std::left << "m(u)";
  for (unsigned int i=0; i<nb_species; i++)
    file << "\t" << std::setw(16) << std::left << fastchem.getSpeciesSymbol(i); //query FastChem for the chemical symbols

  file << "\n";


  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    file << std::setprecision(10) << std::scientific
         << pressure[i]/1.e6 << "\t"
         << temperature[i] << "\t"
         << h_densities[i] << "\t"
         << total_density[i] << "\t"
         << mean_molecular_weights[i];

    for (unsigned int j=0; j<nb_species; j++)
      if (!output_mixing_ratios) file << "\t" << densities[i][j];
      else file << "\t" << densities[i][j] /total_density[i];

    file << "\n";
  }

  file.close();


  std::cout << "FastChem finished! " << std::endl;

  return 0;
}
