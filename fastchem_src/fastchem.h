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


#ifndef _fastchem_h
#define _fastchem_h

#include <vector>
#include <iostream>
#include <string>

#include "options.h"
#include "fastchem_constants.h"
#include "species_struct.h"
#include "input_output_struct.h"
#include "elements/elements.h"
#include "gas_phase/gas_phase.h"
#include "condensed_phase/condensed_phase.h"


namespace fastchem {


//FastChem class
template <class double_type>
class FastChem {
  public:
    FastChem(
      const std::string& model_parameter_file, 
      const unsigned int verbose_level_init);
    FastChem(
      const std::string& element_abundances_file,
      const std::string& species_data_file,
      const std::string& cond_species_data_file,
      const unsigned int verbose_level_init);
    FastChem(
      const std::string& element_abundances_file,
      const std::string& species_data_file,
      const unsigned int verbose_level_init);
    FastChem(const FastChem &obj);

    //function calls to calculate number densities
    unsigned int calcDensities(FastChemInput& input, FastChemOutput& output);

    //public query functions
    std::string getGasSpeciesName(const unsigned int species_index);
    std::string getGasSpeciesSymbol(const unsigned int species_index);
    unsigned int getGasSpeciesIndex(const std::string symbol);

    std::string getElementName(const unsigned int species_index);
    std::string getElementSymbol(const unsigned int species_index);
    unsigned int getElementIndex(const std::string symbol);

    std::string getCondSpeciesName(const unsigned int species_index);
    unsigned int getCondSpeciesIndex(const std::string symbol);
    std::string getCondSpeciesSymbol(const unsigned int species_index);

    unsigned int getGasSpeciesNumber() {return gas_phase.nb_species;}
    unsigned int getElementNumber() {return element_data.nb_elements;}
    unsigned int getMoleculeNumber() {return gas_phase.nb_molecules;}
    unsigned int getCondSpeciesNumber() {return condensed_phase.nb_condensates;}

    double getElementAbundance(const unsigned int species_index);
    std::vector<double> getElementAbundances();

    double getGasSpeciesWeight(const unsigned int species_index);
    double getElementWeight(const unsigned int species_index);
    double getCondSpeciesWeight(const unsigned int species_index);

    //functions to set internal variables during runtime
    //they will override any read-in values
    void setElementAbundances(std::vector<double> abundances);

    void setVerboseLevel(const unsigned int level) { 
      if (level > 4) options.verbose_level = 4; else options.verbose_level = level;}

    void setParameter(const std::string& parameter, const double_type param_value);
    void setParameter(const std::string& parameter, const bool param_value);
    void setParameter(const std::string& parameter, const unsigned int param_value);

  private:
    FastChemOptions<double_type> options;
    
    ElementData<double_type> element_data;
    GasPhase<double_type> gas_phase;
    CondensedPhase<double_type> condensed_phase;

    bool is_initialised = false;
    bool is_busy = false;

    //Initialisation functions
    void init();

    //Functions for the calculations of the number densities
    unsigned int calcDensity(
      const double temperature,
      const double pressure,
      const bool use_previous_solution,
      std::vector<double>& number_densities,
      double& total_element_density, 
      double& mean_molecular_weight,
      std::vector<unsigned int>& element_conserved,
      unsigned int& nb_chemistry_iterations);

    unsigned int equilibriumCondensation(
      const double temperature,
      const double pressure,
      std::vector<double>& number_densities,
      std::vector<double>& number_densities_cond,
      std::vector<double>& element_cond_degree,
      double& total_element_density, 
      double& mean_molecular_weight,
      std::vector<unsigned int>& element_conserved,
      unsigned int& nb_chemistry_iterations,
      unsigned int& nb_cond_iterations,
      unsigned int& nb_combined_iter);

    void rainoutCondensation(FastChemInput& input, FastChemOutput& output);
};


}

#endif
