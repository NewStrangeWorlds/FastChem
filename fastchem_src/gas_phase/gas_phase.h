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

#ifndef _gas_phase_h
#define _gas_phase_h

#include <string>
#include <vector>

#include "../species_struct.h"
#include "../elements/elements.h"
#include "../options.h"

#include "solver.h"


namespace fastchem {


template <class double_type>
class GasPhase {
  public:
    GasPhase(
      FastChemOptions<double_type>& options_,
      ElementData<double_type>& element_data_);
    GasPhase(
      const GasPhase &obj,
      FastChemOptions<double_type>& options_,
      ElementData<double_type>& element_data_);

    std::vector<Molecule<double_type>> molecules;
    std::vector< ChemicalSpecies<double_type>* > species;

    size_t nb_molecules = 0;
    size_t nb_elements = 0;
    size_t nb_species = 0;

    bool is_initialised = false;

    bool calculate(
      const double temperature_gas, const double gas_density, unsigned int& nb_iterations);

    double meanMolecularWeight(const double gas_density);
    double totalElementDensity();

    void reInitialise();
  private:
    FastChemOptions<double_type>& options;
    ElementData<double_type>& element_data;
    std::vector<Element<double_type>>& elements;

    GasPhaseSolver<double_type> solver;
    std::vector<unsigned int> element_calculation_order;

    unsigned int e_ = FASTCHEM_UNKNOWN_SPECIES; //electron element index

    void init();

    bool readSpeciesData(const std::string& file_path);
    void addMolecule(
      const std::string& name,
      const std::string& symbol,
      const std::vector<std::string>& species_elements,
      const std::vector<int>& stoichiometric_coeff,
      const std::vector<double_type>& mass_action_coeff,
      const int charge);
    bool checkForDuplicates();

    unsigned int determineSolverOrder(const Element<double_type>& species);
    void determineSolverOrder();
    void determineElementCalculationOrder();
    
    void setMoleculeAbundances();
    void createMoleculeLists();

    void calculateElementDensities(
      Element<double_type>& species,
      const double_type gas_density,
      bool use_backup_solver,
      double_type& n_major);
    double_type calculateMoleculeDensities(
      Element<double_type>& species, const double_type gas_density);

    void calculateElectronDensities(
      Element<double_type>& species,
      const double_type& old_number_density,
      const double_type gas_density);
    void calculateSinglyIonElectrons(
      Element<double_type>& electron,
      const double_type& old_number_density);
    void calculateMultIonElectrons(
      Element<double_type>& electron,
      const double_type& old_number_density,
      const double_type& gas_density);
};



}

#endif
