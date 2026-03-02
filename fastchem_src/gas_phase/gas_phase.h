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


class GasPhase {
  public:
    GasPhase(
      FastChemOptions& options_,
      ElementData& element_data_);
    GasPhase(
      const GasPhase &obj,
      FastChemOptions& options_,
      ElementData& element_data_);

    std::vector<Molecule> molecules;
    std::vector< ChemicalSpecies* > species;

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
    FastChemOptions& options;
    ElementData& element_data;
    std::vector<Element>& elements;

    GasPhaseSolver solver;
    std::vector<unsigned int> element_calculation_order;

    unsigned int e_ = FASTCHEM_UNKNOWN_SPECIES; //electron element index

    void init();

    bool readSpeciesData(const std::string& file_path);
    void addMolecule(
      const std::string& name,
      const std::string& symbol,
      const std::vector<std::string>& species_elements,
      const std::vector<int>& stoichiometric_coeff,
      const std::vector<double>& mass_action_coeff,
      const int charge);
    bool checkForDuplicates();

    unsigned int determineSolverOrder(const Element& species);
    void determineSolverOrder();
    void determineElementCalculationOrder();
    
    void setMoleculeAbundances();
    void createMoleculeLists();

    void calculateElementDensities(
      Element& species,
      const double gas_density,
      bool use_backup_solver,
      double& n_major);
    double calculateMoleculeDensities(
      Element& species, const double gas_density);

    void calculateElectronDensities(
      Element& species,
      const double& old_number_density,
      const double gas_density);
    void calculateSinglyIonElectrons(
      Element& electron,
      const double& old_number_density);
    void calculateMultIonElectrons(
      Element& electron,
      const double& old_number_density,
      const double& gas_density);
};



}

#endif
