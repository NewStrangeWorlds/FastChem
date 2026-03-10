
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


#ifndef _species_struct_h
#define _species_struct_h


#include <vector>
#include <string>

#include "fastchem_constants.h"


namespace fastchem {


enum PhaseState {gas, liquid, solid, solid_liquid};
const std::vector<std::string> phase_state_string{"gas", "liquid", "solid", "solid/liquid"};


//Class for the basic chemical element data
//Those are read in from the elemental abundance file
//and are later cross-correlated with the elements in the network
//template <class double>
struct ChemicalElementData
{
  std::string symbol;
  std::string name;
  
  double atomic_weight;
  double abundance;   //element abundance
};



//Parent class of all species
struct ChemicalSpecies
{
  std::string symbol;
  std::string name;

  double weight = 0.0;
  int charge = 0;
  PhaseState phase;

  double abundance = 0.0;
  double number_density = 0.0;
  double log_number_density = LOG_DENSITY_FLOOR;
};


//forward declarations
struct Molecule;
struct Condensate;



//Class for the elements
struct Element : public ChemicalSpecies
{
  unsigned int element_data_index;
  unsigned int index;

  unsigned int solver_order = 0;

  std::vector<unsigned int> molecule_list;         //contains the list of molecule indices the element is part of
  std::vector<unsigned int> major_molecules_inc;   //contains the list of 'major' molecules that includes the element
  std::vector<unsigned int> major_molecules_exc;   //contains the list of 'major' molecules that doesn't include the element
  std::vector<unsigned int> minor_molecules;       //contains the list of 'minor' molecules
  
  std::vector<unsigned int> condensate_list;       //contains the list of condensate indices the element is part of

  unsigned int element_conserved = 0;              //check if element is conserved during calculation, for electrons this is charge conservation

  double epsilon = 0.0;
  double number_density_maj = 0.0;
  double number_density_min = 0.0;

  double phi = 0.0;                           //relative abundance in the gas phase (minus condensates)

  double degree_of_condensation = 0.0;
  bool fixed_by_condensation = false;

  void calcMinorSpeciesDensities(
    const std::vector< Molecule > &molecules);
  void calcEpsilon(
    const std::vector< Element > &elements);
  void checkN(
    const double& min_limit, const double log_gas_density);
  bool checkElementConservation(
    const std::vector<Molecule>& molecules,
    const std::vector< Condensate >& condensates,
    const double total_density,
    const double& accuracy);
  bool checkChargeConservation(
    const std::vector<Molecule>& molecules,
    const double& accuracy);
  void calcDegreeOfCondensation(
    const std::vector< Condensate > &condensates,
    const double total_element_density);
  void normalisePhi(const double element_phi_sum);
};



//Class for the molecules
struct Molecule : public ChemicalSpecies
{
  std::vector<unsigned int> element_indices;
  std::vector<int> stoichiometric_vector;

  std::vector<double> mass_action_coeff;
  double mass_action_constant = 0.0;

  double abundance_scaled = 0.0;

  double sigma = 0.0;
  double sum = 0.0;

  void calcMassActionConstant(const double temperature);
  void calcLogNumberDensity(const std::vector< Element >& elements);
  void calcNumberDensity(const std::vector< Element >& elements);
  void checkN(const double& min_limit, const double log_gas_density);
};



//Parent class for condensates
struct Condensate : public ChemicalSpecies
{
  std::vector<unsigned int> element_indices;
  std::vector<int> stoichiometric_vector;

  std::vector<std::vector<double>> fit_coeff;
  std::vector<double> fit_coeff_limits;

  double mass_action_constant = 0.0;

  double log_activity = 0;
  double activity_correction = 0;
  double tau = 0;
  double log_tau = 0;
  double max_number_density = 0;

  bool linear_system_remove = false;
  bool is_calculated = false;

  unsigned int reference_element = FASTCHEM_UNKNOWN_SPECIES;     //the element, the degree of condensation is defined for
  double degree_of_condensation = 0;

  void calcMassActionConstant(const double temperature);
  void calcActivity(
    const double temperature, 
    const std::vector<Element>& elements,
    const bool use_data_validity_limits);
  double calcActivity(
    const double temperature, 
    const std::vector<Element>& elements,
    const std::vector<double> elem_number_densities,
    const bool use_data_validity_limits);
  void findReferenceElement(
    const std::vector<Element>& elements);
  void degreeOfCondensation(
    const std::vector<Element>& elements, 
    const double total_element_density);
  void maxDensity(
    const std::vector< Element >& elements, double total_number_density);
};



}
#endif