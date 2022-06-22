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


#include "fastchem.h"
#include "chemical_element_data.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>


namespace fastchem {



template <class double_type>
void FastChem<double_type>::init()
{
  //read input files
  if (    readElementList() == false
       || readElementAbundances() == false
       || readSpeciesData() == false)
  {
    std::cout << "FastChem data file loading failed!\n";

    is_initialized = false;

    return;
  }


  //setup the element list without the electron
  elements_wo_e.resize(0);
  elements_wo_e.reserve(nb_elements);

  for (auto & i : elements)
    if (i.symbol != "e-") elements_wo_e.push_back(&i);

  elements_wo_e.shrink_to_fit();


  //init accuracy limit based on the numerical precision
  if (std::numeric_limits<double_type>::max_exponent10 > 1000)
  {
    options.element_density_minlimit = 1e-512L;
    options.molecule_density_minlimit = 1e-512L;
  }
  else
  {
    options.element_density_minlimit = 1e-155;
    options.molecule_density_minlimit = 1e-155;
  }

  
  nb_species = nb_elements + nb_molecules;
  
  //setting up the list of all species
  species.reserve(nb_species);

  for (size_t i=0; i<nb_elements; ++i)
    species.push_back(&elements[i]);

  for (size_t i=0; i<nb_molecules; ++i)
    species.push_back(&molecules[i]);

  
  determineElementCalculationOrder();

  setMoleculeAbundances();

  determineSolverOrder();

  for (auto & i : elements) i.calcEpsilon(elements);


  //some diagnostic output in case you'd like to see it
  if (options.verbose_level >= 4)
  {
    std::cout << "\nSpecies list: \n";
    for (size_t i=0; i<molecules.size(); ++i)
    {
      std::cout << "  " << molecules[i].name << "\t" << molecules[i].symbol << "\n";
      
      std::cout << "    lnK coeff: ";
      for (size_t j=0; j<molecules[i].mass_action_coeff.size(); ++j)
        std::cout << molecules[i].mass_action_coeff[j] << "\t";
      std::cout << "\n";
      
      std::cout << "    stoichiometry: ";
      for (size_t j=0; j<molecules[i].stoichiometric_vector.size(); ++j)
        std::cout << molecules[i].stoichiometric_vector[j] << " ";
      std::cout << "\n";

      std::cout << "    elements: ";
      for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
        std::cout << elements[molecules[i].element_indices[j]].symbol << ", index: " << molecules[i].element_indices[j] << "; ";
      std::cout << "\n";

      std::cout << "    charge: " << molecules[i].charge << "\n";
    }


    std::cout << "\nElement list: \n";
    for (size_t i=0; i<elements.size(); ++i)
      std::cout << "  " << elements[i].name << "\t" << elements[i].symbol << "\t" 
                << elements[i].element_data_index << "\t" << elements[i].solver_order << "\t" 
                << elements[i].abundance << "\n";


    std::cout << "\nElement calculation order:\n";
      for (size_t i=0; i<element_calculation_order.size(); ++i)

    std::cout << "  " << elements[element_calculation_order[i]].symbol << "\t" << elements[element_calculation_order[i]].abundance << "\n";

    std::cout << "\n";
  }


  if (options.verbose_level >= 2)
  {
    std::cout << "\nConsidered species in FastChem:\n";

    for (size_t i=0; i<nb_species; ++i)
      std::cout << "  " << species[i]->symbol << "\t" << species[i]->name << "\t" << species[i]->abundance << "\t" << species[i]->weight << "\n";

    std::cout << "\n";
  }


  if (checkForDuplicates() == true)
  {
    is_initialized = false;
    return;
  }


  /*std::fstream file("nu.dat", std::ios::out);
  
  for (auto & i : molecules)
  {
    for (auto & j : i.stoichiometric_vector)
      file << j << "\t";
    
    file << "\n";
  }

  file.close();*/


  e_ = getElementIndex("e-");


  if (options.verbose_level >= 1)
    std::cout << "\nFastChem initialisation summary:\n"
              << "  number of species: " << nb_species
              << "    elements: " << nb_elements
              << "    molecules: " << nb_molecules
              << "    chemical elements: " << nb_chemical_element_data << "\n\n";

  is_initialized = true;
}



//Reinitializes certain internal data after the element abundances were changed
template <class double_type>
void FastChem<double_type>::reInitialiseFastChem()
{
  //reset the calculation order
  element_calculation_order.resize(0);

  //order the elements according to their abundances
  determineElementCalculationOrder();

  //update the abundances of the molecules
  setMoleculeAbundances();

  //update the solver order for the new abundances
  determineSolverOrder();
  
  //recalculate the epsilons
  for (auto & i : elements) i.calcEpsilon(elements);
}


//update the definition of the molecular abundances
template <class double_type>
void FastChem<double_type>::setMoleculeAbundances()
{
  for (auto & i : molecules)
  {
    auto it = std::min_element(
      i.element_indices.begin(),
      i.element_indices.end(), 
      [&] (const size_t a, const size_t b) {
        if (elements[a].abundance != 0) 
          return elements[a].abundance < elements[b].abundance;
        else 
          return false; 
        });

    i.abundance = elements[*it].abundance;
    
    //scaled abundances
    it = std::min_element(
      i.element_indices.begin(),
      i.element_indices.end(), 
      [&] (const size_t a, const size_t b) {
        if (elements[a].abundance != 0) 
          return elements[a].abundance/i.stoichiometric_vector[a] < elements[b].abundance/i.stoichiometric_vector[b];
        else 
          return false; 
        });

    i.abundance_scaled = elements[*it].abundance/i.stoichiometric_vector[*it];
  }

  createMoleculeLists();
}



//looks through the list of all species and tries to fund duplicate entries
template <class double_type>
bool FastChem<double_type>::checkForDuplicates()
{
  //make a copy of the species vector
  auto species_cp = species;
  
  //sort with respect to the species' symbols
  std::sort(
    species_cp.begin(), 
    species_cp.end(),
    [&](ChemicalSpecies<double_type>* a, ChemicalSpecies<double_type>* b) {return a->symbol > b->symbol;} );
  
  //try to find adjacent ones
  auto it = std::adjacent_find(
    species_cp.begin(), 
    species_cp.end(), 
    [&](ChemicalSpecies<double_type>* a, ChemicalSpecies<double_type>* b) {return a->symbol == b->symbol;} );
  
  bool duplicate_exits = (it != species_cp.end());

  if (duplicate_exits)
    std::cout << "Species " << (*it)->symbol << " seems to appear twice in the species data file. Please check!\n";

  return duplicate_exits;
}



template <class double_type>
void FastChem<double_type>::createMoleculeLists()
{
  //reset the lists
  for (auto & i : elements)
  {
    i.major_molecules_inc.resize(0);
    i.major_molecules_inc.reserve(nb_molecules);
    
    i.major_molecules_exc.resize(0);
    i.major_molecules_exc.reserve(nb_molecules);

    i.minor_molecules.resize(0);
    i.minor_molecules.reserve(nb_molecules);
  }


  for (size_t i=0; i<nb_molecules; ++i)
  {

    for (size_t j=0; j<nb_elements; ++j)
    {
       if (elements[j].abundance > molecules[i].abundance)
         elements[j].minor_molecules.push_back(i);
       else
       {

         if (molecules[i].stoichiometric_vector[j] == 0)
           elements[j].major_molecules_exc.push_back(i);
         else
           elements[j].major_molecules_inc.push_back(i);
    
       }
    }
  }


  for (auto & i : elements)
  {
    i.major_molecules_inc.shrink_to_fit();
    i.major_molecules_exc.shrink_to_fit();
    i.minor_molecules.shrink_to_fit();
  }


  if (options.verbose_level >= 4)
  {
    std::cout << "\nMolecule lists for each element: \n";

    for (size_t j=0; j<nb_elements; ++j)
    {
      std::cout << "  element " << elements[j].symbol << "\n";

      std::cout << "    major elements inc:\n";
      for (size_t i=0; i<elements[j].major_molecules_inc.size(); ++i)
        std::cout << "    " << molecules[elements[j].major_molecules_inc[i]].symbol << "\n";

      std::cout << "    major elements exc:\n";
      for (size_t i=0; i<elements[j].major_molecules_exc.size(); ++i)
        std::cout << "    " << molecules[elements[j].major_molecules_exc[i]].symbol << "\n";

      std::cout << "    minor elements:\n";
      for (size_t i=0; i<elements[j].minor_molecules.size(); ++i)
        std::cout << "    " << molecules[elements[j].minor_molecules[i]].symbol << "\n";
    }
  }
}


template class FastChem<double>;
template class FastChem<long double>;
}
