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


#include "fastchem.h"


#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "input_output_struct.h"


namespace fastchem {


template <class double_type>
unsigned int FastChem<double_type>::calcDensities(
  FastChemInput& input, FastChemOutput& output)
{
  if (!is_initialised)
    return FASTCHEM_INITIALIZATION_FAILED;

  if ((input.rainout_condensation || input.equilibrium_condensation) && condensed_phase.is_initialised == false)
  {
    std::cout << " FastChem is unable to perform calculations including condensation. The data for the condensate species has not been properly initialised!\n";
    return FASTCHEM_INITIALIZATION_FAILED;
  }

  if (is_busy == true)
  {
    std::cout << "FastChem can only be called once per instance!\n";
    return FASTCHEM_IS_BUSY;
  }
  
  //check input values
  if (input.temperature.size() != input.pressure.size())
  {
    std::cout << "Temperature and pressure vector do not have the same size!\n";
    return FASTCHEM_WRONG_INPUT_VALUES;
  }


  is_busy = true;

  size_t nb_gridpoints = input.temperature.size();

  output.element_conserved.assign(nb_gridpoints, std::vector<unsigned int>(element_data.nb_elements, 0));
  output.number_densities.assign(nb_gridpoints, std::vector<double>(gas_phase.nb_species, 0.0));
  output.number_densities_cond.assign(nb_gridpoints, std::vector<double>(condensed_phase.nb_condensates, 0.0));
  output.element_cond_degree.assign(nb_gridpoints, std::vector<double>(element_data.nb_elements, 0.0));
  output.total_element_density.assign(nb_gridpoints, 0.0);
  output.nb_chemistry_iterations.assign(nb_gridpoints, 0);
  output.nb_cond_iterations.assign(nb_gridpoints, 0);
  output.nb_iterations.assign(nb_gridpoints, 0);
  output.mean_molecular_weight.assign(nb_gridpoints, 0.0);
  output.fastchem_flag.assign(nb_gridpoints, FASTCHEM_NO_CONVERGENCE);


  if (input.rainout_condensation)
  {
    rainoutCondensation(input, output);
  }
  else
  {
    #ifdef _OPENMP
    unsigned int nb_omp_threads = omp_get_max_threads();

    if (input.temperature.size() < nb_omp_threads)
      nb_omp_threads = input.temperature.size();


    omp_set_num_threads(nb_omp_threads);

    std::vector<FastChem<double_type>> fastchems(nb_omp_threads, *this);

    #pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int i=0; i<input.temperature.size(); i++)
    { //std::cout << i << " of " << input.temperature.size() << "  " << omp_get_thread_num() << "\n";
      if (!input.equilibrium_condensation)
      {
        output.fastchem_flag[i] = 
        fastchems[omp_get_thread_num()].calcDensity(
          input.temperature[i], 
          input.pressure[i]*1e6,
          false, 
          output.number_densities[i],
          output.total_element_density[i], 
          output.mean_molecular_weight[i],
          output.element_conserved[i],
          output.nb_chemistry_iterations[i]);

        output.nb_iterations[i] = output.nb_chemistry_iterations[i];
        output.nb_cond_iterations[i] = 0;
      }
      else
      {
        output.fastchem_flag[i] = 
        fastchems[omp_get_thread_num()].equilibriumCondensation(
          input.temperature[i],
          input.pressure[i]*1e6,
          output.number_densities[i],
          output.number_densities_cond[i],
          output.element_cond_degree[i],
          output.total_element_density[i], 
          output.mean_molecular_weight[i],
          output.element_conserved[i],
          output.nb_chemistry_iterations[i],
          output.nb_cond_iterations[i],
          output.nb_iterations[i]);
      }
    }
    #else
    for (unsigned int i=0; i<input.temperature.size(); i++)
    { 
      if (!input.equilibrium_condensation)
      {
        output.fastchem_flag[i] = calcDensity(
          input.temperature[i],
          input.pressure[i]*1e6,
          false, 
          output.number_densities[i],
          output.total_element_density[i], 
          output.mean_molecular_weight[i],
          output.element_conserved[i],
          output.nb_chemistry_iterations[i]);

        output.nb_iterations[i] = 0;
        output.nb_cond_iterations[i] = 0;
      }
      else
      {
        output.fastchem_flag[i] = equilibriumCondensation(
          input.temperature[i],
          input.pressure[i]*1e6,
          output.number_densities[i],
          output.number_densities_cond[i],
          output.element_cond_degree[i],
          output.total_element_density[i], 
          output.mean_molecular_weight[i],
          output.element_conserved[i],
          output.nb_chemistry_iterations[i],
          output.nb_cond_iterations[i],
          output.nb_iterations[i]);
      }
    }
    #endif
  }

  unsigned int status = *std::max_element(std::begin(output.fastchem_flag), std::end(output.fastchem_flag));

  is_busy = false;

  return status;
}



//Solve the chemistry for a single temperature and a single pressure
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
template <class double_type>
unsigned int FastChem<double_type>::calcDensity(
  const double temperature,
  const double pressure,
  const bool use_previous_solution,
  std::vector<double>& number_densities,
  double& total_element_density, 
  double& mean_molecular_weight,
  std::vector<unsigned int>& element_conserved,
  unsigned int& nb_chemistry_iterations)
{
  for (auto & i : gas_phase.molecules)  i.calcMassActionConstant(temperature, options.logK_limit);

  //this value will be fixed.
  double_type gas_density = pressure/(CONST_K * temperature);


  if (use_previous_solution == true)
  {
   //if we use the previous solution, convert the stored mixing ratios to number densities
   for (auto & i : gas_phase.species)  i->number_density *= gas_density;
  }
  else
  { 
    element_data.init(options.element_density_minlimit);

    //for a fresh start set all species to the minimum value
    for (auto & i : gas_phase.species) i->number_density = options.element_density_minlimit;
    
    //set the initial electron density to 1 (for stability reasons)
    if (element_data.e_ != FASTCHEM_UNKNOWN_SPECIES)
      element_data.elements[element_data.e_].number_density = 1.0;
  }


  //call the main FastChem solver  
  bool fastchem_converged = gas_phase.calculate(
    temperature, gas_density, nb_chemistry_iterations);


  if (!fastchem_converged && options.verbose_level >= 1) 
    std::cout << "Convergence problem in FastChem: Reached maximum number of chemistry iterations :(\n";


  //return output
  number_densities.assign(gas_phase.nb_species, 0.0);

  for (size_t i=0; i<gas_phase.nb_species; ++i)
    number_densities[i] = gas_phase.species[i]->number_density; 


  mean_molecular_weight = gas_phase.meanMolecularWeight(gas_density);
  total_element_density = gas_phase.totalElementDensity();


  for (auto & i : element_data.elements) 
    i.checkElementConservation(
      gas_phase.molecules,
      condensed_phase.condensates,
      total_element_density,
      options.element_conserve_accuracy);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    element_conserved[i] = element_data.elements[i].element_conserved;


  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!fastchem_converged) return_state = FASTCHEM_NO_CONVERGENCE;


  //store the mixing ratios in case we want to use them in the future
  for (auto & i : gas_phase.species) i->number_density /= gas_density;


  return return_state;
}



//calculates condensation using the rainout approximation, see Sect. 3.6 in Paper III
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
template <class double_type>
void FastChem<double_type>::rainoutCondensation(
  FastChemInput& input, FastChemOutput& output)
{
  std::vector<double> original_element_abundance = getElementAbundances();
  std::vector<double> original_element_epsilon(element_data.nb_elements, 0.0);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    original_element_epsilon[i] = element_data.elements[i].epsilon;


  for (unsigned int i=0; i<input.temperature.size(); i++)
  {
    output.fastchem_flag[i] = equilibriumCondensation(
      input.temperature[i],
      input.pressure[i]*1e6,
      output.number_densities[i],
      output.number_densities_cond[i],
      output.element_cond_degree[i],
      output.total_element_density[i], 
      output.mean_molecular_weight[i],
      output.element_conserved[i],
      output.nb_chemistry_iterations[i],
      output.nb_cond_iterations[i],
      output.nb_iterations[i]);

    //if the calculation at a point fails, 
    //we stop the entire calculation along the p-T structure
    if (output.fastchem_flag[i] != FASTCHEM_SUCCESS)
      break;


    //calculate the new effective element abundances in the gas phase
    std::vector<double> element_abundance_cond(element_data.nb_elements, 0.0);

    for (size_t j=0; j<element_data.nb_elements; ++j)
    {
      element_abundance_cond[j] = element_data.elements[j].phi;
      element_data.elements[j].degree_of_condensation = 0.0;

      output.element_cond_degree[i][j] = 
        (original_element_epsilon[j] - element_abundance_cond[j])/original_element_epsilon[j];

      if (original_element_epsilon[j] == 0 || output.element_cond_degree[i][j] < 0) 
        output.element_cond_degree[i][j] = 0;
    }

    setElementAbundances(element_abundance_cond);
    element_data.setRelativeAbundances();
    gas_phase.reInitialise();
    
    for (auto & j : condensed_phase.condensates)
      j.findReferenceElement(element_data.elements);
  }


  //set everything back to their original values
  setElementAbundances(original_element_abundance);
  element_data.setRelativeAbundances();
  gas_phase.reInitialise();
}



//Solve the equilibrium condensation for a single temperature and a single pressure
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
template <class double_type>
unsigned int FastChem<double_type>::equilibriumCondensation(
  const double temperature,
  const double pressure,
  std::vector<double>& number_densities,
  std::vector<double>& number_densities_cond,
  std::vector<double>& element_cond_degree,
  double& total_element_density, 
  double& mean_molecular_weight,
  std::vector<unsigned int>& element_conserved,
  unsigned int& nb_chem_iter,
  unsigned int& nb_cond_iter,
  unsigned int& nb_combined_iter)
{ 

  for (auto & i : gas_phase.molecules) i.calcMassActionConstant(temperature, options.logK_limit);
  for (auto & i : condensed_phase.condensates) i.calcMassActionConstant(temperature);

  //this value will be fixed.
  double_type gas_density = pressure/(CONST_K * temperature);

  //for a fresh start set all species to the minimum value
  for (auto & i : gas_phase.species) i->number_density = options.element_density_minlimit;
    
  //set the initial electron density to 1 (for stability reasons)
  if (element_data.e_ != FASTCHEM_UNKNOWN_SPECIES)
    element_data.elements[element_data.e_].number_density = 1.0;

  element_data.init(options.element_density_minlimit);

  for (auto & c : condensed_phase.condensates)
  {
    c.number_density = 0;
    c.degree_of_condensation = 0;
    c.activity_correction = 0;
    c.is_calculated = false;
  }


  nb_chem_iter = 0;
  nb_cond_iter = 0;

  unsigned int nb_iter = 0;

  //call the main FastChem solver  
  bool fastchem_converged = gas_phase.calculate(
    temperature, gas_density, nb_iter);

  nb_chem_iter += nb_iter;


  total_element_density = gas_phase.totalElementDensity();


  //search for potential condensates
  for (auto & i : condensed_phase.condensates)
  {
    i.calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);
    i.maxDensity(element_data.elements, total_element_density);
  }

  std::vector<Condensate<double_type>*> condensates_act;
  std::vector<Element<double_type>*> elements_cond;

  condensed_phase.selectActiveCondensates(condensates_act, elements_cond);


  //check for the phase rule
  if (elements_cond.size() == element_data.elements_wo_e.size())
    return FASTCHEM_PHASE_RULE_VIOLATION;


  bool cond_converged = false;
  bool combined_converged = false;

  if (condensates_act.size() > 0)
  { 
    std::vector<double_type> number_density_old(element_data.nb_elements, 0.0);

    for (size_t i=0; i<element_data.nb_elements; ++i)
      number_density_old[i] = element_data.elements[i].number_density;

    //run the equilibrium condensation and gas phase iteration
    for (nb_combined_iter=0; nb_combined_iter<options.nb_chem_cond_iter; ++nb_combined_iter)
    {
      condensed_phase.selectActiveCondensates(condensates_act, elements_cond);


      //check for the phase rule
      if (elements_cond.size() == element_data.elements_wo_e.size())
        return FASTCHEM_PHASE_RULE_VIOLATION;


      for (auto & i : condensates_act)
      {
        i->calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);
        i->maxDensity(element_data.elements, total_element_density);
      }

      cond_converged = condensed_phase.calculate(
        condensates_act,
        elements_cond,
        temperature,
        gas_density,
        total_element_density,
        gas_phase.molecules,
        nb_iter);

      nb_cond_iter += nb_iter;

      //gas_phase.reInitialise();
      
      fastchem_converged = gas_phase.calculate(
        temperature,
        gas_density,
        nb_iter);
   
      nb_chem_iter += nb_iter;

      total_element_density = gas_phase.totalElementDensity() + condensed_phase.totalElementDensity();

      for (auto & i : condensed_phase.condensates)
        i.calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);

      combined_converged = true;

      for (auto & i : element_data.elements)
      {
        if (std::fabs((i.number_density - number_density_old[i.index])) > options.chem_accuracy*number_density_old[i.index]*0.1)
          combined_converged = false;

        number_density_old[i.index] = i.number_density;
      }

      if (combined_converged && cond_converged) break;
    }


    //sanity check for the condensate activities
    for (auto & i : condensed_phase.condensates)
      if (i.log_activity > 0.001)
      { 
        cond_converged = false;
      }


    //remove condensates that are not present 
    //i.e. those with an activity smaller than 1
    for (auto & i : condensed_phase.condensates)
      if (i.log_activity < -0.01) i.number_density = 0.0;
    
    //and run the gas phase calculation one last time
    double_type phi_sum = 0;

    for (auto & i : element_data.elements)
    {
      i.calcDegreeOfCondensation(condensed_phase.condensates, total_element_density);
      phi_sum += i.phi;
    }

    for (auto & i : element_data.elements)
      i.normalisePhi(phi_sum);

    fastchem_converged = gas_phase.calculate(
      temperature,
      gas_density,
      nb_iter);
  }
  else
  { //no condensates stable...
    cond_converged = true;
    combined_converged = true;
  }


  if (!fastchem_converged && options.verbose_level >= 1) 
    std::cout << "Convergence problem in FastChem: Reached maximum number of chemistry iterations :(\n";

  if (!cond_converged&& options.verbose_level >= 1) 
    std::cout << "Convergence problem in FastChem: Equilibrium condensation calculation failed :(\n";

  if (!combined_converged && options.verbose_level >= 1) 
    std::cout << "Convergence problem in FastChem: Combined gas-phase & equilibrium condensation calculation failed :(\n";


  //return output
  number_densities.assign(gas_phase.nb_species, 0.0);

  for (size_t i=0; i<gas_phase.nb_species; ++i)
    number_densities[i] = gas_phase.species[i]->number_density; 

  number_densities_cond.assign(condensed_phase.nb_condensates, 0.0);

  for (size_t i=0; i<condensed_phase.nb_condensates; ++i)
    number_densities_cond[i] = condensed_phase.condensates[i].number_density;
  
  element_cond_degree.assign(element_data.nb_elements, 0.0);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    element_cond_degree[i] = element_data.elements[i].degree_of_condensation;

  mean_molecular_weight = gas_phase.meanMolecularWeight(gas_density);
  total_element_density = gas_phase.totalElementDensity() + condensed_phase.totalElementDensity();


  for (auto & i : element_data.elements) 
    i.checkElementConservation(
      gas_phase.molecules,
      condensed_phase.condensates,
      total_element_density,
      options.element_conserve_accuracy);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    element_conserved[i] = element_data.elements[i].element_conserved;


  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!fastchem_converged || !cond_converged || !combined_converged) 
    return_state = FASTCHEM_NO_CONVERGENCE;


  return return_state;
}



template class FastChem<double>;
template class FastChem<long double>;
}


