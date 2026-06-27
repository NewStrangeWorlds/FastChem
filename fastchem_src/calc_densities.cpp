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
#include <iomanip>


#ifdef _OPENMP
#include <omp.h>
#endif

#include "input_output_struct.h"


namespace fastchem {


unsigned int FastChem::calcDensities(
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


    if (thread_copies_.size() != nb_omp_threads)
    {
      thread_copies_.clear();
      thread_copies_.reserve(nb_omp_threads);
      for (unsigned int j = 0; j < nb_omp_threads; ++j)
        thread_copies_.emplace_back(*this);
    }

    #pragma omp parallel for schedule(dynamic, 1) num_threads(nb_omp_threads)
    for (unsigned int i=0; i<input.temperature.size(); i++)
    { std::cout << i << " of " << input.temperature.size() << "  " << omp_get_thread_num() << "\n";
      if (!input.equilibrium_condensation)
      {
        output.fastchem_flag[i] =
        thread_copies_[omp_get_thread_num()].calcDensity(
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
        thread_copies_[omp_get_thread_num()].equilibriumCondensation(
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
unsigned int FastChem::calcDensity(
  const double temperature,
  const double pressure,
  const bool use_previous_solution,
  std::vector<double>& number_densities,
  double& total_element_density, 
  double& mean_molecular_weight,
  std::vector<unsigned int>& element_conserved,
  unsigned int& nb_chemistry_iterations)
{
  for (auto & i : gas_phase.molecules)  i.calcMassActionConstant(temperature);

  //this value will be fixed.
  double gas_density = pressure/(CONST_K * temperature);


  if (use_previous_solution == true)
  {
   //if we use the previous solution, convert the stored mixing ratios to number densities
   for (auto & i : gas_phase.species)
   {
     i->number_density *= gas_density;
     i->log_number_density = safeLog(i->number_density);
   }
  }
  else
  {
    element_data.init(options.element_density_minlimit);

    //for a fresh start set all species to the minimum value
    const double log_min = std::log(options.element_density_minlimit);
    for (auto & i : gas_phase.species)
    {
      i->number_density = options.element_density_minlimit;
      i->log_number_density = log_min;
    }

    //set the initial electron density to 1 (for stability reasons)
    if (element_data.e_ != FASTCHEM_UNKNOWN_SPECIES)
    {
      element_data.elements[element_data.e_].number_density = 1.0;
      element_data.elements[element_data.e_].log_number_density = 0.0;
    }
  }


  //call the main FastChem solver  
  bool fastchem_converged = gas_phase.calculate(
    temperature, 
    gas_density, 
    nb_chemistry_iterations);


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
  for (auto & i : gas_phase.species)
  {
    i->number_density /= gas_density;
    i->log_number_density = safeLog(i->number_density);
  }


  return return_state;
}



//calculates condensation using the rainout approximation, see Sect. 3.6 in Paper III
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
void FastChem::rainoutCondensation(
  FastChemInput& input, FastChemOutput& output)
{
  std::vector<double> original_element_abundance = getElementAbundances();
  std::vector<double> original_element_epsilon(element_data.nb_elements, 0.0);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    original_element_epsilon[i] = element_data.elements[i].epsilon;


  for (unsigned int i=0; i<input.temperature.size(); i++)
  { std::cout << i << " of " << input.temperature.size() << "\n";
    //warm-start every grid point (except the first) from the previous point's
    //converged solution; in a rainout calculation the points are processed in sequence
    //and adjacent layers have nearly identical chemistry.
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
      output.nb_iterations[i],
      i > 0);

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


void FastChem::updatePhi(double total_element_density)
{
  //Guard against invalid total (can occur if Newton diverged and produced NaN)
  if (!std::isfinite(total_element_density) || total_element_density <= 0)
    return;

  double phi_sum = 0;
  for (auto & i : element_data.elements)
  {
    i.calcDegreeOfCondensation(condensed_phase.condensates, total_element_density);
    phi_sum += i.phi;
  }

  if (!std::isfinite(phi_sum) || phi_sum <= 0)
    return;

  for (auto & i : element_data.elements)
    i.normalisePhi(phi_sum);
}


//Solve the equilibrium condensation for a single temperature and a single pressure
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
unsigned int FastChem::equilibriumCondensation(
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
  unsigned int& nb_combined_iter,
  const bool use_previous_solution)
{

  for (auto & i : gas_phase.molecules) i.calcMassActionConstant(temperature);
  for (auto & i : condensed_phase.condensates) i.calcMassActionConstant(temperature);

  //this value will be fixed.
  double gas_density = pressure/(CONST_K * temperature);
  
  
  // if (use_previous_solution)
  // {
  //   //Warm start: reuse the converged solution of the previous (adjacent) grid point.
  //   //Between grid points the densities are kept as mixing ratios (see the end of this
  //   //function), so they are converted back to number densities for the new gas density.
  //   //Adjacent layers have nearly identical chemistry, which makes this a far better
  //   //initial guess than the density floor and lets the gas-phase solver converge
  //   //compositions (e.g. strongly oxygen-dominated, hydrogen-poor rainout layers) that
  //   //it cannot reach from a cold start.
  //   //
  //   //Note: this must restore the densities WITHOUT going through element_data.init(),
  //   //which would reset every element (and the electron) back to the density floor.
  //   //Only the per-grid-point condensation bookkeeping is reset here.
  //   for (auto & i : gas_phase.species)
  //   {
  //     i->number_density *= gas_density;
  //     i->log_number_density = safeLog(i->number_density);
  //   }

  //   for (auto & e : element_data.elements)
  //   {
  //     e.degree_of_condensation = 0.0;
  //     e.phi = e.epsilon;
  //     e.fixed_by_condensation = false;
  //   }
  // }
  // else
  {
    //reset the per-element condensation state and set all element densities to the
    //minimum value
    element_data.init(options.element_density_minlimit);

    //for a fresh start set all species to the minimum value
    const double log_min = std::log(options.element_density_minlimit);
    for (auto & i : gas_phase.species)
    {
      i->number_density = options.element_density_minlimit;
      i->log_number_density = log_min;
    }

    //set the initial electron density to 1 (for stability reasons)
    if (element_data.e_ != FASTCHEM_UNKNOWN_SPECIES)
    {
      element_data.elements[element_data.e_].number_density = 1.0;
      element_data.elements[element_data.e_].log_number_density = 0.0;
    }
  }

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
    temperature,
    gas_density,
    nb_iter);

  nb_chem_iter += nb_iter;


  total_element_density = gas_phase.totalElementDensity();

  //search for potential condensates
  for (auto & i : condensed_phase.condensates)
  {
    i.calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);
    i.maxDensity(element_data.elements, total_element_density);
  }

  std::vector<Condensate*> condensates_act;
  std::vector<Element*> elements_cond;
  
  condensed_phase.selectActiveCondensates(condensates_act, elements_cond, total_element_density);

  bool cond_converged = false;
  bool combined_converged = false;

  bool chem_backup_solver_default = options.chem_use_backup_solver;
  size_t nb_condensed_elements = 0;

  if (options.verbose_level >= 4 && condensates_act.size() > 0)
  {
    std::cout << "  T=" << temperature << " Initial active condensates: " << condensates_act.size() << "\n";
    for (auto & c : condensates_act)
      std::cout << "    " << c->symbol << "  log_act=" << c->log_activity << "\n";
    for (auto & e : elements_cond)
      std::cout << "    elem " << e->symbol << "  eps=" << e->epsilon << "  n=" << e->number_density << "\n";
    std::cout << "  total_element_density=" << total_element_density << "\n";
  }

  if (condensates_act.size() > 0)
  {
    options.chem_use_backup_solver = true;

    std::vector<double> log_density_old(element_data.nb_elements, static_cast<double>(LOG_DENSITY_FLOOR));

    for (size_t i=0; i<element_data.nb_elements; ++i)
      log_density_old[i] = element_data.elements[i].log_number_density;

    //Early termination tracking for divergence and stagnation
    double prev_max_change = std::numeric_limits<double>::max();
    unsigned int nb_growing = 0;
    const unsigned int max_growing = 100;

    double checkpoint_max_change = std::numeric_limits<double>::max();
    unsigned int nb_stagnant_checkpoints = 0;
    const unsigned int checkpoint_interval = 500;
    const unsigned int max_stagnant_checkpoints = 2;

    //run the equilibrium condensation and gas phase iteration
    for (nb_combined_iter=0; nb_combined_iter<options.nb_max_comb_iter; ++nb_combined_iter)
    {
      double total_element_density_old = total_element_density;
      condensed_phase.selectActiveCondensates(condensates_act, elements_cond, total_element_density);

      for (auto & i : condensates_act)
      {
        i->calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);
        i->maxDensity(element_data.elements, total_element_density);
      }

      if (nb_combined_iter >= options.nb_switch_to_joint && condensates_act.size() > 0)
      {
        // Joint Newton step: couples gas-phase element densities and condensate densities
        jointNewtonStep(condensates_act, temperature, gas_density, total_element_density);

        // Recompute activities for convergence check
        for (auto & i : condensed_phase.condensates)
          i.calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);
      }
      else
      {
        cond_converged = condensed_phase.calculate(
          condensates_act,
          elements_cond,
          temperature,
          gas_density,
          total_element_density,
          gas_phase.molecules,
          nb_iter);

        nb_cond_iter += nb_iter;

        for (auto & e : element_data.elements)
        {
          bool is_condensed = false;

          for (auto & c : e.condensate_list)
          {
            if (condensed_phase.condensates[c].log_activity > LOG_ACTIVITY_THRESHOLD)
            {
              is_condensed = true;
              break;
            }
          }

          if (!is_condensed)
            e.fixed_by_condensation = false;
        }
        
        //Recompute phi values after the condensed phase has updated condensate densities,
        //so that the gas-phase solver below uses the current condensate state.
        updatePhi(total_element_density);
        
        fastchem_converged = gas_phase.calculate(
          temperature,
          gas_density,
          nb_iter);

        nb_chem_iter += nb_iter;

        total_element_density = gas_phase.totalElementDensity() + condensed_phase.totalElementDensity();

        for (auto & i : condensed_phase.condensates)
          i.calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);
      }

      //sanity check for the condensate activities
      //Skip condensates with negligible max density (trace elements after rainout):
      //these are handled separately in selectActiveCondensates and their activity
      //is not a reliable convergence indicator.
      for (auto & i : condensed_phase.condensates)
        if (i.log_activity > 0.001
            && !condensed_phase.isGhostCondensate(i, total_element_density))
        {
          cond_converged = false;
        }

      combined_converged = true;

      double rel_total_delta = std::fabs(total_element_density - total_element_density_old)/total_element_density_old;

      if (rel_total_delta > options.chem_accuracy)
        combined_converged = false;

      double max_change = 0;
      for (auto & i : element_data.elements)
      {
        double change = std::fabs(i.log_number_density - log_density_old[i.index]);
        if (change > max_change) max_change = change;
        if (change > options.chem_accuracy)
          combined_converged = false;
        log_density_old[i.index] = i.log_number_density;
      }

      if (combined_converged && cond_converged) break;

      //Divergence detection: max_change growing for too many consecutive iterations
      if (max_change > prev_max_change)
        nb_growing++;
      else
        nb_growing = 0;
      prev_max_change = max_change;

      if (nb_growing >= max_growing)
      {
        if (options.verbose_level >= 2)
          std::cout << "  Combined iteration diverging at iter " << nb_combined_iter
                    << " (max_change=" << max_change << " growing for " << nb_growing << " iters)\n";
        break;
      }

      //Stagnation detection: no improvement over checkpoint intervals
      if (nb_combined_iter > 0 && nb_combined_iter % checkpoint_interval == 0)
      {
        if (max_change >= checkpoint_max_change)
          nb_stagnant_checkpoints++;
        else
          nb_stagnant_checkpoints = 0;
        checkpoint_max_change = max_change;

        if (nb_stagnant_checkpoints >= max_stagnant_checkpoints)
        {
          if (options.verbose_level >= 2)
            std::cout << "  Combined iteration stagnant at iter " << nb_combined_iter
                      << " (max_change=" << max_change << " for " << nb_stagnant_checkpoints << " checkpoints)\n";
          break;
        }
      }
    }


    //sanity check for the condensate activities
    for (auto & i : condensed_phase.condensates)
      if (i.log_activity > 0.001
          && !condensed_phase.isGhostCondensate(i, total_element_density))
      {
        cond_converged = false;
      }

    //remove condensates that are not present
    //i.e. those with an activity smaller than 1.
    //Ghost condensates (limiting element rained out) are never activated, so their
    //density is already zero and they are left untouched here.
    for (auto & i : condensed_phase.condensates)
      if (i.log_activity < LOG_ACTIVITY_THRESHOLD
          && !condensed_phase.isGhostCondensate(i, total_element_density))
        i.number_density = 0.0;

    //recompute total_element_density after zeroing so that the subsequent phi
    //computation is consistent with the condensate densities actually present
    total_element_density = gas_phase.totalElementDensity() + condensed_phase.totalElementDensity();
    
    //and run the gas phase calculation one last time
    updatePhi(total_element_density);

    nb_condensed_elements = 0;
    for (auto & i : element_data.elements)
    {
      if (i.degree_of_condensation == 0)
        i.fixed_by_condensation = false;
      else
        nb_condensed_elements++;
    }

    fastchem_converged = gas_phase.calculate(
      temperature,
      gas_density,
      nb_iter);

    //The alternating gas-condensation iteration converges on per-step changes and can leave
    //a small element-conservation bias for stiff, strongly-condensing mixtures. The joint
    //gas-condensate Newton couples the gas-phase element densities and the condensate
    //densities and removes it - but it is more expensive than the standard solver, so it is
    //only invoked when the cheap final gas-phase solve above actually left elements
    //unconserved. Inert elements (noble gases: no molecules and no condensate) are excluded
    //from this test - they are corrected exactly further below and would otherwise trigger
    //the joint solver needlessly.
    total_element_density = gas_phase.totalElementDensity()
                          + condensed_phase.totalElementDensity();

    bool needs_joint = false;
    for (auto & e : element_data.elements)
    {
      if (e.symbol == "e-" || e.epsilon <= 0) continue;
      if (e.molecule_list.empty() && e.condensate_list.empty()) continue;

      double a = e.number_density;
      for (auto & m : e.molecule_list)
        a += gas_phase.molecules[m].stoichiometric_vector[e.index]
             * gas_phase.molecules[m].number_density;
      for (auto & c : e.condensate_list)
        a += condensed_phase.condensates[c].stoichiometric_vector[e.index]
             * condensed_phase.condensates[c].number_density;

      if (std::fabs(a/(total_element_density*e.epsilon) - 1.0)
            > options.element_conserve_accuracy)
      { needs_joint = true; break; }
    }

    if (needs_joint)
    {
      std::vector<double> log_density_joint(element_data.nb_elements,
                                            static_cast<double>(LOG_DENSITY_FLOOR));

      bool joint_converged = false;
      for (unsigned int k=0; k<options.nb_max_comb_iter; ++k)
      {
        for (auto & e : element_data.elements)
          log_density_joint[e.index] = e.log_number_density;

        jointNewtonStep(condensates_act, temperature, gas_density, total_element_density);

        for (auto & i : condensed_phase.condensates)
          i.calcActivity(temperature, element_data.elements, options.cond_use_data_validity_limits);

        total_element_density = gas_phase.totalElementDensity()
                              + condensed_phase.totalElementDensity();

        double max_change = 0;
        for (auto & e : element_data.elements)
          max_change = std::max(max_change,
            std::fabs(e.log_number_density - log_density_joint[e.index]));

        if (max_change < options.chem_accuracy) { joint_converged = true; break; }
      }

      //if the joint Newton reached a fixed point, the coupled system is converged
      if (joint_converged)
      {
        combined_converged = true;
        cond_converged = true;
      }
    }
  }
  else
  { //no condensates stable...
    cond_converged = true;
    combined_converged = true;
  }


  options.chem_use_backup_solver = chem_backup_solver_default;


  if (!fastchem_converged && options.verbose_level >= 1) 
    std::cout << "Convergence problem in FastChem: Reached maximum number of chemistry iterations :(\n";

  if (!cond_converged && options.verbose_level >= 1)
  {
    std::cout << "Convergence problem in FastChem: Equilibrium condensation calculation failed :(\n";
    std::cout << "  T=" << temperature << " P=" << gas_density << "\n";
    std::cout << "  Active condensates with log_activity > 0.001:\n";
    for (auto & i : condensed_phase.condensates)
      if (i.log_activity > 0.001)
        std::cout << "    " << i.symbol << "  log_act=" << i.log_activity
                  << "  n=" << i.number_density << "\n";
    std::cout << "  Elements involved:\n";
    for (auto & i : element_data.elements)
      if (i.symbol != "e-" && i.epsilon > 0)
        std::cout << "    " << i.symbol << "  eps=" << i.epsilon
                  << "  phi=" << i.phi << "  n=" << i.number_density
                  << "  DOC=" << i.degree_of_condensation << "\n";
  }

  if (!combined_converged && options.verbose_level >= 1) 
    std::cout << "Convergence problem in FastChem: Combined gas-phase & equilibrium condensation calculation failed :(\n";


  //Enforce exact conservation for inert elements, i.e. those that form no molecules
  //and enter no condensate (e.g. the noble gases). For such an element the entire
  //abundance is in the free atom, so the exact result is n = eps * N_total. The
  //gas-phase solve fixes such an element solely through the gas-phase reference (phi
  //together with the major/minor and degree-of-condensation bookkeeping), which can be
  //marginally inconsistent with the atom-based total budget. For an element that forms
  //molecules this small discrepancy is absorbed by its atomic/molecular partition, but
  //an inert element has no such buffer and is left slightly off; because the total
  //budget is fixed, that error is then redistributed as a (uniform) conservation offset
  //onto every other element. The fixed-point passes make the correction self-consistent
  //with the total it feeds back into.
  {
    auto is_inert = [](const Element& e)
      { return e.symbol != "e-" && e.epsilon > 0
               && e.molecule_list.empty() && e.condensate_list.empty(); };

    bool any_inert = false;
    for (auto & e : element_data.elements) if (is_inert(e)) { any_inert = true; break; }

    if (any_inert)
      for (int pass=0; pass<3; ++pass)
      {
        const double tot = gas_phase.totalElementDensity() + condensed_phase.totalElementDensity();
        for (auto & e : element_data.elements)
          if (is_inert(e))
          {
            e.number_density = e.epsilon * tot;
            e.log_number_density = std::log(e.number_density);
          }
      }
  }


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

  //check for the phase rule
  if (nb_condensed_elements == element_data.elements_wo_e.size())
    return_state = FASTCHEM_PHASE_RULE_VIOLATION;

  //store the gas-phase densities as mixing ratios so that the next (adjacent) grid
  //point can warm-start from them (see the use_previous_solution branch above). This
  //must come after all output and conservation quantities have been computed from the
  //absolute number densities.
  for (auto & i : gas_phase.species)
  {
    i->number_density /= gas_density;
    i->log_number_density = safeLog(i->number_density);
  }

  return return_state;
}



}


