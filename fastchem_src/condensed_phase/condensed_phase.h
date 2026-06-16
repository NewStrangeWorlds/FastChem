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

#ifndef _condensed_phase_h
#define _condensed_phase_h

#include <string>
#include <vector>

#include "../elements/elements.h"
#include "../options.h"
#include "../species_struct.h"
#include "solver.h"


namespace fastchem {


class CondensedPhase {
  public:
    CondensedPhase(
      FastChemOptions& options_,
      ElementData& element_data_);
    CondensedPhase(
      const CondensedPhase &obj,
      FastChemOptions& options_,
      ElementData& element_data_);

    std::vector< Condensate > condensates;

    void selectActiveCondensates(
      std::vector<Condensate*>& condensates_act,
      std::vector<Element*>& elements_cond,
      const double total_element_density);

    bool isGhostCondensate(
      const Condensate& condensate,
      const double total_element_density) const;

    size_t nb_condensates = 0;
    size_t nb_elements = 0;

    bool is_initialised = false;

    bool calculate(
      std::vector<Condensate*>& condensates_act,
      std::vector<Element*>& elements_cond,
      const double temperature,
      const double density,
      const double total_element_density,
      std::vector<Molecule>& molecules,
      unsigned int& nb_iterations);

    double totalElementDensity();
  private:
    FastChemOptions& options;
    ElementData& element_data;
    std::vector<Element>& elements;

    CondPhaseSolver solver;

    void init();

    bool readCondensateData(const std::string& species_data_file);
    void addCondensate(
      const std::string name,
      const std::string symbol,
      const std::vector<std::string> species_elements,
      const std::vector<int> stoichiometric_coeff,
      const std::string phase,
      const std::vector<double>& fit_coeff_limits,
      const std::vector<std::vector<double>>& fit_coeff);

    void selectJacobianCondensates(
      const std::vector<Condensate*>& condensates,
      const std::vector<double>& number_density_cond,
      const std::vector<double>& activity_corr,
      std::vector<unsigned int>& condensates_jac,
      std::vector<unsigned int>& condensates_rem);

    double correctValues(
      const Eigen::VectorXdt& result,
      const std::vector<Condensate*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double>& log_activity_corr_old,
      std::vector<double>& log_activity_corr_new,
      const std::vector<double>& log_cond_dens_old,
      std::vector<double>& log_cond_dens_new,
      const std::vector<Element*>& elements,
      const std::vector<double>& log_elem_dens_old,
      std::vector<double>& log_elem_dens_new,
      const double max_change);

    double correctValuesFull(
      const Eigen::VectorXdt& result,
      const std::vector<Condensate*>& condensates,
      const std::vector<double>& log_activity_corr_old,
      std::vector<double>& log_activity_corr_new,
      const std::vector<double>& log_cond_dens_old,
      std::vector<double>& log_cond_dens_new,
      const std::vector<Element*>& elements,
      const std::vector<double>& log_elem_dens_old,
      std::vector<double>& log_elem_dens_new);
};



}

#endif
