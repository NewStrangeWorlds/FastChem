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


template <class double_type>
class CondensedPhase {
  public:
    CondensedPhase(
      FastChemOptions<double_type>& options_,
      ElementData<double_type>& element_data_);
    CondensedPhase(
      const CondensedPhase &obj,
      FastChemOptions<double_type>& options_,
      ElementData<double_type>& element_data_);

    std::vector< Condensate<double_type> > condensates;

    void selectActiveCondensates(
      std::vector<Condensate<double_type>*>& condensates_act,
      std::vector<Element<double_type>*>& elements_cond);
    
    size_t nb_condensates = 0;
    size_t nb_elements = 0;

    bool is_initialised = false;

    bool calculate(
      std::vector<Condensate<double_type>*>& condensates_act,
      std::vector<Element<double_type>*>& elements_cond,
      const double temperature,
      const double density,
      const double total_element_density,
      std::vector<Molecule<double_type>>& molecules,
      unsigned int& nb_iterations);

    double totalElementDensity();
  private:
    FastChemOptions<double_type>& options;
    ElementData<double_type>& element_data;
    std::vector<Element<double_type>>& elements;

    CondPhaseSolver<double_type> solver;

    void init();

    bool readCondensateData(const std::string& species_data_file);
    void addCondensate(
      const std::string name,
      const std::string symbol,
      const std::vector<std::string> species_elements,
      const std::vector<int> stoichiometric_coeff,
      const std::string phase,
      const std::vector<double>& fit_coeff_limits,
      const std::vector<std::vector<double_type>>& fit_coeff);

    void selectJacobianCondensates(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& number_density_cond,
      const std::vector<double_type>& activity_corr,
      std::vector<unsigned int>& condensates_jac,
      std::vector<unsigned int>& condensates_rem);

    void selectJacobianCondensates2(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& number_density_cond,
      const std::vector<double_type>& activity_corr,
      std::vector<unsigned int>& condensates_jac,
      std::vector<unsigned int>& condensates_rem,
      Eigen::MatrixXdt<double_type>& jacobian);

    double_type correctValues(
      const Eigen::VectorXdt<double_type>& result,
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double_type>& activity_corr_old,
      std::vector<double_type>& activity_corr_new,
      const std::vector<double_type>& cond_number_dens_old,
      std::vector<double_type>& cond_number_dens_new,
      const std::vector<Element<double_type>*>& elements,
      const std::vector<double_type>& elem_number_dens_old,
      std::vector<double_type>& elem_number_dens_new,
      const double max_change);

    double_type correctValuesFull(
      const Eigen::VectorXdt<double_type>& result,
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& activity_corr_old,
      std::vector<double_type>& activity_corr_new,
      const std::vector<double_type>& cond_number_dens_old,
      std::vector<double_type>& cond_number_dens_new,
      const std::vector<Element<double_type>*>& elements,
      const std::vector<double_type>& elem_number_dens_old,
      std::vector<double_type>& elem_number_dens_new);

    double_type newtonBacktrack(
      const double_type objective_function_0,
      const Eigen::VectorXdt<double_type>& result,
      const Eigen::VectorXdt<double_type>& scaling_factors,
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double_type>& activity_corr_old,
      std::vector<double_type>& activity_corr_new,
      const std::vector<double_type>& cond_number_dens_old,
      std::vector<double_type>& cond_number_dens_new,
      const std::vector<Element<double_type>*>& elements,
      const std::vector<double_type>& elem_number_dens_old,
      std::vector<double_type>& elem_number_dens_new,
      std::vector<Molecule<double_type>>& molecules,
      const double total_element_density,
      const double temperature,
      const double max_change);
};



}

#endif
