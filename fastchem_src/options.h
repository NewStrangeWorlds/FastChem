/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2024 Daniel Kitzmann, Joachim Stock
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


#ifndef _options_h
#define _options_h

#include <vector>
#include <iostream>
#include <string>
#include <map>


namespace fastchem {


enum class ParameterFloat {
    invalid_parameter,
    cond_tau,
    cond_limit_change,
    chem_accuracy,
    element_conserve_accuracy,
    newton_err,
    cond_accuracy,
    element_minlimit,
    molecule_minlimit,
    cond_trace_density_threshold,
    chem_conv_mixing_ratio_floor,
    newton_lm_lambda
};


enum class ParameterInt {
    invalid_parameter,
    nb_max_cond_iter,
    nb_max_comb_iter,
    nb_max_fastchem_iter,
    nb_max_bisection_iter,
    nb_max_newton_iter,
    nb_max_neldermead_iter,
    nb_switch_to_newton,
    nb_switch_to_joint
};


enum class ParameterBool {
    invalid_parameter,
    cond_reduce_system_size,
    cond_use_svd,
    cond_use_data_validity_limits,
    cond_use_lm
};


struct FastChemOptions{
  FastChemOptions(
    const std::string& parameter_file,
    const unsigned int verbose_level_start)
      : verbose_level(verbose_level_start) {
          parameter_file_loaded = readParameterFile(parameter_file);
          init();
      };
  FastChemOptions(
    const std::string& element_abundances_file_,
    const std::string& species_data_file_,
    const std::string& cond_species_data_file_,
    const unsigned int verbose_level_start)
      : verbose_level(verbose_level_start)
      , element_abundances_file(element_abundances_file_)
      , species_data_file(species_data_file_)
      , condensates_data_file(cond_species_data_file_) {
          init();
      };

  void init();

  unsigned int nb_max_fastchem_iter = 3000;
  unsigned int nb_max_bisection_iter = 3000;
  unsigned int nb_max_newton_iter = 3000;
  unsigned int nb_max_neldermead_iter = 3000;
  unsigned int nb_max_cond_iter = 3000;
  unsigned int nb_max_comb_iter = 3000;
  unsigned int nb_switch_to_newton = 400;
  unsigned int nb_switch_to_joint = 3000;

  double chem_accuracy = 1e-5;
  double newton_err = 1e-5;
  double cond_accuracy = 1e-5;
  double element_conserve_accuracy = 1e-4;
  //Levenberg-Marquardt regularisation added to the diagonal of the (row-scaled)
  //gas-phase Newton Jacobian. Keeps the multidimensional Newton accelerator stable when
  //the Jacobian becomes near-singular for strongly non-solar / rained-out compositions.
  double newton_lm_lambda = 1e-4;
  
  //smallest allowed particle number densities
  double element_density_minlimit = 1e-155;
  double molecule_density_minlimit = 1e-155;
  //species with a mixing ratio (n_i / n_gas) below this value are ignored by the
  //gas-phase convergence test. They carry no mass, but high-order clusters (e.g. C5,
  //C2H4) are proportional to a power of a trace parent density and so amplify its
  //numerical noise; they can oscillate forever and prevent convergence even though
  //they are physically zero. The value is far below any abundance that can affect
  //element conservation, so it does not change the converged result.
  double chem_conv_mixing_ratio_floor = 1e-150;
  //condensates with a max_number_density below this are skipped in the Newton solver
  //(number_density set to 0) to avoid Jacobian ill-conditioning from near-depleted trace
  //elements after rainout
  double condensate_density_threshhold = 1e-100;
  //a condensate whose maximum possible density (limited by its most depleted
  //constituent element) falls below this fraction of the total element density is
  //treated as a "ghost": its limiting element has effectively rained out, so the
  //condensate cannot form. Such condensates are excluded from the active set and the
  //convergence checks. This avoids feeding a degenerate, ill-conditioned system
  //(element densities ~1e-60 against a gas of ~1e18) into the Newton solver, which
  //otherwise produces NaNs in rainout calculations with strongly depleted elements.
  double cond_limiting_density_ratio = 1e-20;

  unsigned int verbose_level = 1;

  bool cond_use_data_validity_limits = true;

  bool chem_use_backup_solver = false;

  bool cond_use_svd = false;
  bool cond_reduce_system_size = true;
  bool cond_use_lm = true;
  double cond_iter_change_limit = 5;
  double cond_tau = 1e-15;

  std::string chemical_element_file = "";
  std::string element_abundances_file = "";
  std::string species_data_file = "";
  std::string condensates_data_file = "";

  bool parameter_file_loaded = false;
  bool readParameterFile(const std::string& model_parameter_file);

  ParameterFloat resolveParameter(const std::string& parameter);
  ParameterBool resolveParameterBool(const std::string& parameter);
  ParameterInt resolveParameterInt(const std::string& parameter);
};


}
#endif
