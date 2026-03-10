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


#ifndef _cond_phase_solver_h
#define _cond_phase_solver_h

#include <vector>
#include <iostream>
#include <string>

#include "../fastchem_constants.h"
#include "../species_struct.h"
#include "../options.h"

#include "../../_ext/Eigen/Dense"

namespace Eigen {

using MatrixXdt = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

using VectorXdt = Matrix<double, Eigen::Dynamic, 1>;
}


namespace fastchem {


//Solver class
class CondPhaseSolver{
  public:
    CondPhaseSolver(FastChemOptions& options_)
      : options(options_)
      {}

    bool newtonStep(
      const std::vector<Condensate*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<Element*>& elements,
      const std::vector<Molecule>& molecules,
      const double total_element_density,
      const std::vector<double>& cond_densities,
      const std::vector<double>& element_densities,
      const std::vector<double>& activity_corr,
      const std::vector<double>& log_cond_densities,
      const std::vector<double>& log_activity_corr,
      Eigen::VectorXdt& result,
      Eigen::VectorXdt& scaling_factors,
      double& objective_function);

    bool newtonStepFull(
      const std::vector<Condensate*>& condensates,
      const std::vector<Element*>& elements,
      const std::vector<Molecule>& molecules,
      const double total_element_density,
      const std::vector<double>& cond_densities,
      const std::vector<double>& element_densities,
      const std::vector<double>& activity_corr,
      const std::vector<double>& log_cond_densities,
      const std::vector<double>& log_activity_corr,
      Eigen::VectorXdt& result,
      Eigen::VectorXdt& scaling_factors,
      double& objective_function);

    double assembleRightHandSide(
      const std::vector<Condensate*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double>& activity_corr,
      const std::vector<double>& number_densities,
      const std::vector<double>& log_number_densities,
      const std::vector<double>& log_activity_corr,
      const std::vector< Element* >& elements,
      const std::vector< Molecule >& molecules,
      const double total_element_density,
      const Eigen::VectorXdt& scaling_factors,
      Eigen::VectorXdt& rhs);

    double assembleRightHandSideFull(
      const std::vector<Condensate*>& condensates,
      const std::vector<double>& activity_corr,
      const std::vector<double>& number_densities,
      const std::vector<double>& log_number_densities,
      const std::vector<double>& log_activity_corr,
      const std::vector< Element* >& elements,
      const std::vector< Molecule >& molecules,
      const double total_element_density,
      const Eigen::VectorXdt& scaling_factors,
      Eigen::VectorXdt& rhs);

    void resetLM();
    void adaptLM(const double objective_function);

    double backtrackStep(
      const double objective_function_0,
      const double objective_function_prev,
      const double objective_function_2prev,
      const double lambda_prev,
      const double lambda_2prev);

    double objectiveFunction(
      const std::vector<Condensate*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double>& activity_corr,
      const std::vector<double>& number_densities,
      const std::vector<double>& log_number_densities,
      const std::vector<double>& log_activity_corr,
      const std::vector< Element* >& elements,
      const std::vector< Molecule >& molecules,
      const double total_element_density,
      const Eigen::VectorXdt& scaling_factors);
  private:
    FastChemOptions& options;

    //Scratch Eigen objects — pre-allocated once, reused across Newton iterations
    Eigen::MatrixXdt scratch_jacobian_;
    Eigen::VectorXdt scratch_rhs_;

    double lm_mu = 0.01;
    double lm_objective_prev = 0;
    bool lm_has_prev_objective = false;

     Eigen::MatrixXdt assemblePerturbedHessian(
      const Eigen::MatrixXdt& jacobian,
      const double perturbation);

    bool solveSystem(
      const Eigen::MatrixXdt& jacobian,
      const Eigen::VectorXdt& rhs,
      Eigen::VectorXdt& result);

    Eigen::VectorXdt assembleJacobian(
      const std::vector<Condensate*>& condensates,
      const std::vector<double>& activity_corr,
      const std::vector<double>& number_densities,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<Element*>& elements,
      const std::vector<Molecule>& molecules,
      const double total_element_density,
      Eigen::MatrixXdt& jacobian);

    Eigen::VectorXdt assembleJacobianFull(
      const std::vector<Condensate*>& condensates,
      const std::vector<double>& activity_corr,
      const std::vector<double>& number_densities,
      const std::vector<Element*>& elements,
      const std::vector<Molecule>& molecules,
      const double total_element_density,
      Eigen::MatrixXdt& jacobian);
};



}
#endif
