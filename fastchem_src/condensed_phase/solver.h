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

template <class double_type>
using MatrixXdt = Matrix<double_type, Eigen::Dynamic, Eigen::Dynamic>;

template <class double_type>
using VectorXdt = Matrix<double_type, Eigen::Dynamic, 1>;
}


namespace fastchem {


//Solver class
template <class double_type>
class CondPhaseSolver{
  public:
    CondPhaseSolver(FastChemOptions<double_type>& options_)
      : options(options_)
      {}
    /*Eigen::MatrixXdt<double_type> assembleJacobian(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& activity_corr,
      const size_t nb_condensates_removed,
      const std::vector<Element<double_type>*>& elements,
      const std::vector<Molecule<double_type>>& molecules);*/
    Eigen::MatrixXdt<double_type> assembleJacobian(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_denities,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<Element<double_type>*>& elements,
      const std::vector<Molecule<double_type>>& molecules);
    Eigen::VectorXdt<double_type> assembleRightHandSide(
      const std::vector<Condensate<double_type>*>& condensates,
      const std::vector<unsigned int>& condensates_jac,
      const std::vector<unsigned int>& condensates_rem,
      const std::vector<double_type>& activity_corr,
      const std::vector<double_type>& number_denities,
      const std::vector< Element<double_type>* >& elements,
      const std::vector< Molecule<double_type> >& molecules,
      const double_type total_element_density, 
      const double_type log_tau);
    std::vector<double_type> solveSystem(
      Eigen::MatrixXdt<double_type>& jacobian,
      Eigen::VectorXdt<double_type>& rhs);
  private:
    FastChemOptions<double_type>& options;
};



}
#endif
