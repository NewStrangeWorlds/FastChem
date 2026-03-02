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


#ifndef _gas_phase_solver_h
#define _gas_phase_solver_h

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
class GasPhaseSolver{
  public:
    GasPhaseSolver(FastChemOptions& options_)
      : options(options_)
      {}

    void intertSol(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double gas_density);
    void linSol(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double gas_density);
    void quadSol(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double gas_density);
    void newtonSol(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double gas_density,
      const bool use_alternative);

    void selectNewtonElements(
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules,
      const std::vector<double>& old_number_densities,
      const double gas_density,
      std::vector<Element*>& newton_elements);

    void newtonSolMult(
      std::vector<Element*>& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double gas_density);

    void newtonSolElectron(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules,
      const double gas_density);

    bool nelderMeadElectron(
      Element& electron,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double initial_solution,
      const double gas_density);


    bool bisection(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules,
      const double gas_density,
      const bool use_all_molecules);

    void backupSol(
      Element& species,
      std::vector<Element>& elements,
      const std::vector<Molecule>& molecules, 
      const double gas_density);

    int order_cation = -999;
    int order_anion = -999;
   
  private:
    FastChemOptions& options;
    
    void logSpaceResidual(
      const Element& species,
      const std::vector<Element>& elements,
      const std::vector<Molecule>& molecules,
      const double gas_density,
      const double y_j,
      double& ln_P,
      double& ln_dP,
      double& R,
      const bool use_all_molecules);

    Eigen::VectorXdt assembleJacobian(
      const std::vector<Element*>& species,
      const std::vector< Element >& elements,
      const std::vector< Molecule >& molecules,
      Eigen::MatrixXdt& jacobian);

    void assembleRightHandSide(
      const std::vector<Element*>& species,
      const std::vector< Element >& elements,
      const std::vector< Molecule >& molecules,
      const double gas_density,
      const Eigen::VectorXdt& scaling_factors,
      Eigen::VectorXdt& rhs);

    double AmCoeffElectron(
      const Element& electron,
      const std::vector<Element>& elements,
      const std::vector<Molecule>& molecules,
      const int order);
};



}
#endif
