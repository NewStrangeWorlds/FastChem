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


#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <limits>

#include "options.h"


namespace fastchem {


void FastChemOptions::init()
{
  element_density_minlimit = 1e-155;
  molecule_density_minlimit = 1e-155;
}



ParameterFloat FastChemOptions::resolveParameter(
  const std::string& parameter)
{
  static const std::map<std::string, ParameterFloat> parameter_strings = 
  {
    { std::string("condTau"), ParameterFloat::cond_tau },
    { std::string("condIterChangeLimit"), ParameterFloat::cond_limit_change },
    { std::string("accuracyChem"), ParameterFloat::chem_accuracy },
    { std::string("accuracyElementConservation"), ParameterFloat::element_conserve_accuracy },
    { std::string("accuracyCond"), ParameterFloat::cond_accuracy },
    { std::string("accuracyNewton"), ParameterFloat::newton_err },
    { std::string("minDensityExponentElement"), ParameterFloat::element_minlimit },
    { std::string("minDensityExponentMolecules"), ParameterFloat::molecule_minlimit },
    { std::string("condTraceDensityThreshold"), ParameterFloat::cond_trace_density_threshold }
  };

  auto itr = parameter_strings.find(parameter);
  
  if( itr != parameter_strings.end()) 
    return itr->second;

  return ParameterFloat::invalid_parameter; 
}



ParameterBool FastChemOptions::resolveParameterBool(
  const std::string& parameter)
{
  static const std::map<std::string, ParameterBool> parameter_strings = 
  {
    { std::string("condReduceSystemSize"), ParameterBool::cond_reduce_system_size },
    { std::string("condUseSVD"), ParameterBool::cond_use_svd },
    { std::string("useCondDataValidityLimit"), ParameterBool::cond_use_data_validity_limits },
    { std::string("condUseLM"), ParameterBool::cond_use_lm }
  };

  auto itr = parameter_strings.find(parameter);
  
  if( itr != parameter_strings.end()) 
    return itr->second;

  return ParameterBool::invalid_parameter; 
}



ParameterInt FastChemOptions::resolveParameterInt(
  const std::string& parameter)
{
  static const std::map<std::string, ParameterInt> parameter_strings = 
  {
    { std::string("nbIterationsCond"), ParameterInt::nb_max_cond_iter },
    { std::string("nbIterationsBisection"), ParameterInt::nb_max_bisection_iter },
    { std::string("nbIterationsChemCond"), ParameterInt::nb_max_comb_iter },
    { std::string("nbIterationsChem"), ParameterInt::nb_max_fastchem_iter },
    { std::string("nbIterationsNelderMead"), ParameterInt::nb_max_neldermead_iter },
    { std::string("nbIterationsNewton"), ParameterInt::nb_max_newton_iter },
    { std::string("nbSwitchToNewton"), ParameterInt::nb_switch_to_newton },
    { std::string("nbSwitchToJoint"), ParameterInt::nb_switch_to_joint }
  };

  auto itr = parameter_strings.find(parameter);
  
  if( itr != parameter_strings.end()) 
    return itr->second;

  return ParameterInt::invalid_parameter; 
}


}
