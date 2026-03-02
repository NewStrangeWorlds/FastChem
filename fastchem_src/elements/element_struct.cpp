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


#include "../species_struct.h"

#include <cmath>
#include <iostream>


namespace fastchem {


//Determination of the correction factor n_j_min, that contains the number density of an element j contained in molecules
//with elements that are less abundant than j
//See Eq. (2.23)
void Element::calcMinorSpeciesDensities(
  const std::vector< Molecule > &molecules)
{
  number_density_min = 0.0;

  for (auto & i : minor_molecules)
    number_density_min += (molecules[i].stoichiometric_vector[index] 
                       + epsilon * molecules[i].sigma) * molecules[i].number_density;
}



void Element::calcEpsilon(const std::vector< Element > &elements)
{ 
  double element_sum = 0.0;

  for (auto & i : elements)
    element_sum += i.abundance;

  epsilon = this->abundance/element_sum;
}



//Check for the number density of elements
void Element::checkN(
  const double& min_limit, const double& gas_density)
{
  if (!this->fixed_by_condensation)
  {
    if (this->log_number_density < static_cast<double>(LOG_DENSITY_FLOOR))
      this->log_number_density = static_cast<double>(LOG_DENSITY_FLOOR);
  }

  const double log_gas = std::log(gas_density);

  if (this->log_number_density > log_gas)
    this->log_number_density = log_gas;

  this->number_density = safeExp(this->log_number_density);
}


//Check for charge conservation
bool Element::checkChargeConservation(
  const std::vector< Molecule >& molecules, const double& accuracy)
{
  //Am I the electron?
  if (this->symbol != "e-") return false;

  //If no ions present, charge conservation is automatically satisfied
  if (molecule_list.size() == 0)
  {
    element_conserved = 1;

    return true;
  }

  bool charge_conserved = false;

  //sum up all positive and negative charges in the network
  double positive_charge = 0;
  double negative_charge = this->number_density;

  for (auto & i : molecule_list)
  {
    if (molecules[i].stoichiometric_vector[index] < 0)
      positive_charge -= molecules[i].stoichiometric_vector[index] * molecules[i].number_density;

    if (molecules[i].stoichiometric_vector[index] > 0)
      negative_charge += molecules[i].stoichiometric_vector[index] * molecules[i].number_density;
  }
  
  if (std::fabs(positive_charge - negative_charge)/std::sqrt(positive_charge*negative_charge) < accuracy 
    || (positive_charge == 0 && negative_charge == 0))
    charge_conserved = true;
  else
    charge_conserved = false;

  element_conserved = charge_conserved;

  return charge_conserved;
}



//Check for element conservation
bool Element::checkElementConservation(
  const std::vector< Molecule >& molecules, 
  const std::vector< Condensate >& condensates,
  const double total_density,
  const double& accuracy)
{
  //electrons are subject to charge conservation
  if (this->symbol == "e-")
    return checkChargeConservation(molecules, accuracy);

  //if the element abundance is 0, it is automatically conserved
  if (this->epsilon == 0)
    return true;

  //sum up the elements contained in each molecule and compare the result to its elemental abundance
  double sum_gas = this->number_density;


  for (auto & i : molecule_list)
    sum_gas += molecules[i].stoichiometric_vector[index] * molecules[i].number_density;

  double sum_cond = 0;

  for (auto & i : condensate_list)
    sum_cond += condensates[i].stoichiometric_vector[index] * condensates[i].number_density;
  
  double sum_total = sum_gas + sum_cond;

  sum_total /= total_density*epsilon;
  
  if (std::fabs(sum_total - 1.0L) < accuracy || molecule_list.size() == 0)
    element_conserved = 1;
  else
    element_conserved = 0;

  return element_conserved;
}



void Element::calcDegreeOfCondensation(
  const std::vector< Condensate > &condensates,
  const double total_element_density)
{
  if (this->symbol == "e-") return;

  double density_cond = 0;

  for (auto & i : condensate_list)
  {
    density_cond += condensates[i].stoichiometric_vector[this->index] * condensates[i].number_density;
  }
    

  degree_of_condensation = density_cond/(epsilon*total_element_density);

  if (degree_of_condensation > 1) 
    degree_of_condensation = 1.0;

  if (this->epsilon == 0) 
    degree_of_condensation = 0;

  phi = this->epsilon * (1.0 - degree_of_condensation);
}


void Element::normalisePhi(const double element_phi_sum)
{

  phi /= element_phi_sum;

}



}
