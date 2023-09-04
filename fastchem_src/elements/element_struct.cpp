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


#include "../species_struct.h"

#include <cmath>
#include <iostream>


namespace fastchem {


//Determination of the correction factor n_j_min, that contains the number density of an element j contained in molecules
//with elements that are less abundant than j
//See Eq. (2.23)
template <class double_type>
void Element<double_type>::calcMinorSpeciesDensities(
  const std::vector< Molecule<double_type> > &molecules)
{
  number_density_min = 0.0;

  for (auto & i : minor_molecules)
    number_density_min += (molecules[i].stoichiometric_vector[index] 
                       + epsilon * molecules[i].sigma) * molecules[i].number_density;
}



template <class double_type>
void Element<double_type>::calcEpsilon(const std::vector< Element<double_type> > &elements)
{ 
  double_type element_sum = 0.0;

  for (auto & i : elements)
    element_sum += i.abundance;

  epsilon = this->abundance/element_sum;
}



template <class double_type>
void Element<double_type>::calcSolverScalingFactor(
  const std::vector<Element<double_type>> &elements, 
  const std::vector<Molecule<double_type>> &molecules,
  const double additional_scaling_factor)
{
  solver_scaling_factor = 0.0;


  for (auto & i : molecule_list)
  {
    if (molecules[i].stoichiometric_vector[index] < 1 || molecules[i].stoichiometric_vector[index] > static_cast<int>(solver_order) )
      continue;

    double_type sum = 0.0;

    if (molecules[i].abundance == this->abundance)
    {
      for (auto & l : molecules[i].element_indices)
      {
        if (l != index)
          sum += molecules[i].stoichiometric_vector[l] * std::log(elements[l].number_density);

      }

      sum += molecules[i].mass_action_constant;
    }


    if (sum > solver_scaling_factor)
      solver_scaling_factor = sum;
  }

  //scale the factor by an additional, optional factor supplied by the user
  solver_scaling_factor -= additional_scaling_factor;
}


//Check for the number density of elements
template <class double_type>
void Element<double_type>::checkN(
  const double_type& min_limit, const double_type& gas_density)
{
  if (this->number_density < min_limit) this->number_density = min_limit;

  if (this->number_density > gas_density) this->number_density = gas_density;
}


//Check for charge conservation
template <class double_type>
bool Element<double_type>::checkChargeConservation(
  const std::vector< Molecule<double_type> >& molecules, const double_type& accuracy)
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
  double_type positive_charge = 0;
  double_type negative_charge = this->number_density;

  for (auto & i : molecule_list)
  {
    if (molecules[i].stoichiometric_vector[index] < 0)
      positive_charge -= molecules[i].stoichiometric_vector[index] * molecules[i].number_density;

    if (molecules[i].stoichiometric_vector[index] > 0)
      negative_charge += molecules[i].stoichiometric_vector[index] * molecules[i].number_density;
  }


  if (std::fabs(positive_charge - negative_charge)/std::sqrt(positive_charge*negative_charge) < accuracy)
    charge_conserved = true;
  else
    charge_conserved = false;

  element_conserved = charge_conserved;


  return charge_conserved;
}



//Check for element conservation
template <class double_type>
bool Element<double_type>::checkElementConservation(
  const std::vector< Molecule<double_type> >& molecules, 
  const std::vector< Condensate<double_type> >& condensates,
  const double_type total_density,
  const double_type& accuracy)
{
  //electrons are subject to charge conservation
  if (this->symbol == "e-")
    return checkChargeConservation(molecules, accuracy);

  //if the element abundance is 0, it is automatically conserved
  if (this->epsilon == 0)
    return true;

  //sum up the elements contained in each molecule and compare the result to its elemental abundance
  double_type sum_gas = this->number_density;


  for (auto & i : molecule_list)
    sum_gas += molecules[i].stoichiometric_vector[index] * molecules[i].number_density;

  double_type sum_cond = 0;

  for (auto & i : condensate_list)
    sum_cond += condensates[i].stoichiometric_vector[index] * condensates[i].number_density;
  
  double_type sum_total = sum_gas + sum_cond;

  sum_total /= total_density*epsilon;
  
  //std::cout << this->symbol << "\t" << sum_gas << "\t" << sum_cond << "\t" << sum_total << "\t" << total_density*epsilon << "\t" << phi << "\t" << epsilon << "\n";

  if (std::fabs(sum_total - 1.0L) < accuracy || molecule_list.size() == 0)
    element_conserved = 1;
  else
    element_conserved = 0;

  //if (!element_conserved) exit(0);
  return element_conserved;
}



template <class double_type>
void Element<double_type>::calcDegreeOfCondensation(
  const std::vector< Condensate<double_type> > &condensates,
  const double_type total_element_density)
{
  if (this->symbol == "e-") return;

  double_type density_cond = 0;

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


template <class double_type>
void Element<double_type>::normalisePhi(const double_type element_phi_sum)
{

  phi /= element_phi_sum;

}



template struct Element<double>;
template struct Element<long double>;
}
