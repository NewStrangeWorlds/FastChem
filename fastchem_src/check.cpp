/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
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


#include <vector>
#include <cmath>


namespace fastchem {


//Check for the number density of elements
template <class double_type>
void FastChem<double_type>::checkN(Element<double_type>& species, const double_type h_density, const unsigned int grid_index)
{
  double_type minlimit = element_density_minlimit;

  if (species.number_density[grid_index] < minlimit) species.number_density[grid_index] = minlimit;


  //an element can not be more abundant than determined by its elemental abundance
  double_type maxlimit = h_density*species.abundance;


  if (species.number_density[grid_index] > maxlimit) species.number_density[grid_index] = maxlimit;
}



//Check for the number density of molecules
template <class double_type>
void FastChem<double_type>::checkN(Molecule<double_type>& species, const double_type h_density, const unsigned int grid_index)
{
  double_type minlimit = molecule_density_minlimit;

  if (species.number_density[grid_index] < minlimit) species.number_density[grid_index] = minlimit;


  //a molecule can not be more abundant than its least abundant element
  double_type maxlimit = h_density*species.abundance_scaled;


  if (species.number_density[grid_index] > maxlimit) species.number_density[grid_index] = maxlimit;
}



//Check for charge conservation
template <class double_type>
bool FastChem<double_type>::checkChargeConservation(const unsigned int grid_index)
{
  bool charge_conserved = false;


  unsigned int e_ = getElementIndex("e-");

  if (e_ == FASTCHEM_UNKNOWN_SPECIES) return true; //no electrons in the system


  //if no ions present, charge conservation is automatically satisfied
  if (elements[e_].molecule_list.size() == 0)
  {
    elements[e_].element_conserved[grid_index] = 1;

    return true;
  }


  //sum up all positive and negative charges in the network
  double_type positive_charge = 0;
  double_type negative_charge = elements[e_].number_density[grid_index];


  for (size_t i=0; i<elements[e_].molecule_list.size(); ++i)
  {
    unsigned int molecule_index = elements[e_].molecule_list[i];

    if (molecules[molecule_index].stoichometric_vector[e_] < 0)
      positive_charge -= molecules[molecule_index].stoichometric_vector[e_] * molecules[molecule_index].number_density[grid_index];

    if (molecules[molecule_index].stoichometric_vector[e_] > 0)
      negative_charge += molecules[molecule_index].stoichometric_vector[e_] * molecules[molecule_index].number_density[grid_index];
  }



  if (verbose_level >= 4)
    std::cout << "charge conservation " << positive_charge << "\t"
                                        << negative_charge << "\t"
                                        << std::fabs(positive_charge - negative_charge)/std::sqrt(positive_charge*negative_charge) << "\n";


  if (std::fabs(positive_charge - negative_charge)/std::sqrt(positive_charge*negative_charge) < accuracy)
    charge_conserved = true;
  else
    charge_conserved = false;

  elements[e_].element_conserved[grid_index] = charge_conserved;


  return charge_conserved;
}



//Check for element conservation
template <class double_type>
bool FastChem<double_type>::checkElementConservation(Element<double_type>& species, const double_type h_density, const unsigned int grid_index)
{

  //electrons are subject to charge conservation
  if (species.symbol == "e-")
    return checkChargeConservation(grid_index);


  //sum up the elements contained in each molecule and compare the result to its elemental abundance
  double_type sum = species.number_density[grid_index];


  for (size_t i=0; i<species.molecule_list.size(); ++i)
    sum += molecules[species.molecule_list[i]].stoichometric_vector[species.index] * molecules[species.molecule_list[i]].number_density[grid_index];

  sum /= chemical_elements[species.element_index].abundance * h_density;


  if (verbose_level >= 4)
    std::cout << "element conservation " << species.symbol << "\t" << std::fabs(sum - 1.0L) << "\t"
              << sum*chemical_elements[species.element_index].abundance * h_density << "\t" << chemical_elements[species.element_index].abundance * h_density << "\n";


  if (std::fabs(sum - 1.0L) < accuracy || species.molecule_list.size() == 0)
    species.element_conserved[grid_index] = 1;
  else
    species.element_conserved[grid_index] = 0;


  return species.element_conserved[grid_index];
}



template class FastChem<double>;
template class FastChem<long double>;

}
