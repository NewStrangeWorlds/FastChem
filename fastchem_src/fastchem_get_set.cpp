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


#include "fastchem.h"
#include "options.h"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <stack>
#include <cctype>
#include <sstream>
#include <regex>



namespace fastchem {


//Query for a species index (both, elements and molecules) with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getGasSpeciesIndex(const std::string symbol)
{
  auto it = std::find_if(
    gas_phase.species.begin(),
    gas_phase.species.end(),
    [&] (const ChemicalSpecies<double_type>* a) {
      return a->symbol == symbol;});

  if (it == gas_phase.species.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - gas_phase.species.begin();
}


//Query for a condensate species index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getElementIndex(const std::string symbol)
{
  auto it = std::find_if(
    element_data.elements.begin(),
    element_data.elements.end(),
    [&] (const Element<double_type> a) {
      return a.symbol == symbol;});

  if (it == element_data.elements.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - element_data.elements.begin();
}


//Query for a condensate species index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getCondSpeciesIndex(const std::string symbol)
{
  auto it = std::find_if(
    condensed_phase.condensates.begin(),
    condensed_phase.condensates.end(),
    [&] (const Condensate<double_type> a) {
      return a.symbol == symbol;});

  if (it == condensed_phase.condensates.end()) 
    return FASTCHEM_UNKNOWN_SPECIES;
  else
    return it - condensed_phase.condensates.begin();
}



//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getGasSpeciesName(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->name;
  else
    return "";
}


//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getElementName(const unsigned int species_index)
{
  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].name;
  else
    return "";
}


//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getCondSpeciesName(const unsigned int species_index)
{
  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].name;
  else
    return "";
}




//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getGasSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->symbol;
  else
    return "";
}



//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getElementSymbol(const unsigned int species_index)
{

  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].symbol;
  else
    return "";
}



//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getCondSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].symbol;
  else
    return "";
}



//Query for a the molecular weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getGasSpeciesWeight(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_species)
    return gas_phase.species[species_index]->weight;
  else
    return 0.;
}


//Query for a the weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getElementWeight(const unsigned int species_index)
{
  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].weight;
  else
    return 0.;
}


//Query for a the weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getCondSpeciesWeight(const unsigned int species_index)
{
  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].weight;
  else
    return 0.;
}


//Get the element abundance for a specific element
template <class double_type>
double FastChem<double_type>::getElementAbundance(const unsigned int species_index)
{
  if (species_index < element_data.nb_elements)
    return element_data.elements[species_index].abundance;
  else
    return 0.;
}



//Get the element abundandes for all elements
template <class double_type>
std::vector<double> FastChem<double_type>::getElementAbundances()
{
  std::vector<double> abundances(element_data.nb_elements, 0.0);

  for (size_t i=0; i<element_data.nb_elements; ++i)
    abundances[i] = element_data.elements[i].abundance;

  return abundances;
}


//Set the element abundances for all elements
template <class double_type>
void FastChem<double_type>::setElementAbundances(std::vector<double> abundances)
{
  element_data.setAbundances(abundances);

  gas_phase.reInitialise();
}



template <class double_type>
std::vector<int> FastChem<double_type>::getGasSpeciesStoichiometry(const unsigned int species_index)
{
  if (species_index < gas_phase.nb_elements)
  {
    std::vector<int> stoichiometry(gas_phase.nb_elements, 0);
    stoichiometry[species_index] = 1;
    
    return stoichiometry;
  }

  if (species_index < gas_phase.nb_species)
    return gas_phase.molecules[species_index - gas_phase.nb_elements].stoichiometric_vector;
  else
    return std::vector<int>();
}



template <class double_type>
std::vector<int> FastChem<double_type>::getCondSpeciesStoichiometry(const unsigned int species_index)
{
  if (species_index < condensed_phase.nb_condensates)
    return condensed_phase.condensates[species_index].stoichiometric_vector;
  else
    return std::vector<int>();
}


//recursive parser
std::map<std::string, int> parseFormula(const std::string& formula, size_t& i) 
{
  auto mergeCounts = [](
    std::map<std::string, 
    int>& target, 
    const std::map<std::string, int>& source, 
    int multiplier) 
    {
      for (const auto& pair : source)
        target[pair.first] += pair.second * multiplier;
    };

  std::map<std::string, int> element_counts;

  while (i < formula.size()) 
  {
    if (std::isupper(formula[i])) 
    {
      //parse element symbol
      std::string element;

      element += formula[i++];
      
      if (i < formula.size() && std::islower(formula[i]))
        element += formula[i++];

      //parse count
      std::string count_str;

      while (i < formula.size() && std::isdigit(formula[i])) 
        count_str += formula[i++];
      
      int count = count_str.empty() ? 1 : std::stoi(count_str);
      element_counts[element] += count;
      } 
      else if (formula[i] == '(')
      {
        ++i; //skip '('
        std::map<std::string, int> innerCounts = parseFormula(formula, i); // Parse inner
        
        if (i >= formula.size() || formula[i] != ')') 
        {
          std::cout << "Mismatched parentheses in chemical formula\n";
          return {};
        }
        
        ++i; //skip ')'

        //parse multiplier after parenthesis
        std::string countStr;
        while (i < formula.size() && std::isdigit(formula[i]))
          countStr += formula[i++];

        int multiplier = countStr.empty() ? 1 : std::stoi(countStr);
        mergeCounts(element_counts, innerCounts, multiplier);
      } 
      else if (formula[i] == ')') 
      {
        break; //let caller handle it
      } 
      else 
      {
        std::cout << "Invalid character in formula: " << formula[i] << "\n";
        return {};
      }
    }

    return element_counts;
}


std::pair<std::map<std::string, int>, std::string> parseFormulaWithCharge(
  const std::string& rawFormula) 
{
  std::string formula = rawFormula;
  std::string charge;

  //Extract trailing charge using regex
  std::smatch match;
  std::regex charge_pattern(R"(^(.*?)(\^{1}\d*[+-]{1}|[+-]{1,2})$)");
  
  if (std::regex_match(formula, match, charge_pattern)) 
  {
    std::string body = match[1];
    std::string charge_part = match[2];

    //Make sure it's actually a charge
    if (!charge_part.empty() && (charge_part[0] == '^' || charge_part.back() == '+' || charge_part.back() == '-')) 
    {
      formula = body;
      charge = charge_part;
      
      //Treat special cases
      if (charge == "^+" || charge == "^1+")
        charge = "+";
      
      if (charge == "^-" || charge == "^1-")
        charge = "-";

      if (charge == "^2+")
        charge = "++";

      if (charge == "^2-")
        charge = "--";
    }
  }
  
  size_t index = 0;
  std::map<std::string, int> counts = parseFormula(formula, index);
  
  return { counts, charge };
}



template <class double_type>
std::string FastChem<double_type>::convertToHillNotation(const std::string& formula) const
{
  if (formula.empty())
    return "";
  
  //first, we treat some special isomer cases
  if (formula == "AlOH")
    return "Al1H1O1_1";

  if (formula == "OAlH")
    return "Al1H1O1_2";

  if (formula == "HCN")
    return "C1H1N1_1";

  if (formula == "HNC")
    return "C1H1N1_2";
  
  if (formula == "FS2F")
    return "F2S2_1";

  if (formula == "SSF2" || formula == "S2F2")
    return "F2S2_2";


  auto parsed_formula = parseFormulaWithCharge(formula);

  auto element_counts = parsed_formula.first;
  auto charge = parsed_formula.second;
  
  if (element_counts.empty())
    return "";

  //sort elements in Hill notation
  //C and H first, then alphabetical order for others
  //use a vector to store pairs of element symbol and count
  std::vector<std::pair<std::string, int>> elements;

  //check if it contains carbon
  bool has_carbon = element_counts.count("C") > 0;

  if (has_carbon) 
  {
    //add C and H first if present
    if (element_counts.count("C")) 
    {
      elements.emplace_back("C", element_counts["C"]);
      element_counts.erase("C");
    }

    if (element_counts.count("H")) 
    {
      elements.emplace_back("H", element_counts["H"]);
      element_counts.erase("H");
    }
  }

  //add remaining elements alphabetically
  for (const auto& e : element_counts) 
    elements.emplace_back(e.first, e.second);
  
  if (elements.size() > 1)
    std::sort(elements.begin() + (has_carbon ? 2 : 0), elements.end());

  //if we only have an element, return its symbol
  if (elements.size() == 1 && elements[0].second == 1 && charge.empty())
    return elements[0].first; 


  //build final formula string
  std::string result;

  for (const auto& e : elements) 
  {
    //std::cout << e.first << "  " << e.second << std::endl;
    result += e.first;
    result += std::to_string(e.second);
  }
  //std::cout << charge << "\n";
  if (!charge.empty())
    result += charge;

  return result;
}




//Set an internal FastChem parameter (for double values)
template <class double_type>
void FastChem<double_type>::setParameter(const std::string& parameter, const double_type value)
{
  auto param = options.resolveParameter(parameter);

  switch (param)
  {
    case ParameterFloat::cond_tau:
      options.cond_tau = value;
      break;
    
    case ParameterFloat::cond_limit_change:
      options.cond_iter_change_limit = value;
      break;

    case ParameterFloat::chem_accuracy:
      options.chem_accuracy = value;
      break;

    case ParameterFloat::element_conserve_accuracy:
      options.element_conserve_accuracy = value;
      break;

    case ParameterFloat::cond_accuracy:
      options.cond_accuracy = value;
      break;

    case ParameterFloat::newton_err:
      options.newton_err = value;
      break;

    case ParameterFloat::element_minlimit:
      options.element_density_minlimit = std::pow(10.0, value);
      break;

    case ParameterFloat::molecule_minlimit:
      options.molecule_density_minlimit = std::pow(10.0, value);
      break;

    case ParameterFloat::logK_limit:
      options.logK_limit = value;
      break;

    case ParameterFloat::additional_scaling_factor:
      options.additional_scaling_factor = value;
      break;
  
    default:
      std::cout << "Unknown parameter \"" << parameter << "\"  with a floatint-point value!\n";
      break;
  }

}


//Set an internal FastChem parameter (for boolean values)
template <class double_type>
void FastChem<double_type>::setParameter(const std::string& parameter, const bool value)
{
  auto param = options.resolveParameterBool(parameter);

  switch (param)
  {
    case ParameterBool::cond_solve_full_system:
      options.cond_solve_full_matrix = value;
      break;

    case ParameterBool::cond_reduce_system_size:
      options.cond_reduce_system_size = value;
      break;

    case ParameterBool::cond_use_full_pivot:
      options.cond_use_full_pivot = value;
      break;

    case ParameterBool::cond_use_svd:
      options.cond_use_svd = value;
      break;

    case ParameterBool::use_scaling_factor:
      options.use_scaling_factor = value;
      break;

     case ParameterBool::cond_use_data_validity_limits:
      options.cond_use_data_validity_limits = value;
      break;
  
    default:
      std::cout << "Unknown parameter \"" << parameter << "\"  with a boolean value!\n";
      break;
  }
}



//Set an internal FastChem parameter (for integer values)
template <class double_type>
void FastChem<double_type>::setParameter(const std::string& parameter, const unsigned int value)
{
  auto param = options.resolveParameterInt(parameter);

  switch (param)
  {
    case ParameterInt::nb_max_bisection_iter:
      options.nb_max_bisection_iter = value;
      break;

    case ParameterInt::nb_max_comb_iter:
      options.nb_chem_cond_iter = value;
      break;

    case ParameterInt::nb_max_cond_iter:
      options.nb_max_cond_iter = value;
      break;

    case ParameterInt::nb_max_fastchem_iter:
      options.nb_max_fastchem_iter = value;
      break;

    case ParameterInt::nb_max_neldermead_iter:
      options.nb_max_neldermead_iter = value;
      break;

    case ParameterInt::nb_max_newton_iter:
      options.nb_max_newton_iter = value;
      break;

    case ParameterInt::nb_switch_to_newton:
      options.nb_switch_to_newton = value;
      break;
  
    default:
      std::cout << "Unknown parameter \"" << parameter << "\"  with an integer value!\n";
      break;
  }
}



template class FastChem<double>;
template class FastChem<long double>;
}
