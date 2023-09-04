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

#ifndef _elements_h
#define _elements_h

#include <string>
#include <vector>

#include "../species_struct.h"

#include "../fastchem_constants.h"


namespace fastchem {


template <class double_type>
class ElementData {
  public:
    ElementData(const std::string& abundance_file, const std::string& chemical_element_file);
    ElementData(const ElementData &obj);

    std::vector<Element<double_type>> elements;
    std::vector<Element<double_type>*> elements_wo_e;

    size_t nb_elements = 0;
    unsigned int e_ = FASTCHEM_UNKNOWN_SPECIES; //electron element index

    bool is_initialised = false;

    unsigned int elementIndex(const std::string& symbol);
    void setAbundances(const std::vector<double>& abundances);
    void setRelativeAbundances();
    void init(double_type initial_density);
  private:
    std::vector<ChemicalElementData> chemical_element_data;
    size_t nb_chemical_element_data = 0;

    unsigned int chemicalElementIndex(const std::string& symbol);

    void setAbundance(const std::string& symbol, const double abundance);
    void add(const std::string& symbol);
    bool readElementList(const std::string& file_path);
    bool readElementAbundances(const std::string& file_path);
};



}

#endif
