
from dataclasses import dataclass

def readElementAbundancesFile(file_path):
  element_abundances = {}
  
  with open(file_path, 'r') as file:
    for line in file:
      if line.startswith('#') or not line.strip():
        continue
      
      parts = line.split()
      element = parts[0]
      abundance = float(parts[1])
      element_abundances[element] = abundance
  
  return element_abundances


def selectElementAbundances(element_abundances, element_list):
  selected_elements = {}
  
  for element in element_list:
    if element in element_abundances:
      selected_elements[element] = element_abundances[element]
    else:
      raise ValueError(f"Element {element} not found in abundance data.")
  
  return selected_elements

@dataclass
class ChemistryData:
    elements : list
    data : list


def readGasSpeciesFile(file_path):
  species_data = []
    
  with open(file_path, 'r') as file:
    for line in file:
      if line.startswith('#') or not line.strip():
        continue
      
      header = line
      data = file.readline()

      element_parts = line.split(":")[1].split("#")
      
      elements = element_parts[0].strip().split()
      elements = elements[::2]
      
      species_data.append(ChemistryData(elements=elements, data=[header, data]))

  return species_data

def readCondSpeciesFile(file_path):
  species_data = []
    
  with open(file_path, 'r') as file:
    for line in file:
      if line.startswith('#') or not line.strip():
        continue
      
      header = line
      element_parts = line.split(":")[1].split("#")
      
      elements = element_parts[0].strip().split()
      elements = elements[::2]
      
      species_data.append(ChemistryData(elements=elements, data=[header]))
      
      for data_line in file:
        if data_line.startswith('#') or not data_line.strip():
          break
        else:
          species_data[-1].data.append(data_line)

  return species_data


def selectSpeciesByElements(species_data, element_list):
  selected_species = []
  
  for species in species_data:
    if all(elem in element_list for elem in species.elements):
      selected_species.append(species)
  
  return selected_species


def writeElementAbundancesFile(file_path, element_abundances):
  with open(file_path, 'w') as file:
    file.write("# Custom Element Abundances\n")
    
    for element, abundance in element_abundances.items():
      file.write(f"{element} {abundance}\n")
  
  return None

def writeSpeciesFile(file_path, species_data):
  with open(file_path, 'w') as file:
    file.write("#logK = a1/T + a2 ln T + a3 + a4 T + a5 T^2 for FastChem:\n")
    file.write("#Includes customly selected elements.\n")
    file.write("#fit coefficients calculated from indicated data source.\n")
    
    for species in species_data:
      for line in species.data:
        file.write(line)
      
      file.write("\n")
  
  return None