 
import pyfastchem
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

#input values for temperature and pressure (only a single values here)
temperature = np.full(1, 1000)
pressure = np.full(1, 1)

#the range of C/O ratios we want to calculate the chemistry for
c_to_o = np.linspace(0.1, 10, 100)


#the chemical species we want to plot later
#note that FastChem uses the Hill notation
plot_species = ['H2O1', 'C1O2', 'C1O1', 'C1H4', 'H3N1', 'C2H2']
#for the plot lables, we therefore use separate strings in the usual notation
plot_species_lables = ['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'C2H2']


#create a FastChem object
#fastchem = pyfastchem.FastChem('input/parameters.dat', 1)
fastchem = pyfastchem.FastChem('input/element_abundances_jo_orich.dat', 'input/logK_full_new.dat', 1)


#allocate the matrix for the number densities
number_densities = np.zeros((c_to_o.size, fastchem.getSpeciesNumber()))


#make a copy of the solar abundances from FastChem
solar_abundances = np.array(fastchem.getElementAbundances())


#create the input and output structures for FastChem
input_data = pyfastchem.FastChemInput()
output_data = pyfastchem.FastChemOutput()

input_data.temperature = temperature
input_data.pressure = pressure


#we need to know the indices for O and C from FastChem
index_C = fastchem.getSpeciesIndex('C')
index_O = fastchem.getSpeciesIndex('O')

#loop over the C/O ratios
for i in range(0, c_to_o.size):
  element_abundances = np.copy(solar_abundances)
  
  #set the C abundance as a function of the C/O ratio
  element_abundances[index_C] = element_abundances[index_O] * c_to_o[i]

  fastchem.setElementAbundances(element_abundances)

  fastchem_flag = fastchem.calcDensities(input_data, output_data)
  print("FastChem reports:", pyfastchem.FASTCHEM_MSG[fastchem_flag])
  
  #copy the FastChem output into the number density matrix
  number_densities[i,:] = np.array(output_data.number_densities[0])



#total gas particle number density from the ideal gas law 
gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


#check the species we want to plot and get their indices from FastChem
plot_species_indices = []
plot_species_symbols = []

for i, species in enumerate(plot_species):
  index = fastchem.getSpeciesIndex(species)

  if index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
    plot_species_indices.append(index)
    plot_species_symbols.append(plot_species_lables[i])
  else:
    print("Species", species, "to plot not found in FastChem")


#and plot...
for i in range(0, len(plot_species_symbols)):
  plt.plot(c_to_o, number_densities[:, plot_species_indices[i]]/gas_number_density)

plt.yscale('log')

plt.ylabel("Mixing ratios")
plt.xlabel("C/O")
plt.legend(plot_species_symbols)

plt.show()
