 
import pyfastchem
import numpy as np
import os
from save_output import saveChemistryOutput, saveMonitorOutput, saveChemistryOutputPandas, saveMonitorOutputPandas
import matplotlib.pyplot as plt
from astropy import constants as const

#input values for temperature (in K) and pressure (in bar)
#we only use a single values here
temperature_single = 1000
pressure_single = 1

#the range of solar metallicities we want to calculate the chemistry for
metallicity = np.logspace(-1, 3, 100)


#define the directory for the output
#here, we currently use the standard one from FastChem
output_dir = '../output'


#the chemical species we want to plot later
#note that the standard FastChem input files use the Hill notation
plot_species = ['H2', 'H2O1', 'C1O2', 'C1O1', 'C1H4', 'H3N1']
#for the plot labels, we therefore use separate strings in the usual notation
plot_species_labels = ['H2', 'H2O', 'CO2', 'CO', 'CH4', 'NH3']


#create a FastChem object
#it needs the locations of the element abundance and equilibrium constants files
#these locations have to be relative to the one this Python script is called from
fastchem = pyfastchem.FastChem(
  '../input/element_abundances/asplund_2009.dat',
  '../input/logK/logK.dat',
  1)


#we could also create a FastChem object by using the parameter file
#note, however, that the file locations in the parameter file are relative
#to the location from where this Python script is called from
#fastchem = pyfastchem.FastChem('../input/parameters_py.dat', 1)


#make a copy of the solar abundances from FastChem
solar_abundances = np.array(fastchem.getElementAbundances())


#allocate the data for the output
nb_points = metallicity.size

number_densities = np.zeros((nb_points, fastchem.getGasSpeciesNumber()))
total_element_density = np.zeros(nb_points)
mean_molecular_weight = np.zeros(nb_points)
element_conserved = np.zeros((nb_points, fastchem.getElementNumber()), dtype=int)
fastchem_flags = np.zeros(nb_points, dtype=int)
nb_iterations = np.zeros(nb_points, dtype=int)
nb_chemistry_iterations = np.zeros(nb_points, dtype=int)
nb_cond_iterations = np.zeros(nb_points, dtype=int)

temperature = np.zeros(nb_points)
pressure = np.zeros(nb_points)


#loop over the metallicity values
#since we constantly change the element abundances, we have to call FastChem separately each time
for i in range(0, nb_points):
  element_abundances = np.copy(solar_abundances)
  
  #scale the element abundances, except those of H and He
  for j in range(0, fastchem.getElementNumber()):
    if fastchem.getElementSymbol(j) != 'H' and fastchem.getElementSymbol(j) != 'He':
      element_abundances[j] *= metallicity[i]

  fastchem.setElementAbundances(element_abundances)
  
  #create the input and output structures for FastChem
  input_data = pyfastchem.FastChemInput()
  output_data = pyfastchem.FastChemOutput()

  input_data.temperature = [temperature_single]
  input_data.pressure = [pressure_single]

  fastchem_flag = fastchem.calcDensities(input_data, output_data)

  #copy the FastChem input and output into the pre-allocated arrays
  temperature[i] = input_data.temperature[0]
  pressure[i] = input_data.pressure[0]

  number_densities[i,:] = np.array(output_data.number_densities[0])

  total_element_density[i] = output_data.total_element_density[0]
  mean_molecular_weight[i] = output_data.mean_molecular_weight[0]
  element_conserved[i,:] = output_data.element_conserved[0]
  fastchem_flags[i] = output_data.fastchem_flag[0]
  nb_iterations[i] = output_data.nb_iterations[0]
  nb_chemistry_iterations[i] = output_data.nb_chemistry_iterations[0]
  nb_cond_iterations[i] = output_data.nb_cond_iterations[0]


#convergence summary report
print("FastChem reports:")
print("  -", pyfastchem.FASTCHEM_MSG[np.max(fastchem_flag)])

if np.amin(output_data.element_conserved) == 1:
  print("  - element conservation: ok")
else:
  print("  - element conservation: fail")


#total gas particle number density from the ideal gas law 
gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


#check if output directory exists
#create it if it doesn't
os.makedirs(output_dir, exist_ok=True)


#save the monitor output to a file
#we add an additional output column for the metallicity
saveMonitorOutput(output_dir + '/monitor.dat', 
                  temperature, pressure, 
                  element_conserved,
                  fastchem_flags,
                  nb_iterations,
                  nb_chemistry_iterations,
                  nb_cond_iterations,
                  total_element_density,
                  mean_molecular_weight,
                  fastchem,
                  metallicity, 'M/H')


#this saves the output of all species
#we add an additional column for the metallicity
saveChemistryOutput(output_dir + '/chemistry.dat', 
                    temperature, pressure,
                    total_element_density, 
                    mean_molecular_weight, 
                    number_densities,
                    fastchem, 
                    None, 
                    metallicity, 'M/H')


#save the monitor output to a file
#here, the data is saved as a pandas DataFrame inside a pickle file
#we add an additional output column for the metallicity
# saveMonitorOutputPandas(output_dir + '/monitor.pkl', 
#                   temperature, pressure, 
#                   element_conserved,
#                   fastchem_flags,
#                   nb_chemistry_iterations,
#                   total_element_density,
#                   mean_molecular_weight,
#                   fastchem,
#                   metallicity, 'M/H')


#this would save the output of all species
#here, the data is saved as a pandas DataFrame inside a pickle file
# saveChemistryOutputPandas(output_dir + '/chemistry.pkl', 
#                     temperature, pressure,
#                     total_element_density, 
#                     mean_molecular_weight, 
#                     number_densities,
#                     fastchem, 
#                     None, 
#                     metallicity, 'M/H')


#check the species we want to plot and get their indices from FastChem
plot_species_indices = []
plot_species_symbols = []

for i, species in enumerate(plot_species):
  index = fastchem.getGasSpeciesIndex(species)

  if index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
    plot_species_indices.append(index)
    plot_species_symbols.append(plot_species_labels[i])
  else:
    print("Species", species, "to plot not found in FastChem")


#and plot...
for i in range(0, len(plot_species_symbols)):
  plt.plot(metallicity, number_densities[:, plot_species_indices[i]]/gas_number_density)

plt.xscale('log')
plt.yscale('log')

plt.ylabel("Mixing ratios")
plt.xlabel("M/H")
plt.legend(plot_species_symbols)

plt.show()

#we could also save the figure as a pdf
#plt.savefig(output_dir + '/fastchem_metallicity_fig.pdf')
