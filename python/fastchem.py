
import pyfastchem
from save_output import saveChemistryOutput, saveMonitorOutput, saveChemistryOutputPandas, saveMonitorOutputPandas
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import constants as const


#some input values for temperature (in K) and pressure (in bar)
temperature = np.full(1000, 500)
pressure = np.logspace(-6, 1, num=1000)


#define the directory for the output
#here, we currently use the standard one from FastChem
output_dir = '../output'


#the chemical species we want to plot later
#note that the standard FastChem input files use the Hill notation
plot_species = ['H2O1', 'C1O2', 'C1O1', 'C1H4', 'H3N1']
#for the plot labels, we therefore use separate strings in the usual notation
plot_species_labels = ['H2O', 'CO2', 'CO', 'CH4', 'NH3']



#First, we have to create a FastChem object
fastchem = pyfastchem.FastChem(
 '../input/element_abundances/asplund_2009.dat', 
 '../input/logK/logK.dat',
 1)


#we could also create a FastChem object by using the parameter file
#note, however, that the file locations in the parameter file are relative
#to the location from where this Python script is called from
#fastchem = pyfastchem.FastChem('../input/parameters_py.dat', 1)



#create the input and output structures for FastChem
input_data = pyfastchem.FastChemInput()
output_data = pyfastchem.FastChemOutput()

input_data.temperature = temperature
input_data.pressure = pressure


#run FastChem on the entire p-T structure
fastchem_flag = fastchem.calcDensities(input_data, output_data)

print("FastChem reports:")
print("  -", pyfastchem.FASTCHEM_MSG[fastchem_flag])

if np.amin(output_data.element_conserved[:]) == 1:
  print("  - element conservation: ok")
else:
  print("  - element conservation: fail")



#check if output directory exists
#create it if it doesn't
os.makedirs(output_dir, exist_ok=True)


#save the monitor output to a file
saveMonitorOutput(output_dir + '/monitor.dat', 
                  temperature, pressure, 
                  output_data.element_conserved,
                  output_data.fastchem_flag,
                  output_data.nb_iterations,
                  output_data.nb_chemistry_iterations,
                  output_data.nb_cond_iterations,
                  output_data.total_element_density,
                  output_data.mean_molecular_weight,
                  fastchem)

#this would save the output of all species
saveChemistryOutput(output_dir + '/chemistry.dat', 
                    temperature, pressure, 
                    output_data.total_element_density, 
                    output_data.mean_molecular_weight, 
                    output_data.number_densities, 
                    fastchem)

#this saves only selected species (here the species we also plot)
saveChemistryOutput(output_dir + '/chemistry_select.dat', 
                    temperature, pressure, 
                    output_data.total_element_density, 
                    output_data.mean_molecular_weight, 
                    output_data.number_densities, 
                    fastchem,
                    plot_species)


#save the monitor output to a file
#here, the data is saved as a pandas DataFrame inside a pickle file
# saveMonitorOutputPandas(output_dir + '/monitor.pkl', 
#                   temperature, pressure, 
#                   output_data.element_conserved,
#                   output_data.fastchem_flag,
#                   output_data.nb_iterations,
#                   output_data.nb_chemistry_iterations,
#                   output_data.nb_cond_iterations,
#                   output_data.total_element_density,
#                   output_data.mean_molecular_weight,
#                   fastchem)

#this would save the output of all species
#here, the data is saved as a pandas DataFrame inside a pickle file
# saveChemistryOutputPandas(output_dir + '/chemistry.pkl',
#                     temperature, pressure,
#                     output_data.total_element_density,
#                     output_data.mean_molecular_weight,
#                     output_data.number_densities,
#                     fastchem)



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


#convert the output into a numpy array
number_densities = np.array(output_data.number_densities)


#total gas particle number density from the ideal gas law 
#used later to convert the number densities to mixing ratios
gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


#and plot...
for i in range(0, len(plot_species_symbols)):
  fig = plt.plot(number_densities[:, plot_species_indices[i]]/gas_number_density, pressure)

plt.xscale('log')
plt.yscale('log')
plt.gca().set_ylim(plt.gca().get_ylim()[::-1])

plt.xlabel("Mixing ratios")
plt.ylabel("Pressure (bar)")
plt.legend(plot_species_symbols)

plt.show()

#we could also save the figure as a pdf
#plt.savefig(output_dir + '/fastchem_fig.pdf')
