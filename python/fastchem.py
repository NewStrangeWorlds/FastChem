 
import pyfastchem
from save_output import saveChemistryOutput, saveMonitorOutput
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

#some input values for temperature and pressure
temperature = np.full(1000, 500)
pressure = np.logspace(-6, 1, num=1000)


#the chemical species we want to plot later
#note that FastChem uses the Hill notation
plot_species = ['H2O1', 'C1O2', 'C1O1', 'C1H4', 'H3N1', 'GRRR']
#for the plot lables, we therefore use separate strings in the usual notation
plot_species_lables = ['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'GRRR']


#create a FastChem object
#fastchem = pyfastchem.FastChem('input/parameters.dat', 1)
fastchem = pyfastchem.FastChem('input/element_abundances_jo_orich.dat', 'input/logK_full_new.dat', 1)


#create the input and output structures for FastChem
input_data = pyfastchem.FastChemInput()
output_data = pyfastchem.FastChemOutput()

input_data.temperature = temperature
input_data.pressure = pressure

#run FastChem on the entire p-T structure
fastchem_flag = fastchem.calcDensities(input_data, output_data)
print("FastChem reports:", pyfastchem.FASTCHEM_MSG[fastchem_flag])


#convert the output into a numpy array
number_densities = np.array(output_data.number_densities)


#total gas particle number density from the ideal gas law 
#used later to convert the number densities to mixing ratios
gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)



#save the monitor output to a file
saveMonitorOutput('output/monitor_test.dat', temperature, pressure, output_data, fastchem)

#this would save the output of all species
#saveChemistryOutput('output/test.dat', temperature, pressure, output_data, fastchem)

#this saves only selected species (here the species we also plot)
saveChemistryOutput('output/test.dat', temperature, pressure, output_data, fastchem, plot_species)



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
  fig = plt.plot(number_densities[:, plot_species_indices[i]]/gas_number_density, pressure)

plt.xscale('log')
plt.yscale('log')
plt.gca().set_ylim(plt.gca().get_ylim()[::-1])

plt.xlabel("Mixing ratios")
plt.ylabel("Pressure (bar)")
plt.legend(plot_species_symbols)
plt.title("T=500 K")

#plt.savefig('eqchem_500K.pdf')
plt.show()
