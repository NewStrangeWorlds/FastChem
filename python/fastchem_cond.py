
import pyfastchem
from save_output import saveChemistryOutput, saveCondOutput, saveMonitorOutput, saveChemistryOutputPandas, saveMonitorOutputPandas, saveCondOutputPandas
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import constants as const


#we read in a p-T structure for a brown dwarf
data = np.loadtxt("../input/example_p_t_structures/Brown_dwarf_Sonora.dat")

#and extract temperature and pressure values
temperature = data[:,1]
pressure = data[:,0]


#define the directory for the output
#here, we currently use the standard one from FastChem
output_dir = '../output'


#the chemical species we want to plot later
#note that the standard FastChem input files use the Hill notation
plot_species = ['H2O1', 'C1O2', 'C1O1', 'C1H4', 'H3N1', 'Fe1S1', 'H2S1']
#for the plot labels, we therefore use separate strings in the usual notation
plot_species_labels = ['H2O', 'CO2', 'CO', 'CH4', 'NH3', 'FeS', 'H2S']

#the default condensate data doesn't use the Hill notation
plot_species_cond = ['Fe(s,l)', 'FeS(s,l)', 'MgSiO3(s,l)', 'Mg2SiO4(s,l)']


#create a FastChem object
fastchem = pyfastchem.FastChem(
  '../input/element_abundances/asplund_2009.dat', 
  '../input/logK/logK.dat',
  '../input/logK/logK_condensates.dat',
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

#use equilibrium condensation
input_data.equilibrium_condensation = True

#this would turn on the rainout condensation approach
#input_data.rainout_condensation = True


#run FastChem on the entire p-T structure
fastchem_flag = fastchem.calcDensities(input_data, output_data)


#convergence summary report
print("FastChem reports:")
print("  -", pyfastchem.FASTCHEM_MSG[fastchem_flag])

if np.amin(output_data.element_conserved[:]) == 1:
  print("  - element conservation: ok")
else:
  print("  - element conservation: fail")


#convert the output into a numpy array
number_densities = np.array(output_data.number_densities)
number_densities_cond = np.array(output_data.number_densities_cond)


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

#this would save the output of all species
saveCondOutput(output_dir + '/condensates.dat', 
               temperature, pressure, 
               output_data.element_cond_degree, 
               output_data.number_densities_cond, 
               fastchem)


#save the monitor output to a file
#here, the data is saved as a pandas DataFrame inside a pickle file
# saveMonitorOutputPandas(output_dir + '/monitor.pkl', 
#                   temperature, pressure, 
#                   output_data.element_conserved,
#                   output_data.fastchem_flag,
#                   output_data.nb_chemistry_iterations,
#                   output_data.total_element_density,
#                   output_data.mean_molecular_weight,
#                   fastchem)

# #this would save the output of all species
# #here, the data is saved as a pandas DataFrame inside a pickle file
# saveChemistryOutputPandas(output_dir + '/chemistry.pkl',
#                     temperature, pressure,
#                     output_data.total_element_density,
#                     output_data.mean_molecular_weight,
#                     output_data.number_densities,
#                     fastchem)


#this would save the condensate output
#the data is saved as a pandas DataFrame inside a pickle file
# saveCondOutputPandas(output_dir + '/condensates.pkl',
#                     temperature, pressure,
#                     output_data.element_cond_degree,
#                     output_data.number_densities_cond,
#                     fastchem)



#total gas particle number density from the ideal gas law 
#used later to convert the number densities to mixing ratios
gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


#check the gas-phase species we want to plot and get their indices from FastChem
plot_species_indices = []
plot_species_symbols = []

for i, species in enumerate(plot_species):
  index = fastchem.getGasSpeciesIndex(species)

  if index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
    plot_species_indices.append(index)
    plot_species_symbols.append(plot_species_labels[i])
  else:
    print("Species", species, "to plot not found in FastChem")


#check the condensates we want to plot and get their indices from FastChem
plot_species_indices_cond = []
plot_species_symbols_cond = []

for i, species in enumerate(plot_species_cond):
  index = fastchem.getCondSpeciesIndex(species)

  if index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
    plot_species_indices_cond.append(index)
    plot_species_symbols_cond.append(plot_species_cond[i])
  else:
    print("Species", species, "to plot not found in FastChem")


#and plot everything...
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

for i in range(0, len(plot_species_symbols)):
  ax1.plot(number_densities[:, plot_species_indices[i]]/gas_number_density, pressure)
  
ax1.set_xlim(1e-10,1.1)
ax1.set_title("Gas phase")

ax1.set(xscale='log', yscale = 'log')
ax1.set_ylim(ax1.get_ylim()[::-1])

ax1.set(xlabel = 'Mixing ratios', ylabel = 'Pressure (bar)')
ax1.legend(plot_species_symbols)


for i in range(0, len(plot_species_symbols_cond)):
  ax2.plot(number_densities_cond[:, plot_species_indices_cond[i]], pressure)

ax2.set_xlim(left=1e5)
ax2.set_title("Condensates")

ax2.set(xscale='log', yscale = 'log')
ax2.set_ylim(ax2.get_ylim()[::-1])

ax2.set(xlabel = 'Number density (cm$^{-3}$)', ylabel = 'Pressure (bar)')
ax2.legend(plot_species_symbols_cond)


plt.show()

#we could also save the figure as a pdf
#plt.savefig(output_dir + '/fastchem_fig.pdf')
