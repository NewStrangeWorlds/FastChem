
import pyfastchem
from save_output import saveChemistryOutput, saveCondOutput, saveMonitorOutput, saveChemistryOutputPandas, saveMonitorOutputPandas, saveCondOutputPandas
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import constants as const


#we read in a p-T structure for a brown dwarf
#data = np.loadtxt("../input/Gliese_229b.dat")

# temperature = np.arange(6000.0, 75, -50.0)
# temperature = np.append(temperature, 75.0)

temperature = np.array([6000.0, 5950.0, 5900.0, 5850.0, 5800.0, 5750.0, 5700.0, 5650.0, 5600.0, 5550.0, 5500.0, 5450.0, 5400.0, 5350.0, 5300.0, 5250.0, 5200.0, 5150.0, 5100.0, 5050.0, 5000.0, 4950.0, 4900.0, 4850.0, 4800.0, 4750.0, 4700.0, 4650.0, 4600.0, 4550.0, 4500.0, 4450.0, 4400.0, 4350.0, 4300.0, 4250.0, 4200.0, 4150.0, 4100.0, 4050.0, 4000.0, 3950.0, 3900.0, 3850.0, 3800.0, 3750.0, 3700.0, 3650.0, 3600.0, 3550.0, 3500.0, 3450.0, 3400.0, 3350.0, 3300.0, 3250.0, 3200.0, 3150.0, 3100.0, 3050.0, 3000.0, 2950.0, 2900.0, 2850.0, 2800.0, 2750.0, 2700.0, 2650.0, 2600.0, 2550.0, 2500.0, 2450.0, 2400.0, 2350.0, 2300.0, 2250.0, 2200.0, 2150.0, 2100.0, 2050.0, 2000.0, 1950.0, 1900.0, 1850.0, 1800.0, 1750.0, 1700.0, 1650.0, 1600.0, 1550.0, 1500.0, 1450.0, 1400.0, 1350.0, 1300.0, 1250.0, 1200.0, 1150.0, 1100.0, 1050.0, 1000.0, 950.0, 900.0, 850.0, 800.0, 750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 475.0, 450.0, 425.0, 400.0, 375.0, 350.0, 325.0, 300.0, 275.0, 250.0, 225.0, 200.0, 175.0, 150.0, 125.0, 100.0, 75.0])

pressure = np.full(temperature.shape, 10**-5.66666666)

#and extract temperature and pressure values

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
  '../input/element_abundances_solar.dat', 
  '../input/logK.dat', 
  '../input/logK_condensates_ge_li.dat',
  1)


#create the input and output structures for FastChem
input_data = pyfastchem.FastChemInput()
output_data = pyfastchem.FastChemOutput()

input_data.temperature = temperature
input_data.pressure = pressure

#use equilibrium condensation
input_data.equilibrium_condensation = True
#fastchem.setParameter('condSolveFullSystem', np.bool_(True))
fastchem.setParameter('minDensityExponentElement', -3900.0)
fastchem.setParameter('accuracyChem', 1e-6)
fastchem.setParameter('accuracyCond', 1e-6)
#fastchem.setParameter('condTau', 1e-30)

#this would turn on the rainout condensation approach
input_data.rainout_condensation = True


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
