
import pyfastchem
from save_output import saveChemistryOutput, saveCondOutput, saveMonitorOutput, saveChemistryOutputPandas, saveMonitorOutputPandas, saveCondOutputPandas
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import constants as const


met = 1000
co = 1.5


temp2 = np.linspace(70, 500, 44)
temp3 = np.linspace(525, 800, 12)
temp4 = np.linspace(900, 6000, 52)

temperatures = temp2
temperatures = np.append(temperatures, temp3)
temperatures = np.append(temperatures, temp4)

#use the pressures from Simon's opacities to avoid interpolation
pressures = np.array([10**-6.0, \
                      10**-5.66666666, 10**-5.33333333, 10**-5.0, \
                      10**-4.66666666, 10**-4.33333333, 10**-4.0, \
                      10**-3.66666666, 10**-3.33333333, 10**-3.0, \
                      10**-2.66666666, 10**-2.33333333, 10**-2.0, \
                      10**-1.66666666, 10**-1.33333333, 10**-1.0, \
                      10**-0.66666666, 10**-0.33333333, 10**0.0, \
                      10**0.33333333, 10**0.66666666, 10**1.0, \
                      10**1.33333333, 10**1.66666666, 10**2.0, \
                      10**2.33333333, 10**2.66666666, 10**3.0])

    
nb_temperatures = temperatures.size
nb_pressures = pressures.size

p_t_grid = np.zeros((nb_temperatures, nb_pressures, 2))

for i in range(nb_temperatures):
  for j in range(nb_pressures):
    p_t_grid[i,j,:] = np.array([temperatures[i], pressures[j]])


p_t_grid_fc = np.reshape(p_t_grid, (nb_temperatures*nb_pressures, 2))

temperature = p_t_grid_fc[:,0]
pressure = p_t_grid_fc[:,1]


#define the directory for the output
#here, we currently use the standard one from FastChem
output_dir = '../output'


#create a FastChem object
#we could also create a FastChem object by using the parameter file
#note, however, that the file locations in the parameter file are relative
#to the location from where this Python script is called from
fastchem = pyfastchem.FastChem('../input/parameters_py1.dat', 1)


#make a copy of the solar abundances from FastChem
solar_abundances = np.array(fastchem.getElementAbundances())


element_abundances = np.copy(solar_abundances)
  
#scale the element abundances, except those of H and He
for j in range(0, fastchem.getElementNumber()):
  if fastchem.getElementSymbol(j) != 'H' and fastchem.getElementSymbol(j) != 'He':
    element_abundances[j] *= met


index_C = fastchem.getElementIndex('C')
index_O = fastchem.getElementIndex('O')

element_abundances[index_C] = element_abundances[index_O] * co
 
fastchem.setElementAbundances(element_abundances)




#create the input and output structures for FastChem
input_data = pyfastchem.FastChemInput()
output_data = pyfastchem.FastChemOutput()


# temperature =  np.atleast_1d(temperature[1])
# pressure =  np.atleast_1d(pressure[1])


input_data.temperature = temperature
input_data.pressure = pressure


#fastchem.setParameter('minDensityExponentElement', -3900.0)
#fastchem.setParameter('condSolveFullSystem', np.bool_(True))
#fastchem.setParameter('condIterChangeLimit', 10.0)
#fastchem.setParameter('condReduceSystemSize', np.bool_(False))
#fastchem.setParameter('condUseSVD', np.bool_(True))
#fastchem.setParameter('condTau', 1e-5)
#fastchem.setParameter('nbSwitchToNewton', 300)
#fastchem.setParameter('accuracyNewton',1e-12)


#use equilibrium condensation
input_data.equilibrium_condensation = False

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

