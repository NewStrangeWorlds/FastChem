
import pyfastchem
import numpy as np
from astropy import constants as const



def saveChemistryOutput(filename, temperature, pressure, fastchem_output, fastchem, output_species=None):
  
  #total gas particle number density from the ideal gas law 
  #used to convert the number densities to mixing ratios
  gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


  #calculate the mixing ratios from the number densities
  mixing_ratios = np.array(fastchem_output.number_densities) / gas_number_density.value[:, None]


  with open(filename, 'w') as file:
    #file header
    file.write('{0:<16}\t{1:<16}\t{2:<16}\t{3:<16}\t{4:<16}'.format('#P (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (g/mol)'))


    if output_species==None: #save all species
      for j in range(fastchem.getSpeciesNumber()):
        file.write('\t{0:<16}'.format(fastchem.getSpeciesSymbol(j)))
    
      file.write('\n')
    
      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}\t{2:1.10e}\t{3.value:1.10e}\t{4:1.10e}'.format(pressure[i], temperature[i], fastchem_output.total_element_density[i], gas_number_density[i], fastchem_output.mean_molecular_weight[i]))
      
        for j in range(fastchem.getSpeciesNumber()):
          file.write('\t{0:1.10e}'.format(mixing_ratios[i, j]))
      
        file.write('\n')
   
    else:  #save only selected species
      for species in output_species:
        if fastchem.getSpeciesIndex(species) != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
          file.write('\t{0:<16}'.format(species))
        else:
          print('Species ', species, ' not found during saving of the chemistry output!')
    
      file.write('\n')
    
      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}\t{2:1.10e}\t{3.value:1.10e}\t{4:1.10e}'.format(pressure[i], temperature[i], fastchem_output.total_element_density[i], gas_number_density[i], fastchem_output.mean_molecular_weight[i]))
      
        for species in output_species:
          species_index = fastchem.getSpeciesIndex(species)
          if species_index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
            file.write('\t{0:1.10e}'.format(mixing_ratios[i, species_index]))
      
        file.write('\n')

  return None



def saveMonitorOutput(filename, temperature, pressure, fastchem_output, fastchem):
  
  #total gas particle number density from the ideal gas law 
  gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


  #convert the 2D element conservation array into a NumPy array
  element_conserved = np.array(fastchem_output.element_conserved)


  with open(filename, 'w') as file:
    #file header
    file.write('{0:<16}{1:<16}{2:<24}{3:<24}{4:<24}{5:<24}{6:<24}{7:<24}{8:<24}'.format('#grid point', 'c_iterations', 'c_convergence', 'elem_conserved', 'P (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (u)'))
    
    for j in range(fastchem.getElementNumber()):
      file.write('{0:<8}'.format(fastchem.getSpeciesSymbol(j)))
    
    file.write('\n')


    output_flags = ['fail', 'ok']
    
    #and the debug output
    for i in range(pressure.size):

      if fastchem_output.fastchem_flag[i] == pyfastchem.FASTCHEM_SUCCESS:
        c_conv = output_flags[1]
      else:
        c_conv = output_flags[0]

      all_elements_conserved = 1

      if np.isin(0, element_conserved[i], assume_unique=True) : all_elements_conserved = 0

      file.write('{0:<16d}{1:<16d}{2:<24s}{3:<24s}{4:<24.10e}{5:<24.10e}{6:<24.10e}{7.value:<24.10e}{8:<24.10e}'.format(i, 
                                                                                                                  fastchem_output.nb_chemistry_iterations[i], 
                                                                                                                  c_conv, 
                                                                                                                  output_flags[all_elements_conserved], 
                                                                                                                  pressure[i], temperature[i], 
                                                                                                                  fastchem_output.total_element_density[i], 
                                                                                                                  gas_number_density[i], 
                                                                                                                  fastchem_output.mean_molecular_weight[i]))
      
      
          
      for j in range(fastchem.getElementNumber()):
        file.write('{0:<8s}'.format(output_flags[element_conserved[i,j]]))
      
      file.write('\n')


  return None