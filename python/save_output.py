import pyfastchem
import numpy as np
from astropy import constants as const


#Saves the FastChem output in an ASCII file
def saveChemistryOutput(file_path,                     #the path to the output file
                        temperature, pressure,         #arrays of temperature and pressure
                        total_element_density,         #array of total element density 
                        mean_molecular_weight,         #array of mean molecular weights 
                        number_densities,              #2D array of number densities
                        fastchem,                      #the FastChem object
                        output_species=None,           #optional array with symbols of species that should be saved
                        additional_columns=None,       #optional, additional columns for the output file
                        additional_columns_desc=None): #the header descriptions of the additional columns
  

  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    add_columns = np.array(additional_columns)
    add_columns_desc = np.array(additional_columns_desc)

    if add_columns.size > add_columns.shape[0]:
      nb_add_columns = add_columns.shape[0]
    else:
      nb_add_columns = 1


  #total gas particle number density from the ideal gas law 
  #used to convert the number densities to mixing ratios
  gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


  #calculate the mixing ratios from the number densities
  mixing_ratios = np.array(number_densities) / gas_number_density.value[:, None]


  with open(file_path, 'w') as file:
    #file header
    file.write('{0:<16}\t{1:<16}\t{2:<16}\t{3:<16}\t{4:<16}'.format('#P (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (g/mol)'))

    #print the description of the additional columns
    #use 'unk' if the their number do not correspond to the number of additonal columns
    for i in range (nb_add_columns):
      if add_columns_desc.size == nb_add_columns:
        if (add_columns_desc.size == 1):
          file.write('\t{0:<16}'.format(additional_columns_desc))
        else:
          file.write('\t{0:<16}'.format(additional_columns_desc[i]))
      else:
        file.write('\t{0:<16}'.format('unk'))


    if output_species is None: #save all species
      for j in range(fastchem.getSpeciesNumber()):
        file.write('\t{0:<16}'.format(fastchem.getSpeciesSymbol(j)))
    
      file.write('\n')
    
      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}\t{2:1.10e}\t{3.value:1.10e}\t{4:1.10e}'.format(pressure[i], temperature[i], total_element_density[i], gas_number_density[i], mean_molecular_weight[i]))


        #print the additional columns
        for j in range(nb_add_columns):
          if nb_add_columns == 1:
            file.write('\t{0:1.10e}'.format(additional_columns[i]))
          else:
            file.write('\t{0:1.10e}'.format(additional_columns[j, i]))

      
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
        file.write('{0:1.10e}\t{1:1.10e}\t{2:1.10e}\t{3.value:1.10e}\t{4:1.10e}'.format(pressure[i], temperature[i], total_element_density[i], gas_number_density[i], mean_molecular_weight[i]))


        #print the additional columns
        for j in range(nb_add_columns):
          if nb_add_columns == 1:
            file.write('\t{0:1.10e}'.format(additional_columns[i]))
          else:
            file.write('\t{0:1.10e}'.format(additional_columns[j, i]))


        for species in output_species:
          species_index = fastchem.getSpeciesIndex(species)
          if species_index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
            file.write('\t{0:1.10e}'.format(mixing_ratios[i, species_index]))
      
        file.write('\n')


  return None




#Saves the FastChem monitor output in an ASCII file
def saveMonitorOutput(filename,                      #the path to the output file
                      temperature, pressure,         #arrays of temperature and pressure
                      element_conserved,             #array if int for the element conservation
                      fastchem_flags,                #array of int for the output flags
                      nb_chemistry_iterations,       #array of int for the number of iterations
                      total_element_density,         #array with total element density
                      mean_molecular_weight,         #array with the mean molecular weight
                      fastchem,                      #the FastChem object
                      additional_columns=None,       #optional, additional columns for the output file
                      additional_columns_desc=None): #the header descriptions of the additional columns


  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    add_columns = np.array(additional_columns)
    add_columns_desc = np.array(additional_columns_desc)

    if add_columns.size > add_columns.shape[0]:
      nb_add_columns = add_columns.shape[0]
    else:
      nb_add_columns = 1


  #total gas particle number density from the ideal gas law 
  gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


  #convert the 2D element conservation array into a NumPy array
  element_conserved = np.array(element_conserved)


  with open(filename, 'w') as file:
    #file header
    file.write('{0:<16}{1:<16}{2:<24}{3:<24}{4:<24}{5:<24}{6:<24}{7:<24}{8:<24}'.format('#grid point', 'c_iterations', 'c_convergence', 'elem_conserved', 'P (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (g/mol)'))

    #print the header description of the additional columns
    #use 'unk' if the their number do not correspond to the number of additonal columns
    for i in range (nb_add_columns):
      if add_columns_desc.size == nb_add_columns:
        if (add_columns_desc.size == 1):
          file.write('{0:<24}'.format(additional_columns_desc))
        else:
          file.write('{0:<24}'.format(additional_columns_desc[i]))
      else:
        file.write('{0:<24}'.format('unk'))

    #header with element symbols
    for j in range(fastchem.getElementNumber()):
      file.write('{0:<8}'.format(fastchem.getSpeciesSymbol(j)))
    
    file.write('\n')


    output_flags = ['fail', 'ok']
    
    #and the debug output
    for i in range(len(fastchem_flags)):

      if fastchem_flags[i] == pyfastchem.FASTCHEM_SUCCESS:
        c_conv = output_flags[1]
      else:
        c_conv = output_flags[0]

      all_elements_conserved = 1

      if np.isin(0, element_conserved[i], assume_unique=True) : all_elements_conserved = 0

      file.write('{0:<16d}{1:<16d}{2:<24s}{3:<24s}{4:<24.10e}{5:<24.10e}{6:<24.10e}{7.value:<24.10e}{8:<24.10e}'.format(i, 
                                                                                                                  nb_chemistry_iterations[i], 
                                                                                                                  c_conv, 
                                                                                                                  output_flags[all_elements_conserved], 
                                                                                                                  pressure[i], temperature[i], 
                                                                                                                  total_element_density[i], 
                                                                                                                  gas_number_density[i], 
                                                                                                                  mean_molecular_weight[i]))

      #print the additional columns
      for j in range(nb_add_columns):
        if nb_add_columns == 1:
          file.write('{0:<24.10e}'.format(additional_columns[i]))
        else:
          file.write('{0:<24.10e}'.format(additional_columns[j, i]))

      #print element conservation info
      for j in range(fastchem.getElementNumber()):
        file.write('{0:<8s}'.format(output_flags[element_conserved[i,j]]))
      
      file.write('\n')


  return None
