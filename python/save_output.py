import pyfastchem
import numpy as np
from astropy import constants as const
import pandas as pd


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
    add_columns = np.atleast_1d(additional_columns)

    if add_columns.size > add_columns.shape[0]:
      nb_add_columns = add_columns.shape[0]
    else:
      nb_add_columns = 1
  

  add_columns_desc = np.atleast_1d(additional_columns_desc)


  if (nb_add_columns > 0 and add_columns_desc.size != nb_add_columns) or (additional_columns_desc is None and nb_add_columns > 0):
    print('Warning from saveChemistryOutput: The number of additional column descriptions does not match the number of data columns.')


  #total gas particle number density from the ideal gas law 
  #used to convert the number densities to mixing ratios
  gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


  #calculate the mixing ratios from the number densities
  mixing_ratios = np.array(number_densities) / gas_number_density.value[:, None]


  with open(file_path, 'w') as file:
    #file header
    file.write('{0:<16}\t{1:<16}\t{2:<16}\t{3:<16}\t{4:<16}'.format('#p (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (u)'))

    #print the description of the additional columns
    #use 'unk' if the their number do not correspond to the number of additonal columns
    for i in range (nb_add_columns):
      if add_columns_desc.size == nb_add_columns:
        file.write('\t{0:<16}'.format(add_columns_desc[i]))
      else:
        file.write('\t{0:<16}'.format('unk'))


    if output_species is None: #save all species
      for j in range(fastchem.getGasSpeciesNumber()):
        file.write('\t{0:<16}'.format(fastchem.getGasSpeciesSymbol(j)))
    
      file.write('\n')
    
      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}\t{2:1.10e}\t{3.value:1.10e}\t{4:1.10e}'.format(pressure[i], 
                                                                                        temperature[i], 
                                                                                        total_element_density[i], 
                                                                                        gas_number_density[i], 
                                                                                        mean_molecular_weight[i]))


        #print the additional columns
        for j in range(nb_add_columns):
          if nb_add_columns == 1:
            file.write('\t{0:1.10e}'.format(additional_columns[i]))
          else:
            file.write('\t{0:1.10e}'.format(additional_columns[j, i]))

      
        for j in range(fastchem.getGasSpeciesNumber()):
          file.write('\t{0:1.10e}'.format(mixing_ratios[i, j]))
      
        file.write('\n')
   
    else:  #save only selected species
      for species in output_species:
        if fastchem.getGasSpeciesIndex(species) != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
          file.write('\t{0:<16}'.format(species))
        else:
          print('Species ', species, ' not found during saving of the chemistry output!')
    
      file.write('\n')
    
      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}\t{2:1.10e}\t{3.value:1.10e}\t{4:1.10e}'.format(pressure[i], 
                                                                                        temperature[i], 
                                                                                        total_element_density[i], 
                                                                                        gas_number_density[i], 
                                                                                        mean_molecular_weight[i]))


        #print the additional columns
        for j in range(nb_add_columns):
          if nb_add_columns == 1:
            file.write('\t{0:1.10e}'.format(additional_columns[i]))
          else:
            file.write('\t{0:1.10e}'.format(additional_columns[j, i]))


        for species in output_species:
          species_index = fastchem.getGasSpeciesIndex(species)
          if species_index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
            file.write('\t{0:1.10e}'.format(mixing_ratios[i, species_index]))
      
        file.write('\n')


  return None



#Saves the FastChem output in an ASCII file
def saveCondOutput(file_path,                     #the path to the output file
                   temperature, pressure,         #arrays of temperature and pressure
                   element_cond_degree,           #degree of condensation for the elements
                   number_densities,              #2D array of number densities
                   fastchem,                      #the FastChem object
                   output_species=None,           #optional array with symbols of species that should be saved
                   additional_columns=None,       #optional, additional columns for the output file
                   additional_columns_desc=None): #the header descriptions of the additional columns
  

  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    add_columns = np.atleast_1d(additional_columns)

    if add_columns.size > add_columns.shape[0]:
      nb_add_columns = add_columns.shape[0]
    else:
      nb_add_columns = 1
  

  add_columns_desc = np.atleast_1d(additional_columns_desc)


  if (nb_add_columns > 0 and add_columns_desc.size != nb_add_columns) or (additional_columns_desc is None and nb_add_columns > 0):
    print('Warning from saveCondOutput: The number of additional column descriptions does not match the number of data columns.')

  element_cond_degree = np.array(element_cond_degree)
  number_densities = np.array(number_densities)


  with open(file_path, 'w') as file:
    #file header
    file.write('{0:<16}\t{1:<16}'.format('#p (bar)', 'T (K)'))

    #print the description of the additional columns
    #use 'unk' if the their number do not correspond to the number of additonal columns
    for i in range (nb_add_columns):
      if add_columns_desc.size == nb_add_columns:
        file.write('\t{0:<16}'.format(add_columns_desc[i]))
      else:
        file.write('\t{0:<16}'.format('unk'))
      
    #for the degree of condensation
    for j in range(fastchem.getElementNumber()):
      file.write('\t{0:<16}'.format(fastchem.getElementSymbol(j)))


    if output_species is None: #save all species
      for j in range(fastchem.getCondSpeciesNumber()):
        file.write('\t{0:<16}'.format(fastchem.getCondSpeciesSymbol(j)))

      file.write('\n')

      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}'.format(pressure[i], temperature[i]))


        #print the additional columns
        for j in range(nb_add_columns):
          if nb_add_columns == 1:
            file.write('\t{0:1.10e}'.format(additional_columns[i]))
          else:
            file.write('\t{0:1.10e}'.format(additional_columns[j, i]))

        #the degree of condensation
        for j in range(fastchem.getElementNumber()):
          file.write('\t{0:1.10e}'.format(element_cond_degree[i, j]))
        
        #the condensates
        for j in range(fastchem.getCondSpeciesNumber()):
          file.write('\t{0:1.10e}'.format(number_densities[i, j]))
      
        file.write('\n')
   
    else:  #save only selected species
      for species in output_species:
        if fastchem.getCondSpeciesIndex(species) != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
          file.write('\t{0:<16}'.format(species))
        else:
          print('Species ', species, ' not found during saving of the chemistry output!')
    
      file.write('\n')
    
      #and the chemistry output
      for i in range(pressure.size):
        file.write('{0:1.10e}\t{1:1.10e}'.format(pressure[i], temperature[i]))


        #print the additional columns
        for j in range(nb_add_columns):
          if nb_add_columns == 1:
            file.write('\t{0:1.10e}'.format(additional_columns[i]))
          else:
            file.write('\t{0:1.10e}'.format(additional_columns[j, i]))

        #the degree of condensation
        for j in range(fastchem.getElementNumber()):
          file.write('\t{0:1.10e}'.format(element_cond_degree[i, j]))

        #the selected condensates
        for species in output_species:
          species_index = fastchem.getCondSpeciesIndex(species)
          if species_index != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
            file.write('\t{0:1.10e}'.format(number_densities[i, species_index]))

        file.write('\n')


  return None




#Saves the FastChem monitor output in an ASCII file
def saveMonitorOutput(file_path,                     #the path to the output file
                      temperature, pressure,         #arrays of temperature and pressure
                      element_conserved,             #array of int for the element conservation
                      fastchem_flags,                #array of int for the output flags
                      nb_iterations,                 #array of int with the number of iterations
                      nb_chemistry_iterations,       #array of int for the number of chemistry iterations
                      nb_condensation_iterations,    #array of int for the number of condensation iterations
                      total_element_density,         #array with total element density
                      mean_molecular_weight,         #array with the mean molecular weight
                      fastchem,                      #the FastChem object
                      additional_columns=None,       #optional, additional columns for the output file
                      additional_columns_desc=None): #the header descriptions of the additional columns


  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    add_columns = np.atleast_1d(additional_columns)

    if add_columns.size > add_columns.shape[0]:
      nb_add_columns = add_columns.shape[0]
    else:
      nb_add_columns = 1


  add_columns_desc = np.atleast_1d(additional_columns_desc)


  if (nb_add_columns > 0 and add_columns_desc.size != nb_add_columns) or (additional_columns_desc is None and nb_add_columns > 0):
    print('Warning from saveMonitorOutput: The number of additional column descriptions does not match the number of data columns.')


  #total gas particle number density from the ideal gas law 
  gas_number_density = pressure*1e6 / (const.k_B.cgs * temperature)


  #convert the 2D element conservation array into a NumPy array
  element_conserved = np.array(element_conserved)


  with open(file_path, 'w') as file:
    #file header
    file.write('{0:<16}{1:<16}{2:<16}{3:<16}{4:<24}{5:<24}{6:<24}{7:<24}{8:<24}{9:<24}{10:<24}'.format(
      '#grid point', 
      '#iterations', 
      '#chem_iter', 
      '#cond_iter',
      'converged', 
      'elem_conserved', 
      'p (bar)', 
      'T (K)', 
      'n_<tot> (cm-3)', 
      'n_g (cm-3)', 
      'm (u)'))

    #print the header description of the additional columns
    #use 'unk' if the their number do not correspond to the number of additonal columns
    for i in range (nb_add_columns):
      if add_columns_desc.size == nb_add_columns:
        file.write('{0:<24}'.format(add_columns_desc[i]))
      else:
        file.write('{0:<24}'.format('unk'))

    #header with element symbols
    for j in range(fastchem.getElementNumber()):
      file.write('{0:<8}'.format(fastchem.getElementSymbol(j)))
    
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

      file.write('{0:<16d}{1:<16d}{2:<16d}{3:<16d}{4:<24s}{5:<24s}{6:<24.10e}{7:<24.10e}{8:<24.10e}{9.value:<24.10e}{10:<24.10e}'.format(
        i, 
        nb_iterations[i],
        nb_chemistry_iterations[i], 
        nb_condensation_iterations[i], 
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



#Saves the FastChem output in a pickle file using pandas
def saveChemistryOutputPandas(file_path,                     #the path to the output file
                              temperature, pressure,         #arrays of temperature and pressure
                              total_element_density,         #array of total element density 
                              mean_molecular_weight,         #array of mean molecular weights 
                              number_densities,              #2D array of number densities
                              fastchem,                      #the FastChem object
                              output_species=None,           #optional array with symbols of species that should be saved
                              additional_columns=None,       #optional, additional columns for the output file
                              additional_columns_desc=None): #the header descriptions of the additional columns
  
  # total gas particle number density from the ideal gas law
  # used to convert the number densities to mixing ratios
  gas_number_density = pressure * 1e6 / (const.k_B.cgs * temperature)

  # calculate the mixing ratios from the number densities
  mixing_ratios = np.array(number_densities) / gas_number_density.value[:, None]


  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    additional_columns = np.atleast_1d(additional_columns)

    if additional_columns.size > additional_columns.shape[0]:
      nb_add_columns = additional_columns.shape[0]
    else:
      nb_add_columns = 1

  add_columns_desc = np.atleast_1d(additional_columns_desc)


  if (add_columns_desc.size != nb_add_columns and nb_add_columns > 0) or (additional_columns_desc is None and nb_add_columns > 0):
    print('Warning from saveChemistryOutputPandas: The number of additional column descriptions does not match the number of data columns.')


  #general column headers
  columns = ['p (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (u)']


  #add the descriptions of the additional columns to the header
  #if their number does not equal the number of additional columns, add 'unk'
  for i in range(nb_add_columns):
    if add_columns_desc.size == nb_add_columns and additional_columns_desc is not None:
      columns.append(add_columns_desc[(i)])
    else:
      columns.append('unk')


  #header for species symbols
  if output_species is None:
    for j in range(fastchem.getGasSpeciesNumber()):
      columns.append(fastchem.getGasSpeciesSymbol(j))
  else:
    select_species_id = []

    for species in output_species:
      if fastchem.getGasSpeciesIndex(species) != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
        columns.append(species)
        select_species_id.append(fastchem.getGasSpeciesIndex(species))
      else:
        print('Species ', species, ' not found during saving of the chemistry output!')


  #combine all general data columns
  if additional_columns is None:
    rows = np.vstack([pressure,
                      temperature,
                      total_element_density,
                      gas_number_density.value,
                      mean_molecular_weight]
                    ).T
  else:
    rows = np.vstack([pressure,
                      temperature,
                      total_element_density,
                      gas_number_density.value,
                      mean_molecular_weight,
                      additional_columns]
                    ).T


  #and now we add the mixing ratios
  if output_species is None:
    data = np.hstack([rows, mixing_ratios])
  else:
    data = np.hstack([rows, mixing_ratios[:,select_species_id]])


  #combine column headers and data into a pandas DataFrame
  df = pd.DataFrame(data=data, columns=columns)


  #and finally save it as a pickle file
  df.to_pickle(file_path)

  return None



#Saves the FastChem output in a pickle file using pandas
def saveCondOutputPandas(file_path,                     #the path to the output file
                         temperature, pressure,         #arrays of temperature and pressure
                         element_cond_degree,           #degree of condensation for the elements
                         number_densities_cond,         #2D array of number densities
                         fastchem,                      #the FastChem object
                         output_species=None,           #optional array with symbols of species that should be saved
                         additional_columns=None,       #optional, additional columns for the output file
                         additional_columns_desc=None): #the header descriptions of the additional columns

  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    additional_columns = np.atleast_1d(additional_columns)

    if additional_columns.size > additional_columns.shape[0]:
      nb_add_columns = additional_columns.shape[0]
    else:
      nb_add_columns = 1

  add_columns_desc = np.atleast_1d(additional_columns_desc)


  if (add_columns_desc.size != nb_add_columns and nb_add_columns > 0) or (additional_columns_desc is None and nb_add_columns > 0):
    print('Warning from saveChemistryOutputPandas: The number of additional column descriptions does not match the number of data columns.')


  #general column headers
  columns = ['P (bar)', 'T (K)']


  #add the descriptions of the additional columns to the header
  #if their number does not equal the number of additional columns, add 'unk'
  for i in range(nb_add_columns):
    if add_columns_desc.size == nb_add_columns and additional_columns_desc is not None:
      columns.append(add_columns_desc[(i)])
    else:
      columns.append('unk')

  #header for the degree of condensation
  for j in range(fastchem.getElementNumber()):
    columns.append(fastchem.getElementSymbol(j))


  #header for species symbols
  if output_species is None:
    for j in range(fastchem.getCondSpeciesNumber()):
      columns.append(fastchem.getCondSpeciesSymbol(j))
  else:
    select_species_id = []

    for species in output_species:
      if fastchem.getCondSpeciesIndex(species) != pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
        columns.append(species)
        select_species_id.append(fastchem.getCondSpeciesIndex(species))
      else:
        print('Species ', species, ' not found during saving of the chemistry output!')


  #combine all general data columns
  if additional_columns is None:
    rows = np.vstack([pressure,
                      temperature]
                    ).T
  else:
    rows = np.vstack([pressure,
                      temperature,
                      additional_columns]
                    ).T


  #and now we add the mixing ratios
  if output_species is None:
    data = np.hstack([rows, element_cond_degree, number_densities_cond])
  else:
    data = np.hstack([rows, element_cond_degree, number_densities_cond[:,select_species_id]])


  #combine column headers and data into a pandas DataFrame
  df = pd.DataFrame(data=data, columns=columns)


  #and finally save it as a pickle file
  df.to_pickle(file_path)

  return None



#Saves the FastChem monitor output in a pickle file using pandas
def saveMonitorOutputPandas(file_path,                     #the path to the output file
                            temperature, pressure,         #arrays of temperature and pressure
                            element_conserved,             #array of int for the element conservation
                            fastchem_flags,                #array of int for the output flags
                            nb_iterations,                 #array of int with the number of iterations
                            nb_chemistry_iterations,       #array of int for the number of chemistry iterations
                            nb_condensation_iterations,    #array of int for the number of condensation iterations
                            total_element_density,         #array with total element density
                            mean_molecular_weight,         #array with the mean molecular weight
                            fastchem,                      #the FastChem object
                            additional_columns=None,       #optional, additional columns for the output file
                            additional_columns_desc=None): #the header descriptions of the additional columns
  
  # total gas particle number density from the ideal gas law
  # used to convert the number densities to mixing ratios
  gas_number_density = pressure * 1e6 / (const.k_B.cgs * temperature)

  #convert the 2D element conservation array into a NumPy array
  element_conserved = np.array(element_conserved)


  #find out how many additional columns and descriptions we have
  nb_add_columns = 0

  if additional_columns is not None:
    additional_columns = np.atleast_1d(additional_columns)

    if additional_columns.size > additional_columns.shape[0]:
      nb_add_columns = additional_columns.shape[0]
    else:
      nb_add_columns = 1

  add_columns_desc = np.atleast_1d(additional_columns_desc)


  if (add_columns_desc.size != nb_add_columns and nb_add_columns > 0) or (additional_columns_desc is None and nb_add_columns > 0):
    print('Warning from saveMonitorOutputPandas: The number of additional column descriptions does not match the number of data columns.')


  #general column headers
  columns = ['iterations', 'chem_iterations', 'cond_iterations', 'c_convergence', 'elem_conserved', 'p (bar)', 'T (K)', 'n_<tot> (cm-3)', 'n_g (cm-3)', 'm (u)']


  #add the descriptions of the additional columns to the header
  #if their number does not equal the number of additional columns, add 'unk'
  for i in range(nb_add_columns):
    if add_columns_desc.size == nb_add_columns and additional_columns_desc is not None:
      columns.append(add_columns_desc[(i)])
    else:
      columns.append('unk')


  nb_elements = fastchem.getElementNumber()

  #header for element symbols
  for j in range(nb_elements):
    columns.append(fastchem.getElementSymbol(j))


  #and the debug output
  all_elements_conserved = np.ones(len(fastchem_flags), dtype=int)

  for i in range(len(fastchem_flags)):
    if np.isin(0, element_conserved[i], assume_unique=True) : all_elements_conserved[i] = 0

  #combine all general data columns
  if additional_columns is None:
    rows = np.vstack([nb_iterations,
                      nb_chemistry_iterations,
                      nb_condensation_iterations,
                      fastchem_flags,
                      all_elements_conserved,
                      pressure,
                      temperature,
                      total_element_density,
                      gas_number_density.value,
                      mean_molecular_weight]
                    ).T
  else:
    rows = np.vstack([nb_iterations,
                      nb_chemistry_iterations,
                      nb_condensation_iterations,
                      fastchem_flags,
                      all_elements_conserved,
                      pressure,
                      temperature,
                      total_element_density,
                      gas_number_density.value,
                      mean_molecular_weight,
                      additional_columns]
                    ).T

  #and now we add the element conservation info
  data = np.hstack([rows, element_conserved])


  #combine column headers and data into a pandas DataFrame
  df = pd.DataFrame(data=data, columns=columns)

  #change some of the datatypes back to int
  df = df.astype({'iterations':'int64'})
  df = df.astype({'chem_iterations':'int64'})
  df = df.astype({'cond_iterations':'int64'})
  df = df.astype({'c_convergence':'int64'})
  df = df.astype({'elem_conserved':'int64'})

  for j in range(nb_elements):
    df = df.astype({fastchem.getElementSymbol(j) :'int64'})

  #and finally save it as a pickle file
  df.to_pickle(file_path)

  return None



