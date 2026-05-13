
import add_functions

#input files
element_abundance_file = 'element_abundances/asplund_2021.dat'
gas_species_file = 'logK/logK_extended.dat'
cond_species_file = 'logK/logK_condensates_extended.dat'

#output files
new_element_abundance_file = 'element_abundances/custom_abundances.dat'
new_gas_species_file = 'logK/logK_custom.dat'
new_cond_species_file = 'logK/logK_custom_condensates.dat'

#The list of elements to consider
#In this example, we select a very small sample
element_list = ['H', 'He', 'C', 'N', 'O']

#The list above would not include ions. For ions, we have to add the e- explicitly as an element:
#element_list = ['e-', 'H', 'He', 'C', 'N', 'O']



#Read in the element abundances from file and check completeness
element_abundances = add_functions.readElementAbundancesFile(element_abundance_file)
select_element_abundances = add_functions.selectElementAbundances(element_abundances, element_list)

#Read in the gas species data and select species based on the element list
gas_data = add_functions.readGasSpeciesFile(gas_species_file)
select_gas_data = add_functions.selectSpeciesByElements(gas_data, element_list)

#Read in the condensate species data and select species based on the element list
cond_data = add_functions.readCondSpeciesFile(cond_species_file)
select_cond_data = add_functions.selectSpeciesByElements(cond_data, element_list)


#Write out the new files with the selected species
add_functions.writeElementAbundancesFile(new_element_abundance_file, select_element_abundances)
add_functions.writeSpeciesFile(new_gas_species_file, select_gas_data)
add_functions.writeSpeciesFile(new_cond_species_file, select_cond_data)


print("The custom network includes:")
print(len(select_element_abundances), "elements")
print(len(select_gas_data), "gas species")
print(len(select_cond_data), "condensate species")

