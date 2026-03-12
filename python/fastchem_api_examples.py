
import pyfastchem


#As a start, we have to create a FastChem object
fastchem = pyfastchem.FastChem(
  '../input/element_abundances/asplund_2021_extended.dat',
  '../input/logK/logK_extended.dat',
  '../input/logK/logK_condensates_extended.dat',
  0)


#First, we want to find the total number of the different species classes
nb_elements = fastchem.getElementNumber()
nb_molecules = fastchem.getMoleculeNumber()
nb_gas_species = fastchem.getGasSpeciesNumber()
nb_cond_species = fastchem.getCondSpeciesNumber()

print('FastChem includes:')
print('  ',nb_elements,'elements')
print('  ',nb_molecules,' molecules/ions')
print('   total number of gas species:', nb_gas_species)
print('  ',nb_cond_species, 'condensates')

print('\n')


#Now we want to find out how many species are available for each element
#For that, we loop through all elements and check the stoichiometry of each species
for element in range(nb_elements):
  nb_species = 0
  
  #look through the gas species
  for j in range(nb_gas_species):
    stoichiometry = fastchem.getGasSpeciesStoichiometry(j)
    if stoichiometry[element] != 0:
      nb_species += 1
  
  #look through the condensates
  for j in range(nb_cond_species):
    stoichiometry = fastchem.getCondSpeciesStoichiometry(j)
    if stoichiometry[element] != 0:
      nb_species += 1

  print(fastchem.getElementSymbol(element),'has', nb_species,'species')

print('\n')


#Next, we want to find all species that contain two carbon atoms
nb_c2_species = 0
c2_species = []

#for this, we need to know the index of carbon first
index_c = fastchem.getElementIndex('C')

for i in range(nb_gas_species):
  stoichiometry = fastchem.getGasSpeciesStoichiometry(i)
  if stoichiometry[index_c] == 2:
    nb_c2_species += 1
    c2_species.append(fastchem.getGasSpeciesSymbol(i))

for i in range(nb_cond_species):
  stoichiometry = fastchem.getCondSpeciesStoichiometry(i)
  if stoichiometry[index_c] == 2:
    nb_c2_species += 1
    c2_species.append(fastchem.getCondSpeciesSymbol(i))

print('Total number of species with 2 C atoms:', nb_c2_species)
print('These are:', end=' ')
for species in c2_species:
  print(species, end=' ')

print('\n')


#Now, we want to find the species with the most carbon atoms 
#in the gas phase and in the condensates
species_index_g = 0
nb_c_atoms_g = 0

for i in range(nb_gas_species):
  stoichiometry = fastchem.getGasSpeciesStoichiometry(i)
  if stoichiometry[index_c] > nb_c_atoms_g:
    nb_c_atoms_g = stoichiometry[index_c]
    species_index_g = i

#find the condensate species with the most carbon atoms
species_index_c = 0
nb_c_atoms_c = 0

for i in range(nb_cond_species):
  stoichiometry = fastchem.getCondSpeciesStoichiometry(i)
  if stoichiometry[index_c] > nb_c_atoms_c:
    nb_c_atoms_c = stoichiometry[index_c]
    species_index_c = i

print("Gas-phase species with the most C atoms:", 
      fastchem.getGasSpeciesSymbol(species_index_g))
print("Condensate species with the most C atoms:", 
      fastchem.getCondSpeciesSymbol(species_index_c))

print('\n')


#What is the heaviest species in the gas phase and in the condensates?
#Find the heaviest gas-phase species
species_index_g = 0
max_mass_g = 0.0

for i in range(nb_gas_species):
  mass = fastchem.getGasSpeciesWeight(i)
  
  if mass > max_mass_g:
    max_mass_g = mass
    species_index_g = i

print("Heaviest gas-phase species:", 
      fastchem.getGasSpeciesSymbol(species_index_g), 
      "with molecular weight", 
      max_mass_g, 'amu')


#Find the heaviest condensate species
species_index_c = 0
max_mass_c = 0.0

for i in range(nb_cond_species):
  mass = fastchem.getCondSpeciesWeight(i)
 
  if mass > max_mass_c:
    max_mass_c = mass
    species_index_c = i

print("Heaviest condensate species:", 
      fastchem.getCondSpeciesSymbol(species_index_c), 
      "with molecular weight", 
      max_mass_c, 'amu')