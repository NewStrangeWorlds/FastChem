
#ifndef _fastchem_h
#define _fastchem_h

#include <vector>
#include <iostream>
#include <string>


namespace fastchem {


//FastChem constants
const unsigned int FASTCHEM_UNKNOWN_SPECIES = 9999999;

const unsigned int FASTCHEM_SUCCESS = 0;
const unsigned int FASTCHEM_NO_PRESSURE_CONVERGENCE = 1;
const unsigned int FASTCHEM_NO_FASTCHEM_CONVERGENCE = 2;
const unsigned int FASTCHEM_NO_CONVERGENCE = 3;
const unsigned int FASTCHEM_INITIALIZATION_FAILED = 4;


//Physical constants
const double CONST_K = 1.3806504e-16;    //Boltzmann's constant in erg K-1



//Class for the basic chemical elements
//Those are read in from the elemental abundance file
//and are later cross-correlated with the elements in the network
template <class double_type>
struct ChemicalElement
{
  std::string symbol;
  std::string name;

  double_type atomic_weight = 0.0;
  double_type abundance = 0.0;
};



//Parent class of all species
template <class double_type>
struct ChemicalSpecies
{
  std::string symbol;
  std::string name;

  double_type molecular_weight = 0.0;
  int charge = 0;

  double_type abundance = 0.0;
  std::vector<double_type> number_density;
};



//Class for the molecules
template <class double_type>
struct Molecule : public ChemicalSpecies<double_type>
{
  std::vector<unsigned int> element_indices;
  std::vector<int> stoichometric_vector;

  std::vector<double_type> mass_action_coeff;
  std::vector<double_type> mass_action_constant;

  double_type abundance_scaled = 0.0;

  double_type sigma = 0.0;
  std::vector<double_type> sum;

  void calcMassActionConstant(const double temperature, const unsigned int grid_index);
};



//Class for the elements
template <class double_type>
struct Element : public ChemicalSpecies<double_type>
{
  unsigned int element_index;
  unsigned int index;

  unsigned int solver_order;

  std::vector<unsigned int> molecule_list;         //contains the list of molecule indices, the element is part of
  std::vector<unsigned int> element_conserved;     //check if element is conserved during calculation, for electrons this is charge conservation
};



//FastChem class
template <class double_type>
class FastChem{
  public:
    FastChem(const std::string& model_parameter_file, const unsigned int verbose_level_init);

    FastChem(const FastChem &obj);

    //function calls to calculate number densities, first one uses single temperatures and pressures,
    //the second takes temperature&pressure vectors as input
    unsigned int calcDensities(const double temperature, const double pressure,
                               std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out);

    unsigned int calcDensities(const std::vector<double>& temperature, const std::vector<double>& pressure,
                               std::vector < std::vector<double> >& density_out,
                               std::vector<double>& h_density_out, std::vector<double>& mean_molecular_weight_out);

    //same as above but with diagnostic output
    unsigned int calcDensities(const double temperature, const double pressure,
                               std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out,
                               std::vector<unsigned int>& element_conserved_out,
                               unsigned int& nb_iterations_out, unsigned int& nb_chemistry_iterations_out);

    unsigned int calcDensities(const std::vector<double>& temperature, const std::vector<double>& pressure,
                               std::vector < std::vector<double> >& density_out,
                               std::vector<double>& h_density_out, std::vector<double>& mean_molecular_weight_out,
                               std::vector< std::vector<unsigned int> >& element_conserved_out,
                               std::vector<unsigned int>& fastchem_flags,
                               std::vector<unsigned int>& nb_iterations_out, std::vector<unsigned int>& nb_chemistry_iterations_out);

    //on special request: a version which uses p_<H> instead of p_gas
    unsigned int calcDensities(const std::vector<double>& temperature, const std::vector<double>& hydrogen_pressure,
                               std::vector < std::vector<double> >& density_out,
                               std::vector<double>& mean_molecular_weight_out,
                               std::vector< std::vector<unsigned int> >& element_conserved_out,
                               std::vector<unsigned int>& fastchem_flags,
                               std::vector<unsigned int>& nb_chemistry_iterations_out);

    //public query functions
    std::string getSpeciesName(const unsigned int species_index);
    std::string getSpeciesSymbol(const unsigned int species_index);
    unsigned int getSpeciesNumber() {return nb_species;}
    unsigned int getSpeciesIndex(const std::string symbol);

    std::string getElementName(const unsigned int species_index);
    std::string getElementSymbol(const unsigned int species_index);
    unsigned int getElementNumber() {return nb_elements;}

    double getSpeciesMolecularWeight(const unsigned int species_index);


    //functions to set internal variables during runtime
    //they will override any read-in values
    void setVerboseLevel(const unsigned int level) { if (level > 4) verbose_level = 4; else verbose_level = level;}
    void setMaxChemistryIter(const unsigned int nb_steps) {nb_max_fastchem_iter = nb_steps;}
    void setMaxPressureIter(const unsigned int nb_steps) {nb_max_pressure_iter = nb_steps;}
    void setMaxNewtonIter(const unsigned int nb_steps) {nb_max_newton_iter = nb_steps;}
    void setMaxBisectionIter(const unsigned int nb_steps) {nb_max_bisection_iter = nb_steps;}
    void setMaxNelderMeadIter(const unsigned int nb_steps) {nb_max_neldermead_iter = nb_steps;}

    void setChemistryAccuracy(const double chem_accuracy) {accuracy = chem_accuracy;}
    void setPressureAccuracy(const double pressure_accuracy) {accuracy_delta = pressure_accuracy;}
    void setNewtonAccuracy(const double newton_accuracy) {newton_err = newton_accuracy;}

    void useScalingFactor(const bool use_switch) {use_scaling_factor = use_switch;}


  private:
    unsigned int nb_chemical_elements = 0;
    unsigned int nb_species = 0;
    unsigned int nb_molecules = 0;
    unsigned int nb_elements = 0;

    unsigned int e_ = FASTCHEM_UNKNOWN_SPECIES; //electron element index

    unsigned int nb_max_fastchem_iter = 300;
    unsigned int nb_max_pressure_iter = 100;
    unsigned int nb_max_bisection_iter = 3000;
    unsigned int nb_max_newton_iter = 20000;
    unsigned int nb_max_neldermead_iter = 3000;

    double_type accuracy = 1e-4;
    double_type accuracy_delta = 1e-4;
    double_type newton_err = 1e-4;


    double_type element_density_minlimit = 1e-300; //smallest allowed particle number densities
    double_type molecule_density_minlimit = 1e-300;

    unsigned int verbose_level = 1;
    bool use_scaling_factor = false;

    bool is_initialized = false;


    std::string chemical_element_file;
    std::string species_data_file;
    std::string element_abundances_file;


    std::vector< ChemicalElement<double_type> > chemical_elements;

    std::vector< ChemicalSpecies<double_type>* > species;
    std::vector< Element<double_type> > elements;
    std::vector< Molecule<double_type> > molecules;

    std::vector<unsigned int> element_calculation_order;


    //Initialisation functions
    void init();

    bool readParameterFile(const std::string& model_parameter_file);

    bool readElementList();
    bool readElementAbundances();
    void setElementAbundance(const std::string symbol, const double abundance);
    bool readSpeciesData();
    void addMolecule(const std::string name, const std::string symbol,
                     const std::vector<std::string> species_elements, const std::vector<int> stoichometric_coeff,
                     const std::vector<double_type> mass_action_coeff, const int charge);
    void addAtom(std::string symbol);

    unsigned int determineSolverOrder(const Element<double_type>& species);
    void determineSolverOrder();
    void determineElementCalculationOrder();
    double_type setInitialHDensity(const double_type total_density, const unsigned int grid_index);


    //Query functions
    unsigned int getChemicalElementIndex(const std::string symbol);
    unsigned int getElementIndex(const std::string symbol);
    unsigned int getMoleculeIndex(const std::string symbol);


    //Functions for the calculations of the number densities
    unsigned int calcDensity(const double temperature, const double pressure, const unsigned int grid_index,
                             std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out,
                             std::vector<unsigned int>& element_conserved_out,
                             unsigned int& nb_pressure_iterations_out, unsigned int& nb_chemistry_iterations_out);

    //On special request: a version using the hydrogen pressure p_<H> instead of p_gas
    unsigned int calcDensity(const double temperature, const double hydrogen_pressure, const unsigned int grid_index,
                             std::vector<double>& density_n_out, double& mean_molecular_weight_out,
                             std::vector<unsigned int>& element_conserved_out,
                             unsigned int& nb_chemistry_iterations_out);


    bool solveFastchem(const double temperature_gas, const double_type h_density, const unsigned int grid_index, unsigned int& nb_iterations);

    void calculateElementDensities(Element<double_type>& species, const double_type h_density, const double_type number_density_min,
                                   const unsigned int grid_index, bool use_backup_solver);
    void calculateMoleculeDensities(const double_type& h_density, const unsigned int grid_index);
    void calculateMinorSpeciesDensities(std::vector<double_type>& number_density_min, const unsigned int grid_index);
    void calculateElectronDensities(const double_type& old_number_density, const double_type& h_density, const unsigned int grid_index);

    bool calcTotalHydrogenDensity(const double temperature_gas, const double pressure, const unsigned int grid_index,
                                  double_type& h_density, double_type& density_iteration_lambda, double_type& density_iteration_error);
    bool calcTotalHydrogenDensityAlt(const double temperature_gas, const double pressure, const unsigned int grid_index,
                                     double_type& h_density, double_type& density_iteration_lambda, double_type& density_iteration_error);
    double calcMeanMolecularWeight(const double total_density, unsigned int grid_index);

    //Check the solution for consistency
    void checkN(Element<double_type>& species, const double_type h_density, const unsigned int grid_index);
    void checkN(Molecule<double_type>& species, const double_type h_density, const unsigned int grid_index);
    bool checkChargeConservation(const unsigned int grid_index);
    bool checkElementConservation(Element<double_type>& species, const double_type h_density, const unsigned int grid_index);

    //FastChem solvers
    void intertSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index);
    void linSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index);
    void quadSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index);
    void newtSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index);
    void newtonSolveAlt(Element<double_type>& species, const double_type h_density, const unsigned int grid_index);

    double_type solverScalingFactor(Element<double_type>& species, const double_type number_density_min, const double_type h_density, const unsigned int grid_index);

    void backupSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index);
    double_type nelderMeadFunction(Element<double_type>& species, const double_type x, const double_type h_density, const unsigned int grid_index);
    double_type bisectionFunction(Element<double_type>& species, const double_type x, const double_type h_density, const unsigned int grid_index);
    bool nelderMeadSimplexSolve(Element<double_type>& species, const double_type initial_solution, const double h_density, const unsigned int grid_index);
    bool bisectionSolve(Element<double_type>& species, const double h_density, const unsigned int grid_index);

};



}

#endif
