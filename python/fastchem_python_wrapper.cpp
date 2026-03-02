
#ifdef _SETUP_PY
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#else
#include "../_deps/pybind11-src/include/pybind11/pybind11.h"
#include "../_deps/pybind11-src/include/pybind11/stl.h"
#endif

#include "../fastchem_src/fastchem.h"
#include "../fastchem_src/input_output_struct.h"
#include "../fastchem_src/fastchem_constants.h"


namespace py = pybind11;

PYBIND11_MODULE(pyfastchem, m) {
    py::class_<fastchem::FastChemInput>(m, "FastChemInput")
        .def(py::init<>())
        .def_readwrite("temperature", &fastchem::FastChemInput::temperature)
        .def_readwrite("pressure", &fastchem::FastChemInput::pressure)
        .def_readwrite("equilibrium_condensation", &fastchem::FastChemInput::equilibrium_condensation)
        .def_readwrite("rainout_condensation", &fastchem::FastChemInput::rainout_condensation);

    py::class_<fastchem::FastChemOutput>(m, "FastChemOutput")
        .def(py::init<>())
        .def_readwrite("number_densities", &fastchem::FastChemOutput::number_densities)
        .def_readwrite("total_element_density", &fastchem::FastChemOutput::total_element_density)
        .def_readwrite("mean_molecular_weight", &fastchem::FastChemOutput::mean_molecular_weight)
        .def_readwrite("number_densities_cond", &fastchem::FastChemOutput::number_densities_cond)
        .def_readwrite("element_cond_degree", &fastchem::FastChemOutput::element_cond_degree)
        .def_readwrite("element_conserved", &fastchem::FastChemOutput::element_conserved)
        .def_readwrite("fastchem_flag", &fastchem::FastChemOutput::fastchem_flag)
        .def_readwrite("nb_chemistry_iterations", &fastchem::FastChemOutput::nb_chemistry_iterations)
        .def_readwrite("nb_cond_iterations", &fastchem::FastChemOutput::nb_cond_iterations)
        .def_readwrite("nb_iterations", &fastchem::FastChemOutput::nb_iterations);

    py::class_<fastchem::FastChem>(m, "FastChem")
        .def(py::init<const std::string &, const unsigned int>())
        .def(py::init<const std::string &, const std::string &, const unsigned int>())
        .def(py::init<const std::string &, const std::string &, const std::string &, const unsigned int>())

        .def("calcDensities", &fastchem::FastChem::calcDensities)
        
        .def("getGasSpeciesNumber", &fastchem::FastChem::getGasSpeciesNumber)
        .def("getElementNumber", &fastchem::FastChem::getElementNumber)
        .def("getMoleculeNumber", &fastchem::FastChem::getMoleculeNumber)
        .def("getCondSpeciesNumber", &fastchem::FastChem::getCondSpeciesNumber)
        
        .def("getGasSpeciesIndex", &fastchem::FastChem::getGasSpeciesIndex)
        .def("getGasSpeciesName", &fastchem::FastChem::getGasSpeciesName)
        .def("getGasSpeciesSymbol", &fastchem::FastChem::getGasSpeciesSymbol)

        .def("getElementIndex", &fastchem::FastChem::getElementIndex)
        .def("getElementName", &fastchem::FastChem::getElementName)
        .def("getElementSymbol", &fastchem::FastChem::getElementSymbol)
        
        .def("getCondSpeciesIndex", &fastchem::FastChem::getCondSpeciesIndex)
        .def("getCondSpeciesName", &fastchem::FastChem::getCondSpeciesName)
        .def("getCondSpeciesSymbol", &fastchem::FastChem::getCondSpeciesSymbol)

        .def("getGasSpeciesWeight", &fastchem::FastChem::getGasSpeciesWeight)
        .def("getElementWeight", &fastchem::FastChem::getElementWeight)
        .def("getCondSpeciesWeight", &fastchem::FastChem::getCondSpeciesWeight)

        .def("getGasSpeciesStoichiometry", &fastchem::FastChem::getGasSpeciesStoichiometry)
        .def("getCondSpeciesStoichiometry", &fastchem::FastChem::getCondSpeciesStoichiometry)

        .def("convertToHillNotation", &fastchem::FastChem::convertToHillNotation)

        .def("setVerboseLevel", &fastchem::FastChem::setVerboseLevel)
        
        .def("getElementAbundance", &fastchem::FastChem::getElementAbundance)
        .def("getElementAbundances", &fastchem::FastChem::getElementAbundances)
        .def("setElementAbundances", &fastchem::FastChem::setElementAbundances)

        .def("setParameter", static_cast<void (fastchem::FastChem::*)(const std::string&, const unsigned int)>(&fastchem::FastChem::setParameter))
        .def("setParameter", static_cast<void (fastchem::FastChem::*)(const std::string&, const bool)>(&fastchem::FastChem::setParameter))
        .def("setParameter", static_cast<void (fastchem::FastChem::*)(const std::string&, const double)>(&fastchem::FastChem::setParameter));

        m.attr("FASTCHEM_UNKNOWN_SPECIES") = py::cast(fastchem::FASTCHEM_UNKNOWN_SPECIES);
        m.attr("FASTCHEM_SUCCESS") = py::cast(fastchem::FASTCHEM_SUCCESS);
        m.attr("FASTCHEM_NO_CONVERGENCE") = py::cast(fastchem::FASTCHEM_NO_CONVERGENCE);
        m.attr("FASTCHEM_INITIALIZATION_FAILED") = py::cast(fastchem::FASTCHEM_INITIALIZATION_FAILED);
        m.attr("FASTCHEM_IS_BUSY") = py::cast(fastchem::FASTCHEM_IS_BUSY);
        m.attr("FASTCHEM_WRONG_INPUT_VALUES") = py::cast(fastchem::FASTCHEM_WRONG_INPUT_VALUES);
        m.attr("FASTCHEM_PHASE_RULE_VIOLATION") = py::cast(fastchem::FASTCHEM_PHASE_RULE_VIOLATION);
        m.attr("FASTCHEM_MSG") = py::cast(fastchem::FASTCHEM_MSG);
}


