
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

    py::class_<fastchem::FastChem<long double>>(m, "FastChem")
        .def(py::init<const std::string &, const unsigned int>())
        .def(py::init<const std::string &, const std::string &, const unsigned int>())
        .def(py::init<const std::string &, const std::string &, const std::string &, const unsigned int>())

        .def("calcDensities", &fastchem::FastChem<long double>::calcDensities)
        
        .def("getGasSpeciesNumber", &fastchem::FastChem<long double>::getGasSpeciesNumber)
        .def("getElementNumber", &fastchem::FastChem<long double>::getElementNumber)
        .def("getMoleculeNumber", &fastchem::FastChem<long double>::getMoleculeNumber)
        .def("getCondSpeciesNumber", &fastchem::FastChem<long double>::getCondSpeciesNumber)
        
        .def("getGasSpeciesIndex", &fastchem::FastChem<long double>::getGasSpeciesIndex)
        .def("getGasSpeciesName", &fastchem::FastChem<long double>::getGasSpeciesName)
        .def("getGasSpeciesSymbol", &fastchem::FastChem<long double>::getGasSpeciesSymbol)

        .def("getElementIndex", &fastchem::FastChem<long double>::getElementIndex)
        .def("getElementName", &fastchem::FastChem<long double>::getElementName)
        .def("getElementSymbol", &fastchem::FastChem<long double>::getElementSymbol)
        
        .def("getCondSpeciesIndex", &fastchem::FastChem<long double>::getCondSpeciesIndex)
        .def("getCondSpeciesName", &fastchem::FastChem<long double>::getCondSpeciesName)
        .def("getCondSpeciesSymbol", &fastchem::FastChem<long double>::getCondSpeciesSymbol)

        .def("getGasSpeciesWeight", &fastchem::FastChem<long double>::getGasSpeciesWeight)
        .def("getElementWeight", &fastchem::FastChem<long double>::getElementWeight)
        .def("getCondSpeciesWeight", &fastchem::FastChem<long double>::getCondSpeciesWeight)

        .def("setVerboseLevel", &fastchem::FastChem<long double>::setVerboseLevel)
        
        .def("getElementAbundance", &fastchem::FastChem<long double>::getElementAbundance)
        .def("getElementAbundances", &fastchem::FastChem<long double>::getElementAbundances)
        .def("setElementAbundances", &fastchem::FastChem<long double>::setElementAbundances)

        .def("setParameter", static_cast<void (fastchem::FastChem<long double>::*)(const std::string&, const unsigned int)>(&fastchem::FastChem<long double>::setParameter))
        .def("setParameter", static_cast<void (fastchem::FastChem<long double>::*)(const std::string&, const bool)>(&fastchem::FastChem<long double>::setParameter))
        .def("setParameter", static_cast<void (fastchem::FastChem<long double>::*)(const std::string&, const long double)>(&fastchem::FastChem<long double>::setParameter));

        m.attr("FASTCHEM_UNKNOWN_SPECIES") = py::cast(fastchem::FASTCHEM_UNKNOWN_SPECIES);
        m.attr("FASTCHEM_SUCCESS") = py::cast(fastchem::FASTCHEM_SUCCESS);
        m.attr("FASTCHEM_NO_CONVERGENCE") = py::cast(fastchem::FASTCHEM_NO_CONVERGENCE);
        m.attr("FASTCHEM_INITIALIZATION_FAILED") = py::cast(fastchem::FASTCHEM_INITIALIZATION_FAILED);
        m.attr("FASTCHEM_IS_BUSY") = py::cast(fastchem::FASTCHEM_IS_BUSY);
        m.attr("FASTCHEM_WRONG_INPUT_VALUES") = py::cast(fastchem::FASTCHEM_WRONG_INPUT_VALUES);
        m.attr("FASTCHEM_PHASE_RULE_VIOLATION") = py::cast(fastchem::FASTCHEM_PHASE_RULE_VIOLATION);
        m.attr("FASTCHEM_MSG") = py::cast(fastchem::FASTCHEM_MSG);
}


