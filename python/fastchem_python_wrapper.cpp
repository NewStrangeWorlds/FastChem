
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
        .def_readwrite("pressure", &fastchem::FastChemInput::pressure);

    py::class_<fastchem::FastChemOutput>(m, "FastChemOutput")
        .def(py::init<>())
        .def_readwrite("number_densities", &fastchem::FastChemOutput::number_densities)
        .def_readwrite("total_element_density", &fastchem::FastChemOutput::total_element_density)
        .def_readwrite("mean_molecular_weight", &fastchem::FastChemOutput::mean_molecular_weight)
        .def_readwrite("element_conserved", &fastchem::FastChemOutput::element_conserved)
        .def_readwrite("fastchem_flag", &fastchem::FastChemOutput::fastchem_flag)
        .def_readwrite("nb_chemistry_iterations", &fastchem::FastChemOutput::nb_chemistry_iterations); 

    py::class_<fastchem::FastChem<long double>>(m, "FastChem")
        .def(py::init<const std::string &, const unsigned int>())
        .def(py::init<const std::string &, const std::string &, const unsigned int>())

        .def("calcDensities", &fastchem::FastChem<long double>::calcDensities)
        
        .def("getSpeciesNumber", &fastchem::FastChem<long double>::getSpeciesNumber)
        .def("getElementNumber", &fastchem::FastChem<long double>::getElementNumber)
        .def("getMoleculeNumber", &fastchem::FastChem<long double>::getMoleculeNumber)
        .def("getSpeciesIndex", &fastchem::FastChem<long double>::getSpeciesIndex)
        .def("getSpeciesName", &fastchem::FastChem<long double>::getSpeciesName)
        .def("getSpeciesSymbol", &fastchem::FastChem<long double>::getSpeciesSymbol)

        .def("getSpeciesMolecularWeight", &fastchem::FastChem<long double>::getSpeciesMolecularWeight)
        
        .def("setVerboseLevel", &fastchem::FastChem<long double>::setVerboseLevel)
        
        .def("getElementAbundances", &fastchem::FastChem<long double>::getElementAbundances)
        .def("setElementAbundances", &fastchem::FastChem<long double>::setElementAbundances)

        .def("setMaxChemistryIter", &fastchem::FastChem<long double>::setMaxChemistryIter)
        .def("setMaxNewtonIter", &fastchem::FastChem<long double>::setMaxNewtonIter)
        .def("setMaxBisectionIter", &fastchem::FastChem<long double>::setMaxBisectionIter)
        .def("setMaxNelderMeadIter", &fastchem::FastChem<long double>::setMaxNelderMeadIter)
        .def("setChemistryAccuracy", &fastchem::FastChem<long double>::setChemistryAccuracy)
        .def("setNewtonAccuracy", &fastchem::FastChem<long double>::setNewtonAccuracy)
        .def("useScalingFactor", &fastchem::FastChem<long double>::useScalingFactor);

        m.attr("FASTCHEM_UNKNOWN_SPECIES") = py::cast(fastchem::FASTCHEM_UNKNOWN_SPECIES);
        m.attr("FASTCHEM_SUCCESS") = py::cast(fastchem::FASTCHEM_SUCCESS);
        m.attr("FASTCHEM_NO_CONVERGENCE") = py::cast(fastchem::FASTCHEM_NO_CONVERGENCE);
        m.attr("FASTCHEM_INITIALIZATION_FAILED") = py::cast(fastchem::FASTCHEM_INITIALIZATION_FAILED);
        m.attr("FASTCHEM_IS_BUSY") = py::cast(fastchem::FASTCHEM_IS_BUSY);
        m.attr("FASTCHEM_WRONG_INPUT_VALUES") = py::cast(fastchem::FASTCHEM_WRONG_INPUT_VALUES);
        m.attr("FASTCHEM_MSG") = py::cast(fastchem::FASTCHEM_MSG);
}


